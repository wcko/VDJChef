#' getTCR
#' @param vdj_dir Path to 10X Genomics Cellranger vdj (or multi) output of VDJ contigs
#' @param contigfile Contig info (10X Genomics VDJ output: filtered_contig_annotations.csv or all_filtered_contig_annotations.csv)
#' @param clonotypefile Clonotype metadata info(10X Genomics VDJ output: clonotypes.csv)
#' @param removeNA Remove any chain with without values
#' @param removeMulti Remove barcodes with more than 2 chains
#' @param filterMulti This option will allow for the selection of the 2
#' corresponding chains with the highest expression for a single barcode.

#' @import dplyr
#' @details Reads output directory of Cellranger VDJ (or multi) to output clonotype metadata (TCR chains, VDJ genes, reads, umis).
#' If extra info on T cell metadata is available (clonotypes.csv), it will be joined to the clonotype metadata.
#' @examples
#' getTCR(vdj_dir="inst/extdata/1kcells-multi/vdj_t/")
#' getTCR(vdj_dir="inst/extdata/1kcells-multi/vdj_t/", filteredcontigs = FALSE, removeMulti = TRUE)
#' getTCR(contigfile="inst/extdata/1kcells-multi/vdj_t/filtered_contig_annotations.csv", clonotypefile="inst/extdata/1kcells-multi/vdj_t/clonotypes.csv")
#' @export
#' @return Data frame of clonotypes for individual cell barcodes
#'
#' Example dataset: https://www.10xgenomics.com/resources/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0
#' https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0

getTCR <- function (vdj_dir, filteredcontigs = TRUE, removeNA = FALSE, removeMulti = FALSE)
{
  # determine if using filtered contigs or all contigs
  contig_file_path <- NULL
  clonotypeinfo_file_path <- NULL
  vdj_path <- vdj_dir
  vdj_files <- list.files(vdj_path)
  if ("TRUE" %in% grepl("filtered_contig_annotations.csv",vdj_files)) {
    contig_file_path <- paste0(vdj_path,vdj_files[grep("filtered_contig_annotations",vdj_files)])
  }
  else {
    contig_file_path <- paste0(vdj_path,vdj_files[grep("all_contig_annotations",vdj_files)])
    #print statement that no filtered contig file found, using all contigs
  }
  if ("TRUE" %in% grepl("clonotypes.csv",vdj_files)) {
    clonotypeinfo_file_path <- paste0(vdj_path,vdj_files[grep("clonotypes.csv",vdj_files)])
  }

  contig_file <- read.csv(contig_file_path)
  df <- contig_file %>% arrange(desc(chain),raw_consensus_id) %>%
    group_by(barcode) %>%
    summarise(chain = paste(chain,collapse = ";"), v_gene = paste(v_gene,collapse = ";"), d_gene = paste(d_gene,collapse = ";"),
              j_gene = paste(j_gene,collapse = ";"), c_gene = paste(c_gene,collapse = ";"), full_length = paste(full_length,collapse = ";"),
              productive = paste(productive, collapse = ";"), cdr3 = paste(cdr3,collapse = ";"), cdr3_nt = paste(cdr3_nt,collapse = ";"),
              reads = paste(reads,collapse = ";"), umis = paste(umis,collapse = ";"), raw_clonotype_id) %>%
    rename(clonotype_id = raw_clonotype_id) %>%
    distinct

  if (exists("clonotypeinfo_file_path")) {
    clonotypeinfo <- read.csv(clonotypeinfo_file_path)
    df <- left_join(df, clonotypeinfo, by="clonotype_id")
  }

  # or may directly pass in the contig file and / or clonotype metadata file, and NOT the vdj directory

  # if removeMulti = TRUE
  # delete rows if num chains > 2

  # if removeNA = TRUE
  # delete rows of NA
  df <- as.data.frame(df)
  rownames(df) <- df$barcode
  df <- df %>% select(-barcode)
  return(df)
}

#' addTCR
#' @param exsc SatijaLab’s Seurat Class or Bioconductor’s ExpressionSet Class
#' @param tcrmeta TCR clonotype metadata with rownames as barcodes (i.e. output of function getTCR() ). Barcodes must match that of @param exsc
#' @import dplyr
#' @details Reads output directory of Cellranger VDJ (or multi) to output clonotype metadata (TCR chains, VDJ genes, reads, umis).
#' If extra info on T cell metadata is available (clonotypes.csv), it will be joined to the clonotype metadata.
#' @examples
#' Example dataset: https://www.10xgenomics.com/resources/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0
#' https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0
#' seu.data <- Read10X_h5(filename = "inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_count_filtered_feature_bc_matrix.h5")
#' seu <- CreateSeuratObject(counts=seu.data)
#' meta <- seu@meta.data
#' meta$barcode <- rownames(meta)
#' eset <- ExpressionSet(as.matrix(seu@assays$RNA@counts))
#' pData(eset) <- seu@meta.data
#' tcrmeta <- getTCR(vdj_dir="inst/extdata/1kcells-multi/vdj_t/")
#' addTCR(seu,tcrmeta)
#' addTCR(eset,tcrmeta)
#' @export
#' @return Seurat Object or ExpressionSet with metadata updated with TCR clonotypes
#'

addTCR <- function (exsc, tcrmeta) {
  if (class(exsc) == "Seurat") {
    exsc <- AddMetaData(exsc,tcrmeta)
  }
  else if (class(exsc) == "ExpressionSet") {
    meta <- pData(exsc)
    meta$barcode <- rownames(meta)
    tcrmeta$barcode <- rownames(tcrmeta)
    meta <- left_join(meta,tcrmeta,by="barcode")
    rownames(meta) <- meta$barcode
    meta <- subset(meta, select = -barcode)
    pData(exsc) <- meta
  }
  else {
    print("Input is neither a Seurat Object or ExpressionSet")
  }
  return(exsc)
}






