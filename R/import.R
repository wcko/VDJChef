#' combineTCR
#'
#' @param df Filtered contig info (10X Genomics VDJ output: filtered_contig_annotations.csv)
#' @param clonotype_info Clonotype metadata info(10X Genomics VDJ output: clonotypes.csv)
#' @param cells The type of T cell - T cell-AB or T cell-GD. Only 1 T cell type can
#' be called at once.
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 2
#' corresponding chains with the highest expression for a single barcode.

#' @import dplyr
#' @importFrom devtools source_url
#' @export
#' @return List of clonotypes for individual cell barcodes
#'
#' Example dataset: https://www.10xgenomics.com/resources/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0
#'
#' https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0
#'
#' seu.data <- Read10X_h5(filename = "inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_count_filtered_feature_bc_matrix.h5")
#' seu <- CreateSeuratObject(counts=seu.data)
#' contigs <- read.csv("inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_filtered_contig_annotations.csv")
#' clonotypes <- read.csv("inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_clonotypes.csv")
#' [to do] set this up so data("contig_list") returns the data built into this package... eventually do it for our own inflamm diseases data
#' combineTCR(list(contigs),list(clonotypes),samples="A")
#' meta <- seu@meta.data
#' meta$barcode <- rownames(meta)
#' test$A
#' left_join(meta,test$A,by="barcode")
#'
contigs <- read.csv("inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_filtered_contig_annotations.csv")
clonotype_info <- read.csv("inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_clonotypes.csv")


combineTCR <- function (df, clonotype_info, cells = "T-AB", removeNA = FALSE,
                        removeMulti = FALSE, filterMulti = FALSE)
{
  # source_url("https://raw.githubusercontent.com/ncborcherding/scRepertoire/master/R/utils.R")
  out <- NULL
  final <- NULL
  tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
  tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
  data1_lines <- c("TCR1", "cdr3", "cdr3_nt")
  data2_lines <- c("TCR2", "cdr3", "cdr3_nt")
  CT_lines <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cellType")

  chain1 <- cellT(cells)[[1]]
  chain2 <- cellT(cells)[[2]]
  cellType <- cellT(cells)[[3]]

  df <- subset(df, chain != "Multi")
  df <- subset(df, chain %in% c(chain1, chain2))
  df <- subset(df, productive %in% c(TRUE, "TRUE", "True", "true"))
  if (nrow(df) == 0) {
    stop("There are 0 contigs after internal filtering -\n check the contig list to see if any issues exist \n        for productive chains",
         call. = FALSE)
  }
  df <- subset(df, cdr3 != "None")
  if (filterMulti == TRUE) {
    out <- filteringMulti(df)
  }
  else {
    out <- df
  }
  data2 <- out
  data_clonotype <- clonotype_info
  data2 <- makeGenesNew(cellType, data2, chain1, chain2)
  unique_df <- unique(data2$barcode)
  Con.df <- data.frame(matrix(NA, length(unique_df), 7))
  colnames(Con.df) <- c("barcode", tcr1_lines, tcr2_lines)
  Con.df$barcode <- unique_df
  Con.df <- parseTCR(Con.df, unique_df, data2)
  Con.df <- assignCT(cellType, Con.df)
  Con.df$cellType <- cells
  Con.df[Con.df == "NA_NA" | Con.df == "NA_NA_NA_NA"] <- NA
  data3 <- merge(data2[, -which(names(data2) %in% c("TCR1","TCR2"))], Con.df, by = "barcode")
  data3 <- data3[, c("barcode", tcr1_lines,
                     tcr2_lines, CT_lines, "raw_clonotype_id")]
  data3 <- as_tibble(data3) %>% left_join(as_tibble(data_clonotype), by = c("raw_clonotype_id" = "clonotype_id"))
  data3 <- as.data.frame(data3)
  final <- data3

  if (removeNA == TRUE) {
    final <- na.omit(final)
  }
  if (removeMulti == TRUE) {
    final <- filter(final, !grepl(";",CTnt))
  }
  return(final)
}

#' Adding clonotype information to a Seurat object or ExpressionSet object
#'
#' @examples
#' #Get combined contigs
#' combined <- combineTCR(contig_list,clonotype_info)
#'
#' sc <- Read10X_h5(filename = "inst/extdata/sc5p_v2_hs_T_1k_multi_5gex_t_count_filtered_feature_bc_matrix.h5")
#' sc <- CreateSeuratObject(counts=sc)
#'
#' sc <- ExpressionSet object ## need to provide example data, along with clonotypes. Jake data?
#'
#' #Using combineExpression()
#' combineExpression(combined,sc)
#'
#' @param df The product of CombineTCR() or CombineBCR() or a list of
#' both c(CombineTCR(), combineBCR())
#' @param sc The seurat or SingleCellExperiment (SCE) object or ExpressionSet to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "Heavy", "Light"
#' @param group.by The column label in the combined contig object in which
#' clonotype frequency will be calculated.
#' @param proportion Whether to use the total frequency (FALSE) or the
#' proportion (TRUE) of the clonotype based on the group.by variable.
#' @param cloneTypes The bins for the grouping based on frequency
#' @param filterNA Method to subset seurat object of barcodes without
#' clonotype information
#' @param addLabel This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after
#' subseting the data.
#' @importFrom dplyr bind_rows %>% summarise
#' @importFrom  rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @export
#' @return Seurat, ExpressionSet, or SingleCellExperiment object with attached clonotype
#' information

combineExpression <- function(df, sc, cloneCall="gene+nt",
                              chain = "both", group.by="none",
                              proportion = TRUE, filterNA = FALSE,
                              cloneTypes=c(Rare = 1e-4, Small = 0.001,
                                           Medium = 0.01, Large = 0.1, Hyperexpanded = 1),
                              addLabel = FALSE) {
  options( dplyr.summarise.inform = FALSE )
  cloneTypes <- c(None = 0, cloneTypes)
  df <- checkList(df)
  cloneCall <- theCall(cloneCall)
  Con.df <- NULL
  meta <- grabMeta(sc)
  cell.names <- rownames(meta)
  if (group.by == "none") {
    for (i in seq_along(df)) {
      if (chain != "both") {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
      data <- data.frame(df[[i]], stringsAsFactors = FALSE)
      data2 <- unique(data[,c("barcode", cloneCall)])
      data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
      if (proportion == TRUE) {
        data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
          summarise(Frequency = n()/nrow(data2))
      } else {
        data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
          summarise(Frequency = n())
      }
      colnames(data2)[1] <- cloneCall
      data <- merge(data, data2, by = cloneCall, all = TRUE)
      # data <- data[,c("barcode", "CTgene", "CTnt",
      #                 "CTaa", "CTstrict", "Frequency")]
      Con.df <- rbind.data.frame(Con.df, data)
    }
  } else if (group.by != "none") {
    data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
    data2 <- na.omit(unique(data[,c("barcode", cloneCall, group.by)]))
    data2 <- data2[data2[,"barcode"] %in% cell.names, ]
    data2 <- as.data.frame(data2 %>% group_by(data2[,cloneCall],
                                              data2[,group.by]) %>% summarise(Frequency = n()))
    colnames(data2)[c(1,2)] <- c(cloneCall, group.by)
    x <- unique(data[,group.by])
    for (i in seq_along(x)) {
      sub1 <- subset(data, data[,group.by] == x[i])
      sub2 <- subset(data2, data2[,group.by] == x[i])
      merge <- merge(sub1, sub2, by=cloneCall)
      if (proportion == TRUE) {
        merge$Frequency <- merge$Frequency/length(merge$Frequency)
      }
      Con.df <- rbind.data.frame(Con.df, merge)
    }
    nsize <- Con.df %>% group_by(Con.df[,paste0(group.by, ".x")])  %>% summarise(n = n())
  }

  Con.df$cloneType <- NA
  for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <-
    paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1],
           ' < X <= ', cloneTypes[x], ')') }
  for (i in 2:length(cloneTypes)) { Con.df$cloneType <-
    ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency
           <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
  # PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt",
  #                             "CTaa", "CTstrict", "Frequency", "cloneType")])
  PreMeta <- unique(Con.df[,]) #specify columns to join
  dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
  PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
  rownames(PreMeta) <- PreMeta$barcode
  if (group.by != "none" && addLabel) {
    location <- which(colnames(PreMeta) == "Frequency")
    colnames(PreMeta)[location] <- paste0("Frequency.", group.by)
  }
  if (inherits(x=sc, what ="Seurat")) {
    if (length(which(rownames(PreMeta) %in%
                     rownames(sc[[]])))/length(rownames(sc[[]])) < 0.01) {
      warning("< 1% of barcodes match: Ensure the barcodes in
            the Seurat object match the
            barcodes in the combined immune receptor list")
    }
    col.name <- names(PreMeta) %||% colnames(PreMeta)
    sc[[col.name]] <- PreMeta
  }
  else if (inherits(x=sc, what="ExpressionSet")) {
    if (length(which(rownames(PreMeta) %in%
                     rownames(pData(sc))))/length(rownames(pData(sc))) < 0.01) {
      warning("< 1% of barcodes match: Ensure the barcodes in
            the Seurat object match the
            barcodes in the combined immune receptor list")
    }
    col.name <- names(PreMeta) %||% colnames(PreMeta)
    sc_metadata <- pData(sc) %>% left_join(PreMeta, by="barcode")
    rownames(sc_metadata) <- sc_metadata$barcode
    pData(sc) <- sc_metadata
  }
  else {
    rownames <- rownames(colData(sc))
    if (length(which(rownames(PreMeta) %in%
                     rownames))/length(rownames) < 0.01) {
      warning("< 1% of barcodes match: Ensure the barcodes
          in the SingleCellExperiment object match the
          barcodes in the combined immune receptor list from
          scRepertoire - most common issue is the addition of the
          prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR()
          functions") }
    colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
    rownames(colData(sc)) <- rownames
  }
  if (filterNA == TRUE) { sc <- filteringNA(sc) }
  return(sc)
}




#Sorting the V/D/J/C gene sequences for T and B cells
#' @importFrom stringr str_c
makeGenesNew <- function(cellType, data2, chain1, chain2) {
  if(cellType %in% c("T-AB", "T-GD")) {
    data2 <- data2 %>%
      mutate(TCR1 = ifelse(chain == chain1,
                           str_c(na.omit(v_gene),  na.omit(j_gene), na.omit(c_gene), sep = "."), NA)) %>%
      mutate(TCR2 = ifelse(chain == chain2,
                           str_c(na.omit(v_gene), na.omit(d_gene),  na.omit(j_gene),  na.omit(c_gene), sep = "."), NA))
  }
  else {
    data2 <- data2 %>%
      mutate(IGKct = ifelse(chain == "IGK",
                            str_c(na.omit(v_gene),  na.omit(j_gene), na.omit(c_gene), sep = "."), NA)) %>%
      mutate(IGLct = ifelse(chain == "IGL",
                            str_c(na.omit(v_gene),  na.omit(j_gene), na.omit(c_gene), sep = "."), NA)) %>%
      mutate(IGHct = ifelse(chain == "IGH",
                            str_c(na.omit(v_gene), na.omit(d_gene),  na.omit(j_gene),  na.omit(c_gene), sep = "."), NA))
  }
  return(data2)

}

modifyClonotypes <- function(data, samples) {
  out <- NULL
  data$raw_clonotype_id[grepl("clonotype",data$raw_clonotype_id)] <- paste0(data$raw_clonotype_id[grepl("clonotype",data$raw_clonotype_id)], "_",
                                                                            data[grepl("clonotype",data$raw_clonotype_id),][[samples]])
  return(data)
}

assignCT <- function(cellType, Con.df) {
  if (cellType %in% c("T-AB", "T-GD")) {
    Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
    Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
    Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
    Con.df$CTstrict <- paste(Con.df$TCR1, Con.df$cdr3_nt1,
                             Con.df$TCR2, Con.df$cdr3_nt2, sep="_")
  } else {
    Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
    Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
    Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_") }
  return(Con.df)
}

cellT <- function(cells) {
  if (cells == "T-AB") {
    chain1 <- "TRA"
    chain2 <- "TRB"
    cellType <- "T-AB"
  } else if (cells == "T-GD") {
    chain1 <- "TRD"
    chain2 <- "TRG"
    cellType <- "T-GD"
  } else if (cells == "B") {
    chain1 <- "IGH"
    chain2 <- "IGL"
    cellType <- "B"
  }
  return(list(chain1, chain2, cellType))
}

filteringMulti <- function(x) {
  x <- x %>%
    group_by(barcode, chain) %>%
    top_n(n = 1, wt = reads)
  table <- subset(as.data.frame(table(x$barcode)), Freq > 2)
  if (nrow(table) > 0) {
    barcodes <- as.character(unique(table$Var1))
    multichain <- NULL
    for (j in seq_along(barcodes)) {
      chain <- x[x$barcode == barcodes[j],] %>%
        group_by(barcode) %>% top_n(n = 2, wt = reads)
      multichain <- rbind(multichain, chain)
    }
    x <- subset(x, barcode %!in% barcodes)
    x <- rbind(x, multichain)
  }
  return(x)
}

parseTCR <- function(Con.df, unique_df, data2) {
  tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
  tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
  data1_lines <- c("TCR1", "cdr3", "cdr3_nt")
  data2_lines <- c("TCR2", "cdr3", "cdr3_nt")
  CT_lines <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cellType")

  for (y in seq_along(unique_df)){
    barcode.i <- Con.df$barcode[y]
    location.i <- which(barcode.i == data2$barcode)
    if (length(location.i) == 2){
      if (is.na(data2[location.i[1],c("TCR1")])) {
        Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines]
        if (is.na(data2[location.i[2],c("TCR2")])) {
          Con.df[y,tcr1_lines]<-data2[location.i[2],data1_lines]
        }
      } else {
        if(!is.na(data2[location.i[1],c("TCR1")])) {
          Con.df[y,tcr1_lines]<-data2[location.i[1],data1_lines]
        }
        if (!is.na(data2[location.i[2],c("TCR2")])) {
          Con.df[y,tcr2_lines]<-data2[location.i[2],data2_lines]
        }
      }
    } else if (length(location.i) == 3) {
      if (is.na(data2[location.i[1],c("TCR1")])) {
        Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines]
        if (is.na(data2[location.i[2],c("TCR1")])) {
          TRdf <- paste(Con.df[y, tcr2_lines],
                        data2[location.i[2], data2_lines],sep=";")
          Con.df[y,tcr2_lines] <- TRdf
          Con.df[y,tcr1_lines] <- data2[location.i[3],data1_lines]
        } else { # if the 2nd location is occupied by TRA
          Con.df[y,tcr1_lines] <- data2[location.i[2],data1_lines]
          if (is.na(data2[location.i[3],c("TCR1")])) {
            TRdf <- paste(Con.df[y, tcr2_lines],
                          data2[location.i[3], data2_lines],sep=";")
            Con.df[y,tcr2_lines] <- TRdf
          } else { # if the 3rd location is occupied by TRA
            TRdf <- paste(Con.df[y, tcr1_lines],
                          data2[location.i[3], data1_lines],sep=";")
            Con.df[y,tcr1_lines] <- TRdf }}
      } else { # if 1st location is occupied by TRA
        Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines]
        if (is.na(data2[location.i[2],c("TCR1")])) {
          if (is.na(data2[location.i[3],c("TCR1")])) { #Two TRB chains
            TRdf <- paste(data2[location.i[2], data2_lines],
                          data2[location.i[3], data2_lines],sep=";")
            Con.df[y,tcr2_lines] <- TRdf
          } else if (!is.na(data2[location.i[3],c("TCR1")])) { # if TRA is on 3rd location
            TRdf <- paste(Con.df[y, tcr1_lines],
                          data2[location.i[3],data1_lines],sep=";")
            Con.df[y,tcr1_lines] <- TRdf
            Con.df[y,tcr2_lines] <- data2[location.i[2],data2_lines]
          }
        } else { # if TRA is on 2nd location
          TRdf <- paste(Con.df[y, tcr1_lines],
                        data2[location.i[2], data1_lines],sep=";")
          Con.df[y,tcr1_lines] <- TRdf
          Con.df[y,tcr2_lines] <- data2[location.i[3],data2_lines]}}
    } else if (length(location.i) == 1) {
      chain.i <- data2$chain[location.i]
      if (chain.i == "TRA"){
        Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines]
      } else {Con.df[y,tcr2_lines] <- data2[location.i[1],data2_lines]}}}
  return(Con.df)
}

cellT <- function(cells) {
  if (cells == "T-AB") {
    chain1 <- "TRA"
    chain2 <- "TRB"
    cellType <- "T-AB"
  } else if (cells == "T-GD") {
    chain1 <- "TRD"
    chain2 <- "TRG"
    cellType <- "T-GD"
  } else if (cells == "B") {
    chain1 <- "IGH"
    chain2 <- "IGL"
    cellType <- "B"
  }
  return(list(chain1, chain2, cellType))
}

grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="ExpressionSet")){
    meta <- pData(sc)
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="SummarizedExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}
