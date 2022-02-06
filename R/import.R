#' combineTCR
#'
#' @param df List of filtered contig annotations from 10x Genomics.
#' @param df_clono List of clonotypes metadata from 10x Genomics.
#' @param samples The labels of samples (required).
#' @param ID The additional sample labeling (optional).
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
#'
combineTCR <- function (df, df_clono, samples = NULL, ID = NULL, cells = "T-AB", removeNA = FALSE,
                               removeMulti = FALSE, filterMulti = FALSE)
{
  source_url("https://raw.githubusercontent.com/ncborcherding/scRepertoire/master/R/utils.R")
  df <- checkList(df)
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
  for (i in seq_along(df)) {
    df[[i]] <- subset(df[[i]], chain != "Multi")
    df[[i]] <- subset(df[[i]], chain %in% c(chain1, chain2))
    df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
    if (nrow(df[[i]]) == 0) {
      stop("There are 0 contigs after internal filtering -\n check the contig list to see if any issues exist \n        for productive chains",
           call. = FALSE)
    }
    df[[i]] <- subset(df[[i]], cdr3 != "None")
    df[[i]]$sample <- samples[i]
    df[[i]]$ID <- ID[i]
    if (filterMulti == TRUE) {
      df[[i]] <- filteringMulti(df[[i]])
    }
  }
  if (!is.null(samples)) {
    out <- modifyBarcodes(df, samples, ID)
  }
  else {
    out <- df
  }
  for (i in seq_along(out)) {
    data2 <- out[[i]]
    data_clonotype <- df_clono[[i]]
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
    if (!is.null(sample) & !is.null(ID)) {
      data3 <- data3[, c("barcode", "sample", "ID", tcr1_lines,
                         tcr2_lines, CT_lines, "raw_clonotype_id")]
    }
    else if (!is.null(sample) & is.null(ID)) {
      data3 <- data3[, c("barcode", "sample", tcr1_lines,
                         tcr2_lines, CT_lines, "raw_clonotype_id")]
    }
    data3 <- as_tibble(data3) %>% left_join(as_tibble(data_clonotype), by = c("raw_clonotype_id" = "clonotype_id"))
    data3 <- modifyClonotypes(data3, "sample")
    data3 <- as.data.frame(data3)
    final[[i]] <- data3
  }
  names <- NULL
  for (i in seq_along(samples)) {
    if (!is.null(sample) & !is.null(ID)) {
      c <- paste(samples[i], "_", ID[i], sep = "")
    }
    else if (!is.null(sample) & is.null(ID)) {
      c <- paste(samples[i], sep = "")
    }
    names <- c(names, c)
  }
  names(final) <- names
  for (i in seq_along(final)) {
    final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),
    ]
    final[[i]] <- final[[i]][rowSums(is.na(final[[i]])) <
                               10, ]
  }
  if (removeNA == TRUE) {
    final <- removingNA(final)
  }
  if (removeMulti == TRUE) {
    final <- removingMulti(final)
  }
  return(final)
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
