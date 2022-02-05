combineTCR_custom <- function (df, samples = NULL, ID = NULL, cells = "T-AB", removeNA = FALSE,
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
    df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE",
                                                 "True", "true"))
    if (nrow(df[[i]]) == 0) {
      stop("There are 0 contigs after internal filtering -\n        check the contig list to see if any issues exist \n        for productive chains",
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
    data2 <- makeGenes(cellType, data2, chain1, chain2)
    unique_df <- unique(data2$barcode)
    Con.df <- data.frame(matrix(NA, length(unique_df), 7))
    colnames(Con.df) <- c("barcode", tcr1_lines, tcr2_lines)
    Con.df$barcode <- unique_df
    Con.df <- parseTCR_custom(Con.df, unique_df, data2)
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

checkList <- function(df) {
  df <- if(is(df)[1] != "list") list(df) else df
  return(df)
}

parseTCR_custom <- function(Con.df, unique_df, data2) {

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
  return(Con.df)}
