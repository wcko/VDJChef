#' exportTCR
#'
#' Exports TCR metadata to csv file
#' Can be used as input to tcrdist python package
#' minimum input is v_gene, j_gene, cdr3 - this is for tcrdist python.
#' because they backcalculate v_gene aa's. Cellranger has the v_gene aa's available.
#' Ideal output is subject, epitope, count, vagene,jagene,cdr3a_aa, cdr3a_nt,vbgene,jbgene,cdr3b_aa,cdr3b_nt,clonotype_id
#'
#' @description Export TCR metadata to csv file
#' @param input Seurat Object or ExpressionSet
#' @param subject (implement - choose which subjects)
#' @param clonotype (implement - choose which clonotypes)
#' @param paired TRUE to pair alpha/beta, FALSE to separate and create two csv's (1 alpha 1 beta)
#' For now it will print out unpaired alpha beta, and paired
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#' ERICA: add in your column of epitope (dextramer column?) and subject
#'
#'
exportTCR <- function (input) {
  df <- input@meta.data %>% filter(nchar(chain)<12) %>% select(barcode,clonotype_id,v_gene,j_gene,cdr3) %>% separate_rows(v_gene,j_gene,cdr3, sep = ';')

  a_df <- df[str_detect(df$v_gene,'TRAV'), ] %>% rename(v_a_gene = v_gene, j_a_gene = j_gene, cdr3_a_aa = cdr3)
  b_df <- df[str_detect(df$v_gene,'TRBV'), ] %>% rename(v_b_gene = v_gene, j_b_gene = j_gene, cdr3_b_aa = cdr3)

  b_df$v_b_gene <- paste0(b_df$v_b_gene,'*01') #needed for python tcrdist3 to infer cdr1,2 sequences
  b_df_out <- b_df
  b_df_out$count <- rep(1,nrow(b_df_out))
  b_df_out <- b_df_out %>% group_by(clonotype_id, v_b_gene, j_b_gene, cdr3_b_aa) %>% summarise(count = sum(count))
  write.csv(b_df_out, file = 'beta-unpaired-tcr.csv',row.names = FALSE)

  a_df$v_a_gene <- paste0(a_df$v_a_gene,'*01') #needed for python tcrdist3 to infer cdr1,2 sequences
  a_df_out <- a_df
  a_df_out$count <- rep(1,nrow(a_df_out))
  a_df_out <- a_df_out %>% group_by(clonotype_id, v_a_gene, j_a_gene, cdr3_a_aa) %>% summarise(count = sum(count))
  write.csv(a_df_out, file = 'alpha-unpaired-tcr.csv',row.names = FALSE)

  #can then add in a count column for rows of matching values (need to determine what to use as the match)

  full_df <- inner_join(a_df, b_df, by = 'clonotype_id') %>% distinct(barcode.x, clonotype_id, v_a_gene, j_a_gene, cdr3_a_aa, v_b_gene, j_b_gene, cdr3_b_aa, .keep_all=TRUE)
  full_df$count <- rep(1,nrow(full_df))
  full_df <- full_df %>% group_by(clonotype_id, v_a_gene, j_a_gene, cdr3_a_aa, v_b_gene, j_b_gene, cdr3_b_aa) %>% summarise(count = sum(count))
  full_df$clonotype_id <- make.unique(full_df$clonotype_id)

  write.csv(full_df, file = 'paired-tcr.csv')
}

