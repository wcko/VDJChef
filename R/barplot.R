#' VDJ_visualize
#'
#' Visualize V, D or J gene frequencies across Patients or Samples in a bar graph.
#'
#' @param df
#' @param gene either v_gene, d_gene or j_gene
#' @param chain either TRA or TRB
#' @param plot only bar
#' @param y.axis
#' @param order
#' @param scale
#' @param split.by Which meta.data to split the data by.
#'
#' @return
#' @export
#'
#' @examples
VDJ_visualize <- function (df, gene = "v_gene", chain = "TRA", plot = "bar", y.axis = "sample",
                           order = "gene", scale = TRUE, split.by = NULL){


  #### subset the seurat object with only the 3 meta.data columns we want
  df2 <- df@meta.data[ , c(gene, "Sample", "Patient", "Sample_Subtype")]

  #### split TRA and TRB in v_gene column
  df2[[gene]] <- strsplit(df2[[gene]], split = ";")

  #### remove TRB data... now doesn't work when not splitting with ";"
  df2[[gene]] <- sapply(df2[[gene]], function(x){x[grepl(chain,x)]})

  #### Remove rows with character(0) and TRAVs with more than one chain
  df2 <- df2[lengths(df2[[gene]]) > 0 & lengths(df2[[gene]]) <= 1,]

  #### Make df2$v_gene into original format
  df2[[gene]] <- paste0(df2[[gene]])

  #### Proportion of barcodes that express the V gene in column 1 within the Sample and the total barcodes in the Sample
  df2$gene_freq <- NA
  for (i in 1:nrow(df2)) {
    df2.row <- df2[i,]
    row.gene <- df2[[gene]][i]
    row.sample <- df2.row$Sample
    gene_freq <- df2 %>% filter(Sample == row.sample) %>% group_by_at(vars(gene)) %>% dplyr::summarize(cell_count = dplyr::n())
    gene_freq$Sample = row.sample
    gene_freq$gene_freq <- gene_freq$cell_count/sum(gene_freq$cell_count)
    df2[i, "gene_freq"] <- gene_freq$gene_freq[gene_freq[[gene]] == row.gene]
    rm(row.gene_freq)
  }

  #### varcount:
  df3 <- df2
  df3$varcount <- NA
  df3  <- df3 %>% group_by_at(vars(c("Patient", "Sample", gene))) %>% dplyr::summarize(count = dplyr::n())
  df3 <- df3[,-which(colnames(df3) == "count")]
  df3$gene_freq <- NA
  for (i in 1:nrow(df3)) {
    df3.row <- df3[i,]
    row.patient <- df3.row$Patient
    row.vgene <- df3.row[[gene]]
    row.sample <- df3.row$Sample
    row.gene_freq <- df2 %>% filter(Sample == row.sample & Patient == row.patient & (!!as.name(gene)) == row.vgene) %>% dplyr::select(gene_freq) %>% unlist() %>% unique() %>% as.numeric()
    df3$gene_freq[i] <- row.gene_freq
  }

  df3 <- df3 %>% group_by_at(vars(c("Patient", gene))) %>% dplyr::summarise(varcount = sum(gene_freq),
                                                                            sd = sd(gene_freq),
                                                                            mean = mean(gene_freq)) %>% full_join(df3,.)
  df3 <- df3 %>% arrange(gene, Patient, Sample)

  #### graphing it in a bar graph
  df4 <- df3[, c(gene, "Sample","Patient", "sd", "mean")]
  df4 <- na.omit(df4)

  ggplot2::ggplot(df4, aes(x = !!as.name(gene), y = mean)) + geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean, ymax = mean + sd),
                  width = 0.2, position = position_dodge(0.9)) +
    theme_classic() + theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1, size = rel(0.5))) + facet_grid(Patient ~.) + geom_hline(yintercept = 0)

}
