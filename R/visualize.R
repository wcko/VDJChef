#' plot_embed_clonotype
#'
#' Visual and highlight clonotypes given an embedding
#'
#' @param input ExpressionSet or Seurat Object
#' @param title Title of the ggplot2 plot
#' @param clonotype_id Clonotype ID among values listed by "clonotype_by" parameter
#' @param clonotype_by the label of the column in metadata that defines clonotypes
#' @param label if true, labels of "color_by" will be visualized on embedding
#' @param color_by What to color points by, either "UMI_sum", or pData categorial variable, ignored if gene is provided
#' @param facet_by What to break the plots by
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param alpha The transparency of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @param theme display theme
#' @param legend_dot_size Size of dot in legend
#' @param xcol pData column to use for x axis
#' @param ycol pData column to use for y axis
#' @param reduction the reduction key for Seurat objects
#' @param text_sizes text sizes for multiple components of the ggplot figure
#' @param shuffle change order of labels
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
plot_embed_clonotype <- function (input, title = "", clonotype_id, clonotype_by, color_by, facet_by = NULL, ncol = NULL, label = FALSE,
                                  size = 1.5, alpha = 0.3, colors = NULL, theme = "classic", legend_dot_size = 1.5,
                                  xcol = "x", ycol = "y", reduction = "umap", text_sizes = c(20, 10, 5, 10), shuffle = F)
{
  # get metadata
  # get metadata, specific to clonotypes
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
    coord <- input@reductions[[reduction]]@cell.embeddings
    tmp[[xcol]] <- coord[,1]
    tmp[[ycol]] <- coord[,2]
  } else if (class(input) == "ExpressionSet") {
    tmp <- pData(input)
  } else {
    print("Input is neither a Seurat or ExpressionSet Object")
  }
  tmp_clonotype <- tmp[which(tmp[[clonotype_by]] == clonotype_id),]

  # shuffle ?
  if (shuffle) {
    tmp_clonotype <- tmp_clonotype[sample(nrow(tmp_clonotype)), ]
  }
  else {
    tmp_clonotype <- tmp_clonotype[order(tmp_clonotype[, color_by]), ]
  }

  # initial plot
  g <- ggplot(tmp_clonotype)
  if (theme == "bw") {
    g <- g + theme_bw()
  }
  else {
    g <- g + theme_classic()
  }

  # title setting
  if (title == "") title <- paste(clonotype_by, clonotype_id, sep = " = ")

  # x and y labs
  if(class(input) == "Seurat"){
    main_lab <- reduction
    g <- g + labs(title = title, x = paste0(reduction,"[1]"), y = paste0(reduction,"[2]"))
  } else {
    g <- g + labs(title = title, x = "tSNE[1]", y = "tSNE[2]")
  }

  # other features
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]))
  g <- g + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  g <- g + guides(colour = guide_legend(title = paste0(color_by, " (Freq.)"), override.aes = list(size = legend_dot_size)))

  # grey background
  g <- g + geom_point(data = tmp, aes_string(x = xcol, y = ycol), shape = 20, col = "gray", size = size,
                      alpha = alpha)

  # subset for clonotype and visualize
  if (class(tmp_clonotype[, color_by]) %in% c("double","integer","numeric")) {
    g <- g + geom_point(data = tmp_clonotype, aes_string(x = xcol, y = ycol, col = color_by),
                        shape = 20, size = size, alpha = alpha)
    g <- g + scale_color_gradientn(colours = c("gray", "blue", "red", "yellow"))
  }
  else {
    color_by_freq <- table(droplevels(factor(tmp_clonotype[, color_by])))
    color_by_freq <- color_by_freq[order(color_by_freq, decreasing = TRUE)]
    label_length <- sapply(names(color_by_freq), nchar)
    label_sep <- max(label_length) - label_length + 2
    label_sep <- sapply(label_sep, function(x) return(paste(rep(" ", x), collapse="")))
    label_color <- paste0(names(color_by_freq), label_sep, "(", color_by_freq,")")
    g <- g + geom_point(data = tmp_clonotype, aes_string(x = xcol, y = ycol, col = color_by),
                        shape = 20, size = size)
    g <- g + scale_colour_discrete(limits = names(color_by_freq), labels = label_color)
    if (!is.null(colors)) {
      g <- g + scale_color_manual(values = c(colors))
    }
    if(label){
      coords <- tmp[,c(xcol,ycol)]
      coords <- aggregate(coords, list(tmp[[color_by]]), median)
      g <- g + geom_text(data = coords, aes(x = x, y=y, label = Group.1))
    }
  }

  # visualize for facets
  if (!is.null(facet_by))
    g <- g + facet_wrap(facets = reformulate(facet_by), drop = TRUE, scales = "free_x", ncol = ncol)

  plot(g)
}

#' get_topclonotypes
#'
#' Return a table of clonal expansion. Default shows the clonal expansion of the top 20 unique clonotypes across all of a grouping variable,
#' or the number of clones of the top 20 clonotypes per each facet by multiple grouping variables
#'
#' @param input Seurat Object or ExpressionSet
#' @param clonotype_by meta.data or pData column name for clonotype ID's
#' @param group_by meta.data or pData column names for grouping specific clonotype frequencies, if multiple columns are entered, columns are merged by paste
#' @param patient_by meta.data or pData column name for patient, ignored if split_by is not NULL
#' @param sample_by meta.data or pData column name for sample, ignored if split_by is not NULL
#' @param ntop the number of top clonotypes to return
#' @param split_by split tables by
#' @param chain_type subset table to chain type
#'
#' @export
#'
#' @examples
#' get_topclonotypes(pfizer, clonotype_by = "CTaa", patient_by = "Patient", sample_by = "Sample") # ClonotypeIDs not listed, as not directly comparable across different patients
#' get_topclonotypes(pfizer, clonotype_by = "CTaa", group_by = "Sample")
#' get_topclonotypes(pfizer[,which(pData(pfizer)$Patient=="VB234")], clonotype_by = "CTaa", group_by="Sample")
#'
get_topclonotypes <- function(input, clonotype_by, group_by = NULL, patient_by = NULL, sample_by = NULL, ntop = 20, split_by = NULL, chain_type = NULL){

  # get metadata, specific to clonotypes
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
  } else if (class(input) == "ExpressionSet") {
    tmp <- pData(input)
  } else {
    print("Input is neither a Seurat or ExpressionSet Object")
  }
  tmp <- tmp[!is.na(tmp[[clonotype_by]]),]

  if(!is.null(chain_type))
    tmp <- tmp[tmp[["chain_summary"]]==chain_type,]


  # merge two grouping columns into one, do not use more than two grouping variables
  if(length(group_by) == 1){
    tmp[["group_clonotype"]] <- tmp[[group_by]]
  } else if(length(group_by) == 2){
    tmp[["group_clonotype"]] <- paste(tmp[[group_by[1]]], tmp[[group_by[2]]], sep = "_")
  } else if(length(group_by) > 2){
    stop("Using more than 2 grouping columns are not suggested")
  }

  # data for heatmap
  if(is.null(split_by)){
    if(is.null(patient_by) || is.null(sample_by)){
      clonotypes_vs_color <- table(tmp[[clonotype_by]], tmp[['group_clonotype']])
      clonotypes_vs_color <- data.frame(rbind(clonotypes_vs_color))
      clonotypes_vs_color$Total <- rowSums(clonotypes_vs_color)
      # clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color$Total, decreasing = TRUE)[1:ntop],-length(clonotypes_vs_color), drop = FALSE]
      clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color$Total, decreasing = TRUE)[1:ntop],, drop = FALSE]
    } else {
      tmp <- split(tmp, list(tmp[[patient_by]]))
      clonotypes_vs_color <- lapply(tmp, function(x){
        x <- table(x[[clonotype_by]], x[[sample_by]])
        x <- data.frame(rbind(x))
        x$Total <- rowSums(x)
        # x <- x[order(x$Total, decreasing = TRUE)[1:ntop],-ncol(x), drop = FALSE]
        x <- x[order(x$Total, decreasing = TRUE)[1:ntop],, drop = FALSE]
      })
      colnames_clono <- unlist(lapply(clonotypes_vs_color, colnames))
      clonotypes_vs_color <- do.call(cbind, clonotypes_vs_color)
      rownames(clonotypes_vs_color) <- NULL
      colnames(clonotypes_vs_color) <- colnames_clono
    }
  } else {
    tmp <- split(tmp, list(tmp[[split_by]]))
    clonotypes_vs_color <- lapply(tmp, function(x){
      x <- table(x[[clonotype_by]], x[['group_clonotype']])
      x <- data.frame(rbind(x))
      x$Total <- rowSums(x)
      # x <- x[order(x$Total, decreasing = TRUE)[1:ntop],-ncol(x), drop = FALSE]
      x <- x[order(x$Total, decreasing = TRUE)[1:ntop],, drop = FALSE]
    })
  }
  return(clonotypes_vs_color)
}

#' plot_heatmap_clonotypes
#'
#' Visualize clonal expansion by heatmap. Default shows the clonal expansion of the top 20 unique clonotypes across all of a grouping variable,
#' or the number of clones of the top 20 clonotypes per each facet by multiple grouping variables
#' Bottom row annotation shows the total number of clones in each grouping variable
#'
#' @param input Seurat Object or ExpressionSet
#' @param clonotype_by meta.data or pData column name for clonotype ID's
#' @param patient_by meta.data or pData column name for patient
#' @param sample_by meta.data or pData column name for sample
#' @param group_by What to color points by, either "UMI_sum", or meta.data or pData column name, ignored if gene is provided
#' @param sort_by The metadata column for sorting
#' @param sort_id The value in the metadata column for sorting
#' @param ntop the number of top clonotypes to return
#' @param Count_limit count limit for plot legend
#'
#' @importFrom ComplexHeatmap Heatmap columnAnnotation anno_text
#' @importFrom circlize colorRamp2
#' @import grid
#'
#' @export
#'
#' @examples
#' plot_heatmap_clonotypes(pfizer, clonotype_by = "CTaa", patient_by = "Patient", sample_by = "Sample") # ClonotypeIDs not listed, as not directly comparable across different patients
#' plot_heatmap_clonotypes(pfizer, clonotype_by = "CTaa", group_by = "Sample_Type")
#' plot_heatmap_clonotypes(pfizer[,which(pData(pfizer)$Patient=="VB234")], clonotype_by = "CTaa", group_by="Sample")
#'
plot_heatmap_clonotypes <- function(input, clonotype_by, patient_by = NULL, sample_by = NULL, group_by = NULL, sort_by = NULL, sort_id = NULL, ntop = 20, Count_limit = NULL) {
  # get metadata, specific to clonotypes
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
  } else if (class(input) == "ExpressionSet") {
    tmp <- pData(input)
  } else {
      print("Input is neither a Seurat or ExpressionSet Object")
  }
  tmp <- tmp[!is.na(tmp[[clonotype_by]]),]

  totalclones <- NULL
  # data for heatmap
  if(is.null(patient_by) || is.null(sample_by)){
    stopifnot(!is.null(group_by))
    clonotypes_vs_color <- table(tmp[[clonotype_by]], tmp[[group_by]])
    clonotypes_vs_color <- data.frame(rbind(clonotypes_vs_color))
    totalclones <- colSums(clonotypes_vs_color)
    clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color$Total, decreasing = TRUE)[1:ntop],-length(clonotypes_vs_color), drop = FALSE]
    heatmap_split <- NULL
  } else {
    tmp.bak <- tmp
    tmp <- split(tmp, list(tmp[[patient_by]]))
    clonotypes_vs_color <- lapply(tmp, function(x){
      if(is.null(sort_by) || is.null(sort_id)){
        x <- table(x[[clonotype_by]], x[[sample_by]])
        x <- data.frame(rbind(x))
        x$Total <- rowSums(x)
      } else {
        sorter_subset <- unique(x[[sample_by]][x[[sort_by]]==sort_id])
        print(sorter_subset)
        x <- table(x[[clonotype_by]], x[[sample_by]])
        x <- data.frame(rbind(x))
        print(head(x[,sorter_subset, drop = FALSE]))
        if(length(sorter_subset) > 0) x$Total <- rowSums(x[,sorter_subset, drop = FALSE])
        else x$Total <- rowSums(x)
      }
      x <- x[order(x$Total, decreasing = TRUE)[1:ntop],-ncol(x), drop = FALSE]
    })
    heatmap_split <- rep(names(clonotypes_vs_color), sapply(clonotypes_vs_color, ncol))
    colnames_clono <- unlist(lapply(clonotypes_vs_color, colnames))
    clonotypes_vs_color <- do.call(cbind, clonotypes_vs_color)
    rownames(clonotypes_vs_color) <- NULL
    colnames(clonotypes_vs_color) <- colnames_clono
    totalclones <- table(tmp.bak[[clonotype_by]], tmp.bak[[sample_by]])
    totalclones <- data.frame(rbind(totalclones))
    totalclones <- colSums(totalclones)
  }

  # heatmap auxiliary tools
  clonotypes_vs_color[is.na(clonotypes_vs_color)] <- 0
  col_fun = colorRamp2(c(max(clonotypes_vs_color), 0), c("red","white"))
  if(!is.null(Count_limit))
    col_fun = colorRamp2(c(Count_limit, 0), c("red","white"))

  # heatmap
  count_heatmap <- Heatmap(clonotypes_vs_color, col = col_fun,
                           row_names_side = "left", column_names_side = "top",
                           column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
                           bottom_annotation =
                             columnAnnotation(Total_Clones_in_Sample = anno_text(totalclones, gp = gpar(fontsize=10), rot=0, just = "top", location=1)),
                           heatmap_legend_param = list(title = ""), border = TRUE, column_split = heatmap_split,
                           column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(paste0(sprintf("%.0f", clonotypes_vs_color[i, j])), x, y, gp = gpar(fontsize = 10))
                           })
  return(count_heatmap)
}




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

## get_silhouette
#' Returns an object 'sil' of class 'silhouette' (from package cluster), which is an n x 3 matrix.
#' For each obs i, the first col is its cluster identity, second col is closest neighbor,
#' and third col is silhouette score (-1,1) where -1: clusters assigned wrong way, 0: clusters indifferent, 1: clusters well apart
#' Requires PCA matrix embedding to have been calculated on the input object
#'
#' @param input Seurat Object or ExpressionSet, where PCA matrix has been calculated and assigned to PCA slot
#' @param pcs Number of principal components to use
#' @param distmetric Distance metric "euclidean", "manhattan"
#' @param anno column name in meta.data or pData column name that refers to the cluster annotations
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @export
#'
#' @examples
#' get_silhouette(pfizer, anno = 'cluster_ID')
#'
get_silhouette <- function(input, pcs = 10, distmetric = 'euclidean', anno = 'seurat_clusters') {
  maxpcs <- ncol(input[['pca']]@cell.embeddings)
  if(maxpcs < pcs){
    print("Maximum number of PCs in embedding is less than that specified. Setting maximum number of PCs to that in embedding")
    pcs <- maxpcs
  }
  x <- as.matrix(Embeddings(object=input[['pca']])[,1:pcs])
  distances <- dist(x, method=distmetric)
  distancemat <- as.matrix(distances)
  input@meta.data[,anno] <- droplevels(input@meta.data[,anno])
  clusterids <- as.numeric(input@meta.data[,anno]) #change factor names to numeric as required for silhouette function
  sil <- silhouette(clusterids,dmatrix=distancemat)
  annonames <- levels(input@meta.data[,anno])
  suppressWarnings(sil[,1:2] <- annonames[match(sil[,1:2],levels(factor(clusterids)),annonames)]) #convert numerics back to names if applicable
  return(sil)
}


