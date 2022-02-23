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
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
    coord <- input@reductions[[reduction]]@cell.embeddings
    tmp[[xcol]] <- coord[,1]
    tmp[[ycol]] <- coord[,2]
  } else {
    tmp <- pData(input)
  }
  tmp_clonotype <- tmp[tmp[[clonotype_by]] %in% clonotype_id,]

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
  g <- g + guides(colour = guide_legend(override.aes = list(size = legend_dot_size)))

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
    color_by_freq <- table(tmp_clonotype[, color_by])
    color_by_freq <- color_by_freq[order(color_by_freq, decreasing = TRUE)]
    label_length <- sapply(names(color_by_freq), nchar)
    label_sep <- max(label_length) - label_length + 2
    label_sep <- sapply(label_sep, function(x) return(paste(rep(" ", x), collapse="")))
    label_color <- paste0(names(color_by_freq), label_sep, color_by_freq)
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
#' Get top clonotypes categorized by a column in metadata
#'
#' @param input ExpressionSet or Seurat Object
#' @param clonotype_by the label of the column in metadata that defines clonotypes
#' @param group_by what to group clonotypes frequencies
#' @param ntop the number of top clonotypes
#' @param split_by split tables by
#' @param chain_type subset table to chain type
#'
#' @export
#'
#' @examples
get_topclonotypes <- function(input, clonotype_by, group_by, ntop = 20, split_by = NULL, chain_type = NULL){

  # get metadata, specific to clonotypes
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
  } else {
    tmp <- pData(input)
  }
  tmp <- tmp[!is.na(tmp[[clonotype_by]]),]

  if(!is.null(chain_type))
    tmp <- tmp[tmp[["chain_summary"]]==chain_type,]

  # data for heatmap
  if(is.null(split_by)){
    clonotypes_vs_color <- table(tmp[[clonotype_by]], tmp[[group_by]])
    clonotypes_vs_color <- data.frame(rbind(clonotypes_vs_color))
    clonotypes_vs_color$Total <- rowSums(clonotypes_vs_color)
    clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color$Total, decreasing = TRUE)[1:ntop],]
  } else {
    tmp <- split(tmp, list(tmp[[split_by]]))
    clonotypes_vs_color <- lapply(tmp, function(x){
      x <- table(x[[clonotype_by]], x[[group_by]])
      x <- data.frame(rbind(x))
      x$Total <- rowSums(x)
      x <- x[order(x$Total, decreasing = TRUE)[1:ntop],]
    })
  }
  return(clonotypes_vs_color)
}

#' plot_heatmap_clonotypes
#'
#' Visual the clonotype expansion in an heatmap
#'
#' @param input ExpressionSet or Seurat Object
#' @param clonotype_by the label of the column in metadata that defines clonotypes
#' @param color_by What to color points by, either "UMI_sum", or pData categorial variable, ignored if gene is provided
#' @param ntop the number of top clonotypes
#' @param Count_limit count limit for plot legend
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @import grid
#'
#' @export
#'
#' @examples
plot_heatmap_clonotypes <- function(input, clonotype_by, color_by, ntop = 20, Count_limit = NULL){

  # get metadata, specific to clonotypes
  if(class(input) == "Seurat"){
    tmp <- input@meta.data
  } else {
    tmp <- pData(input)
  }
  tmp <- tmp[!is.na(tmp[[clonotype_by]]),]

  # data for heatmap
  clonotypes_vs_color <- table(tmp[[clonotype_by]], tmp[[color_by]])
  clonotypes_vs_color <- data.frame(rbind(clonotypes_vs_color))
  clonotypes_vs_color_Total <- rowSums(clonotypes_vs_color)
  clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color_Total, decreasing = TRUE)[1:ntop],]

  # heatmap auxiliary tools
  col_fun = colorRamp2(c(max(clonotypes_vs_color), 0), c("red","white"))
  if(!is.null(Count_limit))
    col_fun = colorRamp2(c(Count_limit, 0), c("red","white"))

  # heatmap
  count_heatmap <- Heatmap(clonotypes_vs_color, col = col_fun,
                           row_names_side = "left", column_names_side = "top",
                           column_names_rot = 45, cluster_columns = FALSE, cluster_rows = FALSE,
                           heatmap_legend_param = list(title = ""), column_title = NULL,
                           column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(paste0(sprintf("%.0f", clonotypes_vs_color[i, j])), x, y, gp = gpar(fontsize = 10))
                           })

  return(count_heatmap)
}
