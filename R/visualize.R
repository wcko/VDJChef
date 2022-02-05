#' plot_embed_clonotype
#'
#' Visual and highlight clonotypes given an embedding
#'
#' @param input ExpressionSet or Seurat Object
#' @param title Title of the ggplot2 plot
#' @param clonotype_id Clonotype ID among values listed by "clonotype_by" parameter
#' @param clonotype_by the label of the column in metadata that defines clonotypes
#' @param color_by
#' @param facet_by
#' @param ncol
#' @param label
#' @param size
#' @param alpha
#' @param colors
#' @param theme
#' @param legend_dot_size
#' @param xcol
#' @param ycol
#' @param text_sizes
#' @param shuffle
#'
#' @return
#' @export
#'
#' @examples
plot_embed_clonotype <- function (input, title = "", clonotype_id, clonotype_by, color_by, facet_by = NULL, ncol = NULL, label = FALSE,
                                  size = 1.5, alpha = 0.3, colors = NULL, theme = "classic", legend_dot_size = 1.5,
                                  xcol = "x", ycol = "y", text_sizes = c(20, 10, 5, 10), shuffle = F)
{
  # get metadata
  tmp <- pData(input)
  tmp_clonotype <- pData(input)[pData(input)[[clonotype_by]] %in% clonotype_id,]

  # shuffle ?
  if (shuffle) {
    tmp_clonotype <- tmp_clonotype[sample(nrow(tmp_clonotype)), ]
  }
  else {
    tmp_clonotype <- tmp_clonotype[order(tmp_clonotype[, color_by]), ]
  }

  # plot and theme setting
  if (title == "") title <- paste(clonotype_by, clonotype_id, sep = " = ")
  g <- ggplot(tmp_clonotype)
  g <- g + labs(title = title, x = "tSNE[1]", y = "tSNE[2]")
  if (theme == "bw") {
    g <- g + theme_bw()
  }
  else {
    g <- g + theme_classic()
  }
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]))
  g <- g + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  g <- g + guides(colour = guide_legend(override.aes = list(size = legend_dot_size)))

  # grey background
  tmp <- pData(input)[c(xcol, ycol)]
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
    # g <- g + scale_color_manual(values = names(color_by_freq), labels = paste(names(color_by_freq), color_by_freq, sep = ", "))
    g <- g + scale_colour_discrete(limits = names(color_by_freq), labels = label_color)
    if (!is.null(colors)) {
      g <- g + scale_color_manual(values = c(colors))
    }
    if(label){
      color_by_levels <- levels(factor(pData(input)[, color_by]))
      coords <- pData(input)[,c(xcol,ycol)]
      coords <- aggregate(coords, list(pData(input)[[color_by]]), median)
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
#' @param input
#' @param clonotype_by
#' @param group_by
#' @param ntop
#' @param split_by
#' @param chain_type
#' @param Count_limit
#'
#' @return
#' @export
#'
#' @examples
get_topclonotypes <- function(input, clonotype_by, group_by, ntop = 20, split_by = NULL, chain_type = NULL, Count_limit = NULL){

  # get metadata, specific to clonotypes
  tmp <- pData(input)
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
#' @param input
#' @param clonotype_by
#' @param color_by
#' @param ntop
#' @param Count_limit
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap_clonotypes <- function(input, clonotype_by, color_by, ntop = 20, Count_limit = NULL){

  # get metadata, specific to clonotypes
  tmp <- pData(input)
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
