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
#' @import Seurat
#'
#' @export
#'
#' @examples
#' # library
#' library(Seurat)
#'
#' # analyzed data with tcr
#' data("seu_analyzed")
#'
#' # get clonotype highlighted on embedding
#' plot_embed_clonotype(seu_analyzed, clonotype_by = "cdr3s_aa", clonotype_id = "TRB:CASNPSTGGGGEQYF;TRA:CAVMDSNYQLIW",
#'                      color_by = "seurat_clusters")
#' plot_embed_clonotype(seu_analyzed, clonotype_by = "cdr3s_aa", clonotype_id = "TRB:CASNPSTGGGGEQYF;TRA:CAVMDSNYQLIW",
#'                      color_by = "seurat_clusters", reduction = "tsne", label = T)
#'
plot_embed_clonotype <- function (input, title = "", clonotype_id, clonotype_by, color_by, facet_by = NULL, ncol = NULL, label = FALSE,
                                  size = 1.5, alpha = 0.3, colors = NULL, theme = "classic", legend_dot_size = 1.5,
                                  xcol = "x", ycol = "y", reduction = "umap", text_sizes = c(20, 10, 5, 10), shuffle = F)
{
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
