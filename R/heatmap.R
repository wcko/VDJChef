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
#' plot_heatmap_clonotypes(seu_analyzed, clonotype_by = "cdr3s_aa", group_by = "seurat_clusters")
#' plot_heatmap_clonotypes(seu_analyzed, clonotype_by = "cdr3s_aa", group_by = "seurat_clusters", ntop = 10)
#'
#' # some additional uses
#' plot_heatmap_clonotypes(pfizer, clonotype_by = "CTaa", patient_by = "Patient", sample_by = "Sample") # ClonotypeIDs not listed, as not directly comparable across different patients
#' plot_heatmap_clonotypes(pfizer, clonotype_by = "CTaa", group_by = "Sample_Type")
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
    total <- rowSums(clonotypes_vs_color)
    clonotypes_vs_color <- clonotypes_vs_color[order(total, decreasing = TRUE)[1:ntop],, drop = FALSE]
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
