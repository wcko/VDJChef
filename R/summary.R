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
#' # get clonotypes ordered by frequency
#' get_topclonotypes(seu_analyzed, clonotype_by = "cdr3s_aa", group_by = "seurat_clusters")
#' get_topclonotypes(seu_analyzed, clonotype_by = "cdr3s_aa", group_by = "seurat_clusters", ntop = 10)
#'
#' # additional commands with possible uses
#' get_topclonotypes(seu_analyzed, clonotype_by = "CTaa", patient_by = "Patient", sample_by = "Sample")
#' get_topclonotypes(seu_analyzed, clonotype_by = "CTaa", group_by = "Sample")
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
      clonotypes_vs_color <- clonotypes_vs_color[order(clonotypes_vs_color$Total, decreasing = TRUE)[1:ntop],, drop = FALSE]
    } else {
      tmp <- split(tmp, list(tmp[[patient_by]]))
      clonotypes_vs_color <- lapply(tmp, function(x){
        x <- table(x[[clonotype_by]], x[[sample_by]])
        x <- data.frame(rbind(x))
        x$Total <- rowSums(x)
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
      x <- x[order(x$Total, decreasing = TRUE)[1:ntop],, drop = FALSE]
    })
  }
  return(clonotypes_vs_color)
}






