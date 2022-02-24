#' plot_tcr_violin_seurat - modification of function from Garber Lab SignallingSingleCell package for input of Seurat Object/TCR Clonotype data 
#' to compare gene expression of Expanded clones vs background Non-expanded clones
#' 
#' @description Create violin plot.
#' @param input SatijaLab’s Seurat Class, with normalized expression values in assay data slot, and TCR Clonotype ID's in meta.data. Supports RNA and ADTs via Seurat::FetchData()
#' @param title Title of the graph. Would be the gene name followed by clone's ID if not specified
#' @param gene Feature for which to plot the expression level.
#' @param clone specific name of expanded clone of interest to compare to the non-expanded clones
#' @param clonotype_id meta.data column name for clonotype ID's
#' @param threshold number of clones above which is considered expanded (e.g. for n=5, 5 clones and below are Non-expanded)
#' @param log_scale If true, transform UMIs by log2(UMI + 1).
#' @param colors What colors to utilize for categorical data. Be sure it is of the proper length.
#' @param facet_by a vector with one or two meta.data column variables. If two, the first variable as columns and the second as rows.
#' @param spread e.g. Healthy category is unique in Disease and Skin. To use Healthy only as skin but not Disease, that is adding Healthy skin to each disease, spread = c("Disease", "Healthy").
#' @param text_sizes a vector of title_size, axis_title, axis_text, legend_title, legend_text, facet_text, number_label_text_size, defaults too c(20,10,5,10,5,5,2)
#' @param theme the plot theme. Default to be "classic" if not set to "bw".
#' @param number_labels show the total cell numbers and cell fraction with non-zero expression values under each bar.
#' @param plot_mean plot the mean value as black dot with second y-axis on the right.
#' @param size the size of dots.
#' @param sig the number of digits after the decimal point for cell fraction value.
#' @param contour_line_width the thickness of the violin contour line
#' @details
#' Utilize information stored in meta.data to control the plot display. Each point_by as a dot with a bar showing the weighted mean of all point_by dots.
#' color_by is set to display Expanded and Non-expanded clones
#' @examples
#' plot_tcr_violin_seurat(ex_sc, gene = "GNLY", clone = "clonotype1_Control-BCH001", facet_by = c("hash.ID"))
#' plot_tcr_violin_seurat(ex_sc, gene = "GNLY", clone = "clonotype1_Control-BCH001", threshold=1, facet_by = c("hash.ID", "CellTypeSuperCluster"), log_scale = F)
#' @export
plot_tcr_violin_seurat <- function (input, title = "", gene, clone, clonotype_id="clonotype_id", threshold=5, log_scale = F,
                         colors = NULL, facet_by = NULL, spread = NULL, jitter_pts = T,
                         plot_mean = T, plot_mean_dot_size = 2, size = 1, sig = 3, number_labels = T,
                         text_sizes = c(20, 10, 5, 10, 5, 5, 2), alpha = 0.5, theme = "classic",
                         contour_line_width = 0.3)
{
  df <- input@meta.data[,colnames(input@meta.data) %in% c(clonotype_id, facet_by), drop = F]
  df <- df[!is.na(df$clonotype_id),]
  clonotypedf <- df %>% 
    count(clonotype_id) %>% 
    mutate(expanded_status = ifelse(n > threshold, "expanded", "non-expanded")) %>% 
    select(-n)
  df$barcodes <- rownames(df)
  df <- inner_join(df, clonotypedf, by=clonotype_id)
  df.sub <- df[df$clonotype_id==clone,]
  df.sub <- bind_rows(df.sub, df[df$expanded_status=='non-expanded',])
  df <- df.sub
  df <- bind_cols(df, raw=FetchData(input,gene,cells=df$barcodes)[,1])

  color_by <- "expanded_status"
  colnames(df) <- gsub("-", "", colnames(df))
  gene <- gsub("-", "", gene)
  if (any(!is.null(spread))) {
    others <- setdiff(unique(df[,spread[1]]), spread[2])
    ind <- which(df[, spread[1]] == spread[2])
    rmdf <- df[ind,]
    df <- df[-ind,]
    for (i in 1:length(others)) {
      rmdf[,spread[1]] <- others[i]
      df <- rbind(df, rmdf)
    }
  }
  if (log_scale == T) {df$plot <- log2(df$raw + 1)}else{df$plot <- df$raw}
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols <- gg_color_hue(length(unique(df[, color_by])))
  
  g <- ggplot(df)
  if (all(!is.null(colors))) {
    g <- g + scale_color_manual(values = c(colors))
    g <- g + scale_fill_manual(values = c(colors))
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }else{
    g <- g + theme_classic()
  }
  if (title == "") title <- paste0(gene,"-",clone)
  g <- g + labs(title = title, y = gene)
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
  if (jitter_pts == T) g <- g + geom_jitter(aes_string(x = color_by, y = "plot", col = color_by), width = 0.2, size = size)
  g <- g + geom_violin(aes_string(x = color_by, y = "plot", fill = color_by), col = "black", trim = T, scale = "width", alpha = alpha, size=contour_line_width) 
  if (number_labels == T) {
    g <- g + stat_summary(aes_string(x = color_by, y = "raw"), fun.data = function(x) {return(c(y = -max(df$plot)/25, label = length(x)))}, colour = "black",
                          geom = "text", size = text_sizes[7])
    g <- g + stat_summary(aes_string(x = color_by, y = "raw"), fun.data = function(x) {return(c(y = -max(df$plot)/10, label = round(mean(as.numeric(x > 0)), sig)))}, colour = "black",
                          geom = "text", size = text_sizes[7])
  }
  if (plot_mean == TRUE) {
    scale <- max(df$plot)/max(tapply(df$raw, INDEX = as.list(df[, colnames(df) %in% c(color_by, facet_by), drop = F]), FUN=mean), na.rm = T)
    g <- g + suppressWarnings(stat_summary(aes_string(x = color_by, y = "raw"), fun.y = function(x) mean(x)*(scale * 0.5), colour = "black", geom = "point", size = plot_mean_dot_size))
    g <- g + scale_y_continuous(sec.axis = sec_axis(~./(scale * 0.5), name = "Mean Expression"))
  }
  if (length(facet_by) == 1) {
    g <- g + facet_grid(facets = reformulate(facet_by), scales = "free_x", space = "free_x")
  }else if (length(facet_by) == 2) {
    g <- g + facet_grid(facets = reformulate(facet_by[1], facet_by[2]), scales = "free_x", space = "free_x")
  }else if (length(facet_by) > 2) {stop("Parameter facet_by needs to be a string with equal or less than two variables.")}
  if (!is.null(facet_by)) g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
  return(g)
}