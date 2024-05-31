# Visualization functions.

#' Generate a grouped heatmap for a specific list of genes and corresponding
#' expression data. Just a wrapper around pheatmap::pheatmap.
#' 
#' @param assay The assay to include. Should be a 2D table with a column of gene symbols and
#'              at least one column of corresponding gene expression values for that symbol.
#'              Can also include a 'type column', or a column including the categories under
#'              which each gene falls (if you're trying to split a heatmap based on specific
#'              categories of genes). Ideally a data.table or a data.frame.
#' @param assay_data_columns A vector of names specifying the data columns within assay.
#'                           The default (NULL) just takes all columns besides the symbol
#'                           and type columns.
#' @param type_column The name of the column of assay specifying the category that each gene falls into.
#'                    If NULL, then no types are assumed and the heatmap won't be split.
#'                    Default is "TYPE".
#' @param symbol_column The column of assay specifying gene symbols. Default is "SYMBOL".
#' @param gap_length The space between splits on the split heatmap. Default is 3.
#' @param cluster_rows Whether to cluster heatmap rows within splits. Default is FALSE.
#' @param cluster_cols Whether to cluster columns across all heatmap entries, disregarding
#'                     splits. Default is FALSE.
#' @param border_color The color of border to use on heatmaps. Default is NA, or no color.
#' @param ... Other arguments to be passed into pheatmap::pheatmap.
#' 
#' @return Heatmap object as output by pheatmap::pheatmap.
#'
#' @export
create_split_heatmap <- function(assay,
                                 assay_data_columns = NULL,
                                 type_column = "TYPE",
                                 symbol_column = "SYMBOL",
                                 gap_length = 3,
                                 cluster_rows = FALSE,
                                 cluster_cols = FALSE,
                                 border_color = NA,
                                 ...) {
  # We assume a properly formatted assay is input already. No NAs, no nothing.
  if (any(is.na(assay))) {
    warning("There are NAs in the assay. This function might fail; check results carefully.")
  }

  # We assume there is a type column. If there is no type, no worries; just alert us.
  if (is.na(type_column)) {
    warning("No type column set. No splits will be made in the heatmap; just make sure that's what you want.")
  }

  # Set the data columns if not specified already.
  if (is.null(assay_data_columns)) {
    assay_data_columns <- setdiff(colnames(assay),
                                  c(type_column, symbol_column))
  }

  if (is.null(type_column)) {
    # Skip the ordering and splitting.
    heatmaps <- list(data.table::as.data.table(assay))
  } else { 
    # Order the assay
    assay <- data.table::as.data.table(assay[order(assay[[type_column]]),])

    # Split the genes of interest into tables based on their type.
    heatmaps <- split(assay, by = type_column)
  }

  # Determine the clustered orders of each row, find the lengths, and offset
  # them.
  maps_row_orders <- lapply(heatmaps, function (heatmap_data) {
    # Set up the data.
    assay <- as.data.frame(heatmap_data[, assay_data_columns])
    rownames(assay) <- heatmap_data[[symbol_column]]
    inds <- NULL

    # Determine indexing order.
    inds <- NULL
    if (dim(assay)[2] == 1) {
      inds <- if (cluster_rows) order(as.vector(unlist(assay)),
                                      decreasing = TRUE) else
        seq_len(dim(assay)[1])
    } else {
      map <- pheatmap::pheatmap(assay, cluster_rows = cluster_rows,
                                cluster_cols = FALSE, silent = TRUE, ...)
      inds <- if (cluster_rows) map$tree_row$order else seq_len(dim(assay)[1])
    }

    # Return index.
    return(inds)
  })
  maps_row_lengths <- lapply(maps_row_orders, length)
  maps_row_orders_offset <- lapply(seq_len(length(maps_row_orders)),
    function(index) {
      offsets <- as.vector(unlist(maps_row_lengths))
      return(maps_row_orders[[index]] + sum(offsets[1:index-1]))
    })
  maps_row_orders_offset <- as.vector(unlist(maps_row_orders_offset))

  # Reorder the assay.
  assay <- as.data.frame(assay[maps_row_orders_offset,])
  rownames(assay) <- assay[[symbol_column]]

  # Create gap columns.
  if (is.null(type_column)) {
    gaps <- numeric()
  } else {
    gaps <- which(diff(as.numeric(as.factor(assay[[type_column]]))) != 0)
    gaps <- as.vector(sapply(gaps, function(x) { rep(x, gap_length) }))
  }

  # Plot and return the heatmap.
  assay_plot <- as.data.frame(assay[,assay_data_columns])
  rownames(assay_plot) <- rownames(assay)
  colnames(assay_plot) <- if (is.numeric(assay_data_columns)) colnames(assay)[assay_data_columns]
                          else assay_data_columns
  return(pheatmap::pheatmap(assay_plot,
                            cluster_rows = FALSE,
                            cluster_cols = cluster_cols,
                            gaps_row = gaps,
                            border_color = border_color,
                            ...))
}

#' A simple function for creating dot plots for gene ontology (GO) or Ingenuity Pathway Analysis
#' (IPA) plots. Code taken from enrichplot: https://github.com/YuLab-SMU/enrichplot/blob/devel/R/dotplot.R#L200-L207
#' 
#' @param tbl A data.frame or similar containing gene symbols, z-scores, gene set ratios,
#'            and -log10(p) values. Both of these should be output or obtainable from GO enrichment
#'            analysis or IPA results.
#' @param z_score_colname Name of the column containing results for z-scores.
#' @param ratio_colname Name of the column containing gene set ratios.
#' @param neg_log_p_colname Name of the column containing -log10(p) values.
#' @param x_axis_label_colname Name of the column containing gene symbols for plotting along the
#'                             x-axis. This will create a horizontal plot. Note that either this
#'                             or y_axis_label_colname should be set but not both.
#' @param y_axis_label_colname Name of the column containing gene symbols for potting along the
#'                             y-axis. This will create a vertical plot. Note that either this
#'                             or y_axis_label_colname should be set but not both.
#' 
#' @returns A ggplot2-compatible dotplot summarizing the input data.
#' 
#' @export
ipa_dotplot <- function(tbl,
                        z_score_colname,
                        ratio_colname,
                        neg_log_p_colname,
                        x_axis_label_colname = NULL,
                        y_axis_label_colname = NULL) {
    # Validate input.
    if (is.null(x_axis_label_colname) == is.null(y_axis_label_colname)) {
      stop(print("You can only set x- or y-axis label names for IPA plots, not both."))
    } else if (!is.null(x_axis_label_colname)) {
      label_colname <- x_axis_label_colname
    } else {
      label_colname <- y_axis_label_colname
    }

    # Determine if the plot is horizontal or vertical.
    is_horizontal <- is.null(y_axis_label_colname)

    # Arrange so the table descends by z-score
    tmp <- dplyr::arrange(tbl, desc(!!sym(z_score_colname)))
    tmp[[label_colname]] <- ordered(tmp[[label_colname]],
                                           levels=rev(tmp[[label_colname]]))
    
    # Form and return the dot plot.
    if (is_horizontal) {
      return(ggplot2::ggplot(tmp,
                            ggplot2::aes(x = !!sym(ifelse(is_horizontal, label_colname, z_score_colname)),
                                          y = !!sym(ifelse(is_horizontal, z_score_colname, label_colname)),
                                          size = !!sym(ratio_colname),
                                          color = !!sym(neg_log_p_colname))) + 
              ggplot2::geom_point() + 
              ggplot2::theme_bw() + # can remove this for default ggplot theme
              ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    } else {
      return(ggplot2::ggplot(tmp,
                            ggplot2::aes(x = !!sym(ifelse(is_horizontal, label_colname, z_score_colname)),
                                          y = !!sym(ifelse(is_horizontal, z_score_colname, label_colname)),
                                          size = !!sym(ratio_colname),
                                          color = !!sym(neg_log_p_colname))) + 
              ggplot2::geom_point() + 
              ggplot2::theme_bw() + # can remove this for default ggplot theme
              ggplot2::theme(axis.text.y = element_text(angle = 0, vjust = 1)))
    }
}

#' Plot a nice-looking volcano plot with minimal effort. Thanks, Elvis.
#' 
#' @param df The data frame to input. This would ideally be the output of a program like
#'           DESeq2, containing logarithmic fold changes and adjusted p-values.
#' @param logfc_column The column name containing the logarithmic fold-changes.
#' @param logfc_cutoff The logarithmic fold-change cutoff that determines "high" or "low" values.
#'                     This will be used on both positive and negative ends, so this cutoff should be positive.
#' @param padj_column The column name containing adjusted p-values.
#' @param padj_cutoff The default p-value cutoff. Default is 0.05.
#' @param point_shape The shape of the points to plot. Default is 21, corrsponding to circles.
#'                    See ggplot2 for more information.
#' @param hline_linetype The type of boundary to plot between significant and nonsignificant points.
#'                       Default is a dashed line.
#' @param vlines_linetype The type of boundary to plot between high/low points and points with lower
#'                        logarithmic fold changes. Default is a dashed line.
#' 
#' @return A ggplot2-compatible volcano plot.
#' 
#' @export
volcano_plot <- function(df,
                         logfc_column,
                         logfc_cutoff,
                         padj_column,
                         padj_cutoff = 0.05,
                         point_shape = 21,
                         hline_linetype = "dashed",
                         vlines_linetype = "dashed") {
  `%>%` <- magrittr::`%>%`
  if (logfc_cutoff < 0) {
    warning("logfc_cutoff is less than 0, making positive.")
  }
  logfc_cutoff <- abs(logfc_cutoff)
  df <- df %>% dplyr::mutate(gene_type = case_when(!!sym(logfc_column) >= logfc_cutoff & !!sym(padj_column) < padj_cutoff ~ "up",
                                                  !!sym(logfc_column) <= -logfc_cutoff & !!sym(padj_column) < padj_cutoff ~ "down",
                                                  !!sym(logfc_column) > 0 & !!sym(logfc_column) < logfc_cutoff & !!sym(padj_column) < padj_cutoff ~ "warm",
                                                  !!sym(logfc_column) < 0 & !!sym(logfc_column) > -logfc_cutoff & !!sym(padj_column) < padj_cutoff ~ "cold",
                                                  TRUE ~ "ns"))
  plt <- df %>% ggplot2::ggplot(ggplot2::aes(x = !!sym(logfc_column),
                                             y = -log10(!!sym(padj_column)),
                                             fill = gene_type,
                                             size = gene_type,
                                             alpha = gene_type)) + 
                ggplot2::geom_point(shape = point_shape, color = "black") + 
                ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = hline_linetype) + 
                ggplot2::geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = vlines_linetype) + 
                ggplot2::scale_fill_manual(values = c("up" = "red", "down" = "blue", "ns" = "gray", "warm" = "pink", "cold" = "skyblue")) + 
                ggplot2::scale_size_manual(values = c("up" = 2, "down" = 2, "ns" = 1, "warm" = 1, "cold" = 1)) + 
                ggplot2::scale_alpha_manual(values = c("up" = 1, "down" = 1, "ns" = 0.5, "warm" = 0.5, "cold" = 0.5))

  return(plt)
}