# Cell imputation functions.

#' Call cell types in a Seurat object using adobo, a Python package that calls
#' cell types using a marker gene list.
#' 
#' @param st_obj Seurat object to process.
#' @param gene_symbols A vector of gene symbols corresponding to the marker
#'                     genes for each cluster.
#' @param cell_types A vector of cell types corresponding to the marker genes
#'                   for each cluster.
#' @param cluster_slot The slot name in the Seurat dataset corresponding to
#'                     cell identities or cluster names. The default behavior
#'                     is Idents(st_obj).
#' @param symbol_switch Whether to switch ENSEMBL to MGI symbols for "mouse"
#'                      or "human". Default is to not switch (None). Input
#'                      gene format must be gene symbols in that case.
#' @param new_idents Whether to set the adobo output results as new cell
#'                   idents. Default is FALSE.
#' 
#' @return The modified Seurat object with adobo clusters.
#' 
#' @export
adobo_impute <- function(st_obj,
                         gene_symbols,
                         cell_types,
                         cluster_slot = NULL,
                         out_cluster_slot = "adobo_clusters",
                         symbol_switch = NULL,
                         new_idents = FALSE) {

  # Format objects.
  cluster_slot <- as.character(cluster_slot)
  gene_symbols <- as.vector(unlist(gene_symbols))
  cell_types <- as.vector(unlist(cell_types))

  # Prepare the Seurat object. Simply cast out the scaled data.
  # Note that if cluster_slot is null, idents are used.
  st_obj_temp <- st_obj
  clusters <- if (is.null(cluster_slot)) {
    Seurat::Idents(st_obj_temp)
  } else {
    unlist(st_obj_temp[[cluster_slot]])
  }
  st_obj_temp[["seurat_clusters"]] <- clusters
  st_obj_temp <- Seurat::DietSeurat(st_obj_temp, scale.data = FALSE,
                                    dimreducs = c("pca", "umap"))

  # Prep work.
  uuid <- uuid::UUIDgenerate()
  fn_base <- sprintf("%s/tmp/%s", impute_globals$base_dir, uuid)
  dir.create(fn_base, recursive = TRUE)
  fn_h5ad <- sprintf("%s/adobo.h5ad", fn_base)

  # Convert arguments to Python.
  ad_temp <- sceasy::convertFormat(st_obj_temp,
                                   from = "seurat",
                                   to = "anndata",
                                   outFile = fn_h5ad,
                                   assay = Seurat::DefaultAssay(st_obj_temp),
                                   drop_single_values = FALSE)
  marker_list <- data.table::data.table(symbol = gene_symbols,
                                        type = cell_types)
  adobo <- reticulate::import_from_path("adobo_wrapper",
    path = file.path(impute_globals$base_dir, "python"),
    convert = TRUE,
    delay_load = FALSE)

  # Perform the imputation and convert accordingly.
  # print(marker_list)
  # flush.console()
  res <- adobo$impute(ad_temp, marker_list, symbol_switch)
  fn_marker_tbl <- data.table::as.data.table(res)

  # Use the returned mapping to impute Seurat clusters.
  orig_clusters <- as.vector(unlist(clusters))
  st_obj[[out_cluster_slot]] <- dmSC::apply_mapper(orig_clusters,
                                                   mapper_table = fn_marker_tbl)
  if (new_idents) {
    Seurat::Idents(st_obj) <- st_obj[[out_cluster_slot]]
  }

  # Remove the temporary directory.
  unlink(fn_base, recursive = TRUE, force = TRUE)

  return(st_obj)
}

# Make the new environment and make the data folder if not present.
impute_globals <- new.env()
impute_globals$base_dir <- normalizePath(".")