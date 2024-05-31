#' Default mapping table creation
#'
#' Creates a default mapping table mapping between MGI symbols and ENSEMBL IDs.
#' Keys are in the first column, values are in the second. One-to-many and
#' one-to-one relations are tossed out to avoid problems with
#' preprocessing_function.
#'
#' Currently creates a one-to-one mapping table based on Mus Musculus genes and
#' ENSEMBL version 102.
#'
#' Note that we could allow one-to-one and many-to-one mappings, but NOT one-to-
#' many. However, the assumption toward the end of preprocessing_function would
#' be broken and we would need to write code that pools values with the same
#' ENSEMBL ID.
#'
#' @param mart The mart to be used for fetching ENSEMBL IDs and MGI symbols.
#'             Default is NULL. @param to_ensembl Whether ENSEMBL IDs should be
#' the keys (FALSE) or values (TRUE). Default is FALSE. @returns A list
#'                  containing the following: * A one-to-one mapping table based
#' on the description above ('mapper_table'), * A data.table containing the
#'          duplicated keys that were removed ('duplicate_keys'), * A data.table
#'            containing the duplicated values that were removed
#'          ('duplicate_values').
#' @export
create_default_mapper_table = function(mart=NULL, to_ensembl=FALSE) {
    # If no mart supplied, use the default mart.
    if (is.null(mart)) { mart = biomaRt::useEnsembl(biomart="genes",
                                                    dataset="mmusculus_gene_ensembl",
                                                    version=102) }
    
    # Make the mapping table. Ensure the direction of the mapping is
    # set based on to_ensembl.
    mapper_table = biomaRt::getBM(c('ensembl_gene_id', 'mgi_symbol'), mart=mart)
    mapper_table[mapper_table == ''] = NA
    mapper_table = mapper_table[,order(colnames(mapper_table))]
    if (to_ensembl) { mapper_table = mapper_table[, c(2,1)] }
    colnames(mapper_table) = paste0('V', 1:2)
    
    # Make mapping table one-to-one.
    mapper_deduplication_list = eliminate_mapper_all_duplicates(mapper_table)
    
    # Return the mapping table and duplicated keys/values.
    return(mapper_deduplication_list)
}

#' Eliminate duplicates in mapping table keys
#'
#' Eliminates duplicate values in the 'keys' (first column) of a mapper_table,
#' returning both the de-duplicated mapper_table and a data.frame consisting
#' of the removed duplicate keys (and the number of occurrences for each
#' duplicate key).
#' 
#' @param mapper_table The mapper_table whose keys are to be de-duplicated.
#' @returns A list containing the following:
#'          * A deduplicated mapping table based on the description above
#'            ('mapper_table'),
#'          * A data.table containing the duplicated keys that were removed
#'            ('duplicate_keys').
#' @export
eliminate_mapper_key_duplicates = function (mapper_table) {
    `%>%` = magrittr::`%>%`
    interestVar = colnames(mapper_table)[1]
    otherVar = colnames(mapper_table)[2]
    groupedDataInfo = mapper_table %>% dplyr::group_by(.data[[interestVar]]) %>%
                          dplyr::summarize(n_distinct = dplyr::n_distinct(.data[[otherVar]]))
    groupedDuplicateInfo = groupedDataInfo[groupedDataInfo[['n_distinct']] > 1,]
    duplicate_keys = as.vector(unlist(groupedDuplicateInfo[, interestVar]))
    return(list(mapper_table=data.table::as.data.table(mapper_table[!(mapper_table[[interestVar]] 
                                                                    %in% duplicate_keys),]),
                duplicate_keys=data.table::as.data.table(groupedDuplicateInfo)))
}

#' Eliminate duplicates in mapping table values
#'
#' Eliminates duplicate values in the 'values' (second column) of a mapper_table,
#' returning both the de-duplicated mapper_table and a data.frame consisting
#' of the removed duplicate values (and the number of occurrences for each
#' duplicate value).
#'
#' At the end of the day, analogous to
#' `eliminate_mapper_key_duplicates(rev(mapper_table))`.
#' 
#' @param mapper_table The mapper_table whose values are to be de-duplicated.
#' @returns A list containing the following:
#'          * A deduplicated mapping table based on the description above
#'            ('mapper_table'),
#'          * A data.table containing the duplicated values that were removed
#'            ('duplicate_values').
#' @export
eliminate_mapper_value_duplicates = function (mapper_table) {
    temp = eliminate_mapper_key_duplicates(rev(mapper_table))
    return(list(mapper_table=temp[['mapper_table']],
                duplicate_values=temp[['duplicate_keys']]))
}

#' Eliminate all duplicates in mapping table
#'
#' Eliminates duplicate values in all columns (both keys and values) of a
#' mapper_table, returning both the de-duplicated mapper_table and data.frames
#' consisting of the removed duplicate keys/values (and the number of occurrences
#' for each duplicate key/value).
#'
#' 
#' @param mapper_table The mapper_table whose values are to be de-duplicated.
#' @returns A list containing the following:
#'          * A deduplicated mapping table based on the description above
#'            ('mapper_table'),
#'          * A data.table containing the duplicated keys that were removed
#'            ('duplicate_keys').
#'          * A data.table containing the duplicated values that were removed
#'            ('duplicate_values').
#' @export
eliminate_mapper_all_duplicates = function (mapper_table) {
    temp1 = eliminate_mapper_key_duplicates(mapper_table)
    temp2 = eliminate_mapper_value_duplicates(mapper_table)
    duplicate_keys = temp1[['duplicate_keys']]
    duplicate_values = temp2[['duplicate_values']]
    return(list(mapper_table=data.table::as.data.table(dplyr::inner_join(temp1[['mapper_table']],
                                                                        temp2[['mapper_table']], 
                                          by=dplyr::join_by(V1, V2))),
                duplicate_keys=duplicate_keys,
                duplicate_values=duplicate_values))
}

#' Apply ENSEMBL naming to Seurat objects
#' 
#' Given a Seurat object with a specific name and assay, this function renames
#' the genes in the assay (assuming the gene names are MGI symbols) by
#' replacing them with their ENSEMBL IDs (ENSMUSG...).
#' 
#' @param st The Seurat object. Make sure that the default assay is the one that
#'           is to be processed.
#' @param mart If mapper_table is NULL, the mart to be used for fetching ENSEMBL IDs
#'             and MGI symbols. Default is NULL. If mapper_table is not NULL, this
#'             parameter will be ignored.
#' @param mapper_table The mapping table to use. Default is NULL; if a table is not
#'                    provided, then one will be created using the 
#'                    `create_default_mapper_table` function.
#' @returns A Seurat object with renamed symbols in a new assay with 'ENSEMBL'
#'          prepended to the old assay name. This new assay will be the
#'          default.
#' @export
apply_ensembl_seurat = function(st, mart=NULL, mapper_table=NULL) {
    # Get the assay information.
    ensembl_assay = Seurat::GetAssayData(st, slot="data")
    default_name = Seurat::DefaultAssay(st)

    # Apply the ENSEMBL row name mapping.
    variable_features = Seurat::GetAssay(st, assay=default_name)@var.features
    num_variable_features = length(variable_features)
    composite_names = c(rownames(ensembl_assay), variable_features)
    composite_names = apply_mapper(composite_names,
                                  to_ensembl=TRUE,
                                  mart=mart,
                                  mapper_table=mapper_table)

    # From the composite array, break it back up and eliminate nas.
    ensembl_rownames = composite_names[1:(length(composite_names)-num_variable_features)]
    ensembl_var_features = if (num_variable_features > 0) {
        composite_names[(length(composite_names)-num_variable_features):length(composite_names)]
    } else {
        character(0)
    }
    na_ensembl_rownames_mask = !is.na(ensembl_rownames)
    ensembl_rownames = ensembl_rownames[na_ensembl_rownames_mask]
    ensembl_var_features = ensembl_var_features[!is.na(ensembl_var_features)]

    # Complete the construction of the object.
    ensembl_assay = ensembl_assay[na_ensembl_rownames_mask,]
    rownames(ensembl_assay) = ensembl_rownames    
    ensembl_assay_name = paste('ENSEMBL', default_name, sep='')
    st[[ensembl_assay_name]] = Seurat::CreateAssayObject(data=ensembl_assay)
    Seurat::DefaultAssay(st) = ensembl_assay_name
    st[[ensembl_assay_name]]@var.features = ensembl_var_features

    return(st)
}

#' Standard preprocessing function for single-cell data
#' 
#' DEPRECATED: Recommended to just do preprocessing yourself so you have
#' more granular control over parameters. See the source code of this function
#' for an idea/starting point if still interested.
#' 
#' TO DO:
#' * Make mitochondrial gene identification case-insensitive.
#' * Manual assay specification in addition to just using the default assay.
#' * Preserve variable features in the ENSEMBL assay since we just change
#'   the gene names.
#' * It seems that this function does not work when multiple assays are present.
#'   Fix that.
#' 
#' Preprocess a Seurat single-cell dataset using the following techniques:
#' * Remove genes that are expressed in fewer than 5% of cells.
#' * Remove genes that are not expressed or have no variance in expression across
#'   all cells.
#' * Remove cells with UMI counts, feature counts, mitochondrial gene percent, and
#'   ribosomal gene percent outside 3 mean absolute deviations from the median.
#' * Remove cells with no ENSEMBL ID counterpart.
#' 
#' Note that in the Seurat dataset, mitochondrial genes must begin with 'mt-' and
#' ribosomal genes must begin with 'Rp'. The genes also cannot be ENSEMBL IDs to
#' begin with. The gene name requirements are case-sensitive. Some mitochondrial
#' genes may begin with CAPITAL 'MT-' based on the source; this will cause problems.
#' Also ensure that the assay of interest in the Seurat object is marked as the 
#' default assay, since this will be the assay that will be used.
#' 
#' @param st The Seurat single-cell file to be preprocessed. Read above for how to
#'           ensure this file is structured.
#' @param mapper_table The mapping table to use. Default is NULL; if a table is not
#'                    provided, then one will be created using the 
#'                    `create_default_mapper_table` function.
#' @param mart If mapper_table is NULL, the mart to be used for fetching ENSEMBL IDs
#'             and MGI symbols. Default is NULL. If mapper_table is not NULL, this
#'             parameter will be ignored.
#' 
#' @returns The preprocessed Seurat single-cell file. A new assay will be created with
#'          the key 'ENSEMBL' appended onto the end of the default assay of the input
#'          Seurat file.
#' 
#' @export
preprocessing_function = function (st, mart=NULL, mapper_table=NULL) {

    warning('DEPRECATED: Recommended to just do preprocessing yourself so you have
             more granular control over parameters. See the source code of this function
             for an idea/starting point if still interested.')

    # Default assay name.
    default_name = Seurat::DefaultAssay(st)
    
    # Identify mitochondrial and ribosomal features.
    st[['percent.mito']] = Seurat::PercentageFeatureSet(st, pattern = "^mt-")
    st[['percent.ribo']] = Seurat::PercentageFeatureSet(st, pattern = "^Rp")
    
    # Save pre-preprocessing seurat object.
    stPre = st
    
    # Remove all zero genes and zero-variance genes.
    assay = Seurat::GetAssayData(st)
    nonZeroRowInfo = dplyr::group_by(Matrix::summary(assay), i)
    nonZeroRowSummary = dplyr::summarize(nonZeroRowInfo, n = dplyr::n(),
                                         n_distinct = dplyr::n_distinct(j))

    # A: Keep the genes (rows) with nonzero expression in more than 5% of cells.
    # B: Keep the genes (rows) with nonzero variance in gene expression.
    # Think about it: you have a nonzero variance if...
    # B1) You have a row with more than one distinct (nonzero) value,
    # B2) You have a row with one distinct value, but zero values also exist.
    expressionPcts = nonZeroRowSummary[['n']] / ncol(assay)
    distinctValues = nonZeroRowSummary[['n_distinct']]
    zerosInRow = nonZeroRowSummary[['n']] != ncol(assay)
    validIndices = nonZeroRowSummary[['i']][(expressionPcts > 0.05) &
                                            (distinctValues > 1 | zerosInRow)]

    # Subset the dataset.
    genesToKeep = rownames(assay)[validIndices]
    st_temp = Seurat::DietSeurat(st, features=genesToKeep, assays=default_name)
    st[[default_name]] = st_temp[[defaultName]]
    
    # Put all objects in the Seurat data.table.
    ulMedian = function(x) { return(stats::median(unlist(x))) }
    ulMad = function(x) { return(stats::mad(unlist(x))) }

    # Calculate mean absolute deviation and medians for all values.
    medCounts = ulMedian(st[['nCount_RNA']])
    medFeatures = ulMedian(st[['nFeature_RNA']])
    medMito = ulMedian(st[['percent.mito']])
    medRibo = ulMedian(st[['percent.ribo']])
    madCounts = ulMad(st[['nCount_RNA']])
    madFeatures = ulMad(st[['nFeature_RNA']])
    madMito = ulMad(st[['percent.mito']])
    madRibo = ulMad(st[['percent.ribo']])
    st = BiocGenerics::subset(st, !(nCount_RNA > medCounts + 3 * madCounts | nCount_RNA < medCounts - 3 * madCounts |
                                    nFeature_RNA > medFeatures + 3 * madFeatures | nFeature_RNA < medFeatures - 3 * madFeatures |
                                    percent.mito > medMito + 3 * madMito | percent.mito < medMito - 3 * madMito |
                                    percent.ribo > medRibo + 3 * madRibo | percent.ribo < medRibo - 3 * madRibo))
    
    # Annotate genes with ensembl.
    # Make the mapper from MGI to Ensembl version 102 (GRCm38.p6)
    # Map gene names in the single-cell assay using the unprocessed data.
    # We assume that we have no MGI symbols that have mapped to the same ENSEMBL ID.
    st = apply_ensembl_seurat(st, mapper_table=mapper_table, mart=mart)

    # Normalize and scale the data as a convenience measure. Raw count data will always
    # exist in the counts slot.
    st = Seurat::ScaleData(Seurat::NormalizeData(st, verbose=FALSE),
                           verbose=FALSE)
    
    # Return!
    return(list(before=stPre, after=st))
}

#' Map using a mapper_table
#'
#' Given a list of gene symbols or ENSEMBL IDs and a mapping table describing
#' the appropriate relationships between them, converts the entries of `to_map`
#' into their corresponding values.
#'
#' Note that `to_map` must contain entries that fall within the keys of the
#' mapping table. These keys will then be turned into their corresponding values
#' within the mapping table. There may be duplicate keys or values in the mapping
#' table, but a mapping table containing a one-to-one mapping between keys and
#' values is recommended for compatibility with other preprocessing functions.
#'
#' Note that this can be used outside the specific context of ENSEMBL IDs and MGI
#' symbols.
#' 
#' UPDATE: The GetBM function may already carry out these functions. The
#' functionality of the method may be replaced in the future.
#'
#' @param to_map A vector containing values to map from.
#' @param mapper_table The mapping table to use. Default is NULL; if a table is not
#'                    provided, then one will be created using the 
#'                    `create_default_mapper_table` function.
#' @param mart If mapper_table is NULL, the mart to be used for fetching ENSEMBL IDs
#'             and MGI symbols. Default is NULL. If mapper_table is not NULL, this
#'             parameter will be ignored. Uses package biomaRt.
#' @param to_ensembl If mapper_table is NULL, whether ENSEMBL IDs should be the keys
#'                  (FALSE) or values (TRUE) of the created mapping table. Default
#'                  is FALSE. If mapper_table is not NULL, this parameter will be
#'                  ignored.
#' @returns A vector with mapped values (or 'NA' in the corresponding slot if the
#'          element in to_map is not in the mapping table).
#' @export
apply_mapper = function(to_map, to_ensembl=FALSE, mart=NULL, mapper_table=NULL) {
    # Prepare to_map.
    to_map = unlist(to_map)

    # Prepare the mapper table.
    if (is.null(mapper_table)) { 
        mapper_table = create_default_mapper_table(mart=mart,
                                                   to_ensembl=to_ensembl)[['mapper_table']]
    }

    # If column names aren't identical, reset them and issue a warning.
    # If there are too many columns, error.
    colnms <- colnames(mapper_table)
    ideal_colnms <- c('V1', 'V2')
    if (length(ideal_colnms) != length(colnms)) {
        stop(print('There must only be two columns in the mapper table.'))
    } else if (!identical(ideal_colnms, colnms)) {
        warning('Mapper table column names must be the vector c(\'V1\', \'V2\') exactly. 
                 Autoreassigning mapper table column names.')
        colnames(mapper_table) <- ideal_colnms
    }
    
    # Do the mapping.
    to_map = data.table::data.table(V1 = to_map)
    colnms = colnames(mapper_table)
    joined_table = dplyr::left_join(to_map, mapper_table, by=dplyr::join_by(V1))
    return(joined_table[['V2']])
}

#' Obtain a single-cell reference matrix compatible with CIBERSORTx from a Seurat object.
#' 
#' Note that this function is basically just calling GetAssayData and then reassigning
#' the column names of the resultant table.
#' 
#' @param st_obj The Seurat object to analyze.
#' @param cluster_colname The name of the column in st_obj[[]] that contains cluster
#'                        information. If NULL, information obtained from
#'                        Seurat::Idents(st_obj). Default is NULL.
#' @param assay The name of the assay to pull from. If NULL, gets the default assay
#'              from the Seurat object. Default is NULL.
#' 
#' @returns A data.frame representing the single-cell reference matrix. Careful, since
#'          these can get pretty large in size!
#' 
#' @export
seurat_to_sc_ref <- function(st_obj, cluster_colname = NULL, ...) {
    assay_data <- as.data.frame(Seurat::GetAssayData(st_obj, ...))
    if (is.null(cluster_colname)) { 
        colnames(assay_data) <- Seurat::Idents(st_obj)
    } else { 
        colnames(assay_data) <- unlist(st_obj[[cluster_colname]])
    }
    return(assay_data)
}

#' Impute cell types for Seurat object tgt_st with Seurat object ref_st.
#' 
#' @param tgt_st Target Seurat object to impute cell types against.
#' @param ref_st Reference Seurat object against which to impute cell types.
#' @param out_nm Column to create imputed clusters under in tgt_st.
#' @param tgt_clusters_nm Column in tgt_st containing clusters. If NULL, Seurat::Idents(tgt_st) is used.
#' @param ref_clusters_nm Column in ref_st containing clusters. If NULL, Seurat::Idents(ref_st) is used.
#' 
#' @returns Target Seurat object with extra column containing imputed clusters.
#' @export
impute_st <- function(tgt_st, ref_st, out_nm = "singler_clusters", tgt_clusters_nm = NULL, ref_clusters_nm = NULL) {
    tgt_st_copy <- tgt_st
    
    # Convert to single-cell experiments.
    tgt_sce <- Seurat::as.SingleCellExperiment(tgt_st)
    ref_sce <- Seurat::as.SingleCellExperiment(ref_st)
    
    # Get the labels (reference categories) and clusters (target categories).
    labels <- if(is.null(ref_clusters_nm))  { Seurat::Idents(ref_st) } else { as.vector(unlist(ref_st[[ref_clusters_nm]])) }
    clusters <- if(is.null(tgt_clusters_nm)) { Seurat::Idents(tgt_st) } else { as.vector(unlist(tgt_st[[tgt_clusters_nm]])) }

    # Impute based on cell type, fine-grained.
    pred_singler <- SingleR::SingleR(test = tgt_sce,
                                     ref = ref_sce,
                                     assay.type.test = "logcounts",
                                     assay.type.ref = "logcounts",
                                     labels = labels,
                                     clusters = clusters,
                                     de.method = "wilcox")
    
    # Map the old to the new clusters.

    mapper <- pred_singler$pruned.labels
    names(mapper) <- rownames(pred_singler)
    mapped_clusters <- as.vector(unlist(purrr::map(clusters, function (x) { mapper[x] })))
    tgt_st_copy[[out_nm]] <- mapped_clusters
    
    return(tgt_st_copy)
}