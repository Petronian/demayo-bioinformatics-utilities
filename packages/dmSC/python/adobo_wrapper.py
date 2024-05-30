# Script to impute clusters using the panglao method.
# Need to ensure adobo is installed.

import adobo as ad
import numpy as np
import scanpy as sc
import pandas as pd
import os
import sys

# Symbol switch for 1d arrays.
def _symbol_switch_arr(arr, species):
    v = (os.path.dirname(ad.IO.__file__), species)
    if species == 'human':
        fn = '%s/data/%s.gencode_v32.genes.txt' % v
    elif species == 'mouse':
        fn = '%s/data/%s.gencode_v23.genes.txt' % v
        
    i = pd.Series(np.asarray(arr).flatten())
    if (np.sum(i.str.match('^ENS.*[0-9]{3,}\\.[0-9]+$'))/len(i)) > 0.99:
        i = pd.Series(np.array(i.str.split('\\.', expand=True).to_list())[:, 0])
    
    gs = pd.read_csv(fn, sep='\t', header=None)
    gs.index = gs.loc[:, 0]
    gs = gs[gs.index.isin(i)]
    missing = i[np.logical_not(i.isin(gs.index))]
    gs = pd.concat([gs, pd.DataFrame({0: missing, 1: ['NA']*len(missing)})])
    gs.index = gs.iloc[:, 0]
    gs = gs.reindex(i)
    return (gs[[1]].values).flatten()

def impute(sc_obj, marker_table, symbol_switch_key = None):
  # Variable names.
  marker_symbol_nm = "symbol"
  # marker_celltype_nm = "type"

  # Note that the data are normalized but not scaled!
  sc_ad = ad.data.dataset(raw_mat=pd.DataFrame(sc_obj.X.T.todense(),
                                               index=sc_obj.var_names,
                                               columns=sc_obj.obs_names),
                          sparse=True, desc='Adobo prediction')
  marker_table[marker_symbol_nm] = marker_table[marker_symbol_nm].str.upper()

  # Load the sc_obj data into the sc_ad object.
  sc_ad.norm_data = {'seurat-external': {'data': sc_ad.count_data,
                     'log': True,
                     'method': 'seurat-external',
                     'clusters': {'seurat-external-clusters': {'membership': sc_obj.obs[['seurat_clusters']].squeeze()}}}}
  
  # If symbol switch, do the ENSEMBL symbol switch:
  if symbol_switch_key is not None:
    # Capitalize all markers in marker table since adobo uses this for symbol
    # switching. If not done, adobo fails to recognize the genes here.
    marker_table[marker_symbol_nm] = marker_table[marker_symbol_nm].str.upper()
    ad.preproc.symbol_switch(sc_ad, species = symbol_switch_key.lower())
    marker_table[marker_symbol_nm] = _symbol_switch_arr(
      marker_table[[marker_symbol_nm]], symbol_switch_key.lower())

  # Capitalize one more time.
  marker_table[marker_symbol_nm] = marker_table[marker_symbol_nm].str.upper()
  marker_table.to_csv("./temp.csv")

  # Do the cell type imputation.
  ad.bio.cell_type_predict(sc_ad, cell_type_markers=marker_table)
  
  # Create the prediction mapper table and then write it to disk.
  predictionTable = sc_ad.norm_data['seurat-external']['clusters']['seurat-external-clusters']['cell_type_prediction']
  mapperTable = pd.DataFrame(predictionTable['cell type'])
  mapperTable.index.name = 'V1'
  mapperTable.columns = ['V2']
  mapperTable = mapperTable.reset_index().loc[:, ['V1', 'V2']].astype(str)
  return mapperTable