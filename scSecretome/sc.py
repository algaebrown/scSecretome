import numpy as np
import pandas as pd
# import annotation
from scSecretome.annot import *
import scanpy as sc
from sklearn import preprocessing

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

sc.settings.set_figure_params(dpi=80)

def read_raw_counts(tsv, ignore_lines = 1, transpose = True):
    '''
    read .tsv raw counts with columns = cells; row = genes
    works specially for nature 2019 data (ignore_line = 1), nature 2019 mice stroma (ignore_line = 0)
    
    input: .tsv file
    ingnore_lines: does not count header. ignore the second line and so on. (natuer 2019 human has library line)
    return anndata fitting scanpy
    '''
    # read raw counts data
    adata = sc.read(tsv, cache = True)
    
    # remove first row "Library", save as obs
    if ignore_lines > 0:
        for i in range(ignore_lines):
            name = adata.obs.index[i]
            print(name)
            adata.var[name] = adata.X[i]
        # transpose to make var = genes; obs = cells;
    if transpose:
        adata = adata[ignore_lines:, :].copy().T 
    
    return(adata)

def benchmarking(adata, species = 'human'):
    # number of UMI counts per cell
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    # number of genes per cell (non-zero)
    adata.obs['n_genes'] = np.count_nonzero(adata.X, axis = 1)
    
    # find mitochondiral genes and count
    if species == 'human':
        mito_genes = adata.var_names.str.startswith('MT-')
        if mito_genes.sum() == 0:
            # some dataset has different naming which is weird!
            mito_genes = adata.var_names.str.startswith('MT.')
    else:
        mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    
    # % housekeepers present
    housekeeper = housekeeping()
    if species != 'human':
        converter = human_mouse_homolog()
        housekeeper = converter.loc[converter['Human'].isin(housekeeper), 'Mouse'].tolist()
    # not all housekeepers are present
    hk = list(set(housekeeper).intersection(set(adata.var.index)))
    adata.obs['percent_housekeeper'] = np.count_nonzero(adata[:, hk].X, axis=1) / len(housekeeper)

def lineage_score(adata, ln, ln_name):
    '''
    calculate lineage score by 1. min-max scalar 2. sum 3. normalize by No. genes
    
    ln: lineage specific genes
    adata: single cell objects
    ln_name: lineage name, str
    
    output: lineage score in adata obs
    '''
    # min max scalar
    min_max_scaler = preprocessing.MinMaxScaler()
    scaled_lineage_matrix = min_max_scaler.fit_transform(adata[:, ln].X)
            
    # some lineage might have more genes, others don't. Therefore normalize by how many genes they have
    adata.obs[ln_name] =  np.sum(scaled_lineage_matrix, axis=1)/scaled_lineage_matrix.shape[1]
    
def lineage_calling(adata, species = 'human'):
    '''
    finding lineage-specific genes for normalized data
    sum expression value of a class of gene
    
    input: scanpy anndata
    output: None. information are in adata.obs
    '''
        
    # % blood lineage specific genes showing up
    for g in read_haemapedia(species = species).groupby('Lineage'):
        ln = list(set(g[1]['Gene Symbol'].values).intersection(set(adata.var.index)))
        
        # LINEAGE CALLING 
        if len(ln) > 0:
            # g[0]: name of lineage, g[1] Gene symbol dataframe
            
            # normalize by max and min of each lineage specific genes
            lineage_score(adata, ln, g[0])
        else:
            print('no blood lineage specific genes detected, check species')
    
    # lineage specific genes: mice stroma
    stroma = niche_specific_genes()
    if species == 'human':
        converter = human_mouse_homolog()
        converter.drop_duplicates(subset = ['Mouse'], inplace = True) # mouse can have multiple human homolog
        converter.set_index('Mouse', inplace = True)
        stroma['Gene Symbol'] = stroma['Gene Symbol'].map(converter.loc[stroma['Gene Symbol'], 'Human'])
    for g in stroma.groupby('Lineage'):
        ln = list(set(g[1]['Gene Symbol'].values).intersection(set(adata.var.index)))
        
        
        if len(ln) > 0:
            lineage_score(adata, ln, g[0])
        else:
            print('no stroma lineage specific genes detected, check species')

