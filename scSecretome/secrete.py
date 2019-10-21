from scSecretome.annot import *
from scSecretome.sc import *
import scanpy as sc

# read secreted protein files
secreted_protein = read_HPA('protein_class_Predicted_secrete.tsv', return_all = True)
membrane_protein = read_HPA('protein_class_Predicted_membrane.tsv', return_all = True)

def merge_subset(sc_data, HPA_subset, name):
    '''
    Find with secreted/membrane protein that is variable in single cell data set. print top 20
    input:
    - sc_data: anndata from scanpy
    - HPA_subset: read_HPA(return_all = True)
    
    return: merged dataframe, top 20 variable secreted-membrane genes
    '''
    # call gene annotation and save in anndata.var
    sc_data.var.loc[sc_data.var.index.isin(HPA_subset.index), name] = True
    sc_data.var[name].fillna(False, inplace = True)
    
    # rank secreted protein by dispersion_norm
    
    secreted_merge = sc_data.var.loc[sc_data.var[name]].sort_values(by = 'dispersions_norm', ascending = False)
    top_twenty = HPA_subset.loc[secreted_merge.iloc[:20].index]
    print(top_twenty['Gene description'])
    
    return(top_twenty.index.tolist())
    # for each bin of mean expression, highly variable genes are selected.

import anndata
def recluster(sc_data, z_score_cutoff, name):
    '''
    re-cluster results using differential secreted protein
    - sc_data: original single cell data
    - merged: protein list with 'dispersions_norm'
    - z_score_cutoff: gene with dispersion norm above this cutoff will be selected
    '''
    # select the secreted proteins with high z-score
    
    subset = sc_data.var.loc[(sc_data.var['dispersions_norm'] >= z_score_cutoff) & sc_data.var[name]].index # select variable gene
               
    exp = sc_data[:, subset]
    # wipe out old annotations on genes
    secreted_exp = anndata.AnnData(exp.X, obs = exp.obs[['cluster', 'louvain']])
    
    # calculate the plotting
    sc.tl.pca(secreted_exp, svd_solver='arpack')
    sc.pp.neighbors(secreted_exp)
    
    # UMAP
    sc.tl.umap(secreted_exp)
    sc.pl.umap(secreted_exp, color = ['cluster', 'louvain']) # see if old structure is ruined
    
    # force directed
    sc.tl.draw_graph(secreted_exp)
    sc.pl.draw_graph(secreted_exp, color=['cluster', 'louvain'])
    
    # PAGA
    sc.tl.paga(secreted_exp)
    sc.pl.paga(secreted_exp, color=['cluster', 'louvain'])
    
    return(secreted_exp)