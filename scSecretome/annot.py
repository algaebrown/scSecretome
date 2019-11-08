import pandas as pd
def read_haemapedia(species = 'human'):
    '''
    read haemopedia results, return dataframe with Gene symbol + Lineage as columns
    species = 'human' or 'murine'
    
    Return: dataframe, columns = ['Gene Symbol', 'Lineage']
    possible values in Lineage 'Multi Potential Progenitor', 'Erythrocyte Lineage',
       'Megakaryocyte Lineage', 'Basophil Lineage', 'Eosinophil Lineage',
       'Neutrophil Lineage', 'Macrophage Lineage',
       'Dendritic Cell Lineage', 'B Cell Lineage', 'T Cell Lineage',
       'NK Cell Lineage']
    '''
    if species == 'human':
        df = pd.read_excel('~/haemopedia/lineage_specific_genes_in_human.xlsx')
        df = df[['Human Gene Symbol', 'Lineage']]
        
    else:
        df = pd.read_excel('~/haemopedia/lineage_specific_gene.xlsx')
        df = df[['GeneSymbol', 'Lineage']]
    
    # rename columns
    df.columns = ['Gene Symbol', 'Lineage']
    return df

def niche_specific_genes(filename = '~/haemopedia/Tikhonova2019_microenv_parsed.xlsx'):
    '''
    read cluster specific genes from Tikhonova 2019, nature, The bone marrow microenvironment at single-cell resolution
    the original excel is big, so we have a parsed version
    '''
    df = pd.read_excel(filename, header = 0)
    df = df[['Gene Symbol', 'Lineage']]
    return(df)

def housekeeping(file = '/home/hsher/tnbc_scrnaseq/data/housekeepers.txt'):
    '''
    retrun housekeeping gene from a https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/
    downloaded from  downloaded from https://github.com/Michorlab/tnbc_scrnaseq
    
    input: filename
    output: list of gene symbols that are house keeping (98 genes in total)
    
    '''
    df = pd.read_csv(file, header = None)
    return(df[0].values)


def read_HPA(tsv_file, return_all = False):
    '''
    extract gene symbols from human protein altas
    
    input: HPA tsv file at ~/human_protein_atlas
    output: list of gene symbols
    
    for example: 
    input: protein_class_Blood.tsv
    output: ~30 blood group antigen proteins's gene symbol
    '''
    base = '~/human_protein_atlas/'
    df = pd.read_csv(base+tsv_file, header = 0, index_col = 0, sep = '\t')
    if return_all:
        return(df[["Gene description", "RNA cell line specificity", "RNA blood lineage specificity score"]])
    return(df.index.values)

def human_mouse_homolog(filename = '~/HMD_HumanPhenotype.rpt'):
    df = pd.read_csv(filename, header = None, sep = '\t')
    
    # save only useful information
    df = df[[0,4]]
    df.columns = ['Human', 'Mouse']
    return(df)

def pellin_lineage():
    '''
    return lineage specific genes used by Pellin et. al (2019)
    See supplementary figure 4 caption
    '''
    genes = ['HLF', 'PLEK', 'HBB', 'CLC', 'ELANE', 'SAMHD1', 'MPO', 'IRF8', 'DNTT', 'CD34', 'CD164']
    lineage = ['P', 'Meg', 'E', 'BEMP','N', 'M', 'undifferentiated granulocytes', 'DC', 'Ly', 'Other', 'Other']
    df = pd.DataFrame(columns = ['Gene Symbol', 'Lineage'])
    df['Gene Symbol'] = genes
    df['Lineage'] = lineage
    return(df)