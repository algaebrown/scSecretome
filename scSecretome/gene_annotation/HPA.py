import pandas as pd
def read_HPA(tsv_file):
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
    return(df.index.values)