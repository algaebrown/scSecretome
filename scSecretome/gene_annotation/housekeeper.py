def housekeeping(file = '/home/hsher/tnbc_scrnaseq/data/housekeepers.txt'):
    '''
    retrun housekeeping gene from a https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/
    downloaded from  downloaded from https://github.com/Michorlab/tnbc_scrnaseq
    
    input: filename
    output: list of gene symbols that are house keeping (98 genes in total)
    
    '''
    df = pd.read_csv(file, header = None)
    return(df[0].values)
    