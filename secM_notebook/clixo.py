import networkx as nx
import pandas as pd
import ddot
from ddot import Ontology

net = nx.read_edgelist('/home/hsher/secM/string_subnetwork_secM') 
df = nx.convert_matrix.to_pandas_dataframe(net, weight = 'combined_score')

print(df.sum(axis = 1))
print(df.shape)

ont = Ontology.run_clixo(df, alpha = 0.2, beta = 0.9, output = '/home/hsher/secM/secM_ontology', square = True, square_names = df.index.tolist())

# put to ndex
ndex_server = 'http://public.ndexbio.org'
ndex_user, ndex_pass = 'b101102109', '5mY87$747l'

url, _ = ont.to_ndex(ndex_server=ndex_server, ndex_user=ndex_user, ndex_pass=ndex_pass)
print(url)
