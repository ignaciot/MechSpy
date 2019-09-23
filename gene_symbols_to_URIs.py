import pandas as pd
import pickle
import mech_utils

gene_uris = []

G = mech_utils.get_graph()
entrez_df = pd.read_csv(mech_utils.ENTREZ_MAPPINGS, sep="\t", na_filter=False, \
                        usecols=[0, 2], names=['symbol', 'entrez_id'])
entrez2_df = pd.read_csv(mech_utils.ALIASES, sep="\t", na_filter=False, \
                        usecols=[10, 18], names=['alias', 'entrez_id'])

with open("./affy_genes/hgu133plus2_all_gene_symbols.txt", 'r') as fd:
    for line in fd.readlines():
        symbol = line[1:-2]     # get rid of double quotes from R output and newlines
        if symbol != '<NA>':
            entrez_id = None
            entrez_row = entrez_df[(entrez_df.symbol=='%s' % symbol)]
            if len(entrez_row) > 0:
                entrez_id = entrez_row.entrez_id.values[0]
            else:
                entrez2_row = entrez2_df[(entrez2_df.alias=='%s' % symbol)]
                if len(entrez2_row) > 0:
                    entrez_id = entrez2_row.entrez_id.values[0]
            if entrez_id:
                gene_uri = "<http://purl.uniprot.org/geneid/%s>" % entrez_id
                if gene_uri in G.nodes:
                    gene_uris.append(gene_uri)


pickle.dump(gene_uris, open("./affy_genes/hgu133plus2_all_URIs.pkl", 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

