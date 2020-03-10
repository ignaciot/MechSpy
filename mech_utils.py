import numpy as np
import networkx as nx
import os
import pandas as pd

MAX_NODES = 3
NON_CLOSED_MAX_NODES = 4
DIRECTED_NX_GRAPH = './knowledge_graphs/Directed_KG_triples_no-mesh.gpickle'
NON_CLOSED_DIRECTED_NX_GRAPH = './knowledge_graphs/non-closed_Directed_KG_triples_no-mesh.gpickle'
#NON_CLOSED_DIRECTED_NX_GRAPH = "../pheknowlator/PheKnowLator_full_NotClosed_NoOWLSemantics_Triples_Integers.gpickle"

ENTREZ_MAPPINGS = "./entrez/gene_symbol_to_entrez_id.tsv"
ALIASES = "./entrez/gene_with_protein_product.txt"

CUO_REL = "<http://purl.obolibrary.org/obo/RO_0002411>" # causally upstream of

PATHS_NP = 'numpy_dumps/mech2gene_paths.npy'
GENE_LABELS_NP = 'numpy_dumps/mech2gene_gene_labels.npy'
MECH_LABELS_NP = 'numpy_dumps/mech2gene_mech_labels.npy'

ONTO_LABELS = "labels/all_labels.tsv"
#ONTO_LABELS = "../pheknowlator/PheKnowLator_full_NotClosed_NoOWLSemantics_NodeLabels.txt"
TISSUE_SPECIFICITY = "./entrez/genes_active_in_tissue.csv"

PATH_COUNT_LIMIT = 10000

GRCH38_GENES = "./variants/refseq_hg38_genes.bed"
GRCH38_SNPS = "./variants/parsed_locii/pathogenic_snps.bed"
GIGGLE_INDICES = "./variants/locii_sort_b"

NODE_BLACKLIST = [  '<http://purl.obolibrary.org/obo/GO_0003674>',
                    '<http://purl.obolibrary.org/obo/GO_0008150>',
                    '<http://purl.obolibrary.org/obo/GO_0005575>',
                    '<http://purl.obolibrary.org/obo/GO_0000989>',
                    '<http://purl.obolibrary.org/obo/GO_0001134>',
                    '<http://purl.obolibrary.org/obo/GO_0006357>',
                    '<http://purl.obolibrary.org/obo/GO_0016591>',
                    '<http://purl.obolibrary.org/obo/GO_0090575>',
                    '<http://purl.obolibrary.org/obo/GO_0045892>',
                    '<http://purl.obolibrary.org/obo/GO_0045893 >',
                    '<http://purl.obolibrary.org/obo/GO_0061586 >',
                    '<http://purl.obolibrary.org/obo/GO_0010621 >',
                    '<http://purl.obolibrary.org/obo/GO_0001077>',
                    '<http://purl.obolibrary.org/obo/GO_0005515>',
                    '<http://purl.obolibrary.org/obo/GO_0008134>',
                    '<http://purl.obolibrary.org/obo/GO_0006366>',
                    '<https://reactome.org/content/detail/R-HSA-1643685>',
                    '<http://www.w3.org/2002/07/owl#Class>',
                    '<http://www.w3.org/2002/07/owl#NamedIndividual>',
                    '<http://purl.obolibrary.org/obo/go.owl>',
                    '<http://www.geneontology.org/formats/oboInOwl#SubsetProperty>',
                 ]

NODE_INDEX_FILENAME = 'embeddings/node_idx_for_node2vec.txt'
EMBEDDINGS_FILENAME = 'embeddings/node2vec_embeddings_no-mesh_q3_d32.txt'
EMBEDDINGS_NUMPY = "%s.npy" % EMBEDDINGS_FILENAME


# Initialize lookup structures
####################################################################

NODE_INDEX = dict()
with open(NODE_INDEX_FILENAME, 'r') as fd:
    for line in fd.readlines():
        chunks = line.split('\t')
        NODE_INDEX[chunks[1][:-1]] = int(chunks[0])

EMB = None
if os.path.isfile(EMBEDDINGS_NUMPY):
    EMB = np.load(EMBEDDINGS_NUMPY)
else:
    is_header = True
    emb_size = 0
    with open(EMBEDDINGS_FILENAME, 'r') as fd:
        for line in fd.readlines():
            chunks = line.split(' ')
            if is_header:
                is_header = False
                emb_size = int(chunks[1])
                EMB = np.zeros((int(chunks[0]), emb_size))
                next
            else:
                node_id = int(chunks[0])
                for i in list(range(1, emb_size+1)):
                    EMB[node_id, i-1] = float(chunks[i])
    np.save(EMBEDDINGS_NUMPY, EMB)



# Helper functions:
###################################################################

def get_mech_step_weights(nr_m, nr_t, current_t, strict=False):
    weights = np.zeros(nr_m)
    mech_splits = np.array_split(list(range(nr_m)), nr_t)
    penalty = 1.0 / nr_t
    split_nr = 0
    for split in mech_splits:
        for position in split:
            if strict:
                if current_t == split_nr:
                    weights[position] = 1
                else:
                    weights[position] = 0
            else:
                weights[position] = 1 - penalty * abs(current_t - split_nr) / 2
        split_nr += 1
    return weights

def get_graph():
    G = nx.read_gpickle(DIRECTED_NX_GRAPH)
    # Get rid of some root/hub nodes for path intersections we don't care about:
    for bl_node in NODE_BLACKLIST:
        if G.has_node(bl_node):
            G.remove_node(bl_node)
    return G

def get_non_closed_graph():
    G = nx.read_gpickle(NON_CLOSED_DIRECTED_NX_GRAPH)
    # Get rid of some root/hub nodes for path intersections we don't care about:
#    for bl_node in NODE_BLACKLIST:
#        if G.has_node(bl_node):
#            G.remove_node(bl_node)
    return G

def get_variants(gene_name):
    gene_df = pd.read_csv(GRCH38_GENES, sep="\t", na_filter=False, \
                          usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'name'])
    # account for multiple isoforms
    this_gene_df = gene_df[gene_df['name'] == gene_name]
    variants = []
    if len(this_gene_df) > 0:
        start = this_gene_df.start.min()
        end = this_gene_df.end.max()
        chrom = this_gene_df.chrom.unique()[0]
        snp_df = pd.read_csv(GRCH38_SNPS, sep="\t", na_filter=False, \
                            usecols=[0, 1, 2, 4], names=['chrom', 'start', 'end', 'rs'])
        variants = list(snp_df[(snp_df['chrom'] == chrom) & (snp_df['start'] >= start) & (snp_df['end'] <= end)].rs)

    return variants


def get_node_embedding(node):
    node_idx = NODE_INDEX[node]
    return EMB[node_idx]

