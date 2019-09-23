from progressbar import ProgressBar
pbar = ProgressBar()

KG_TRIPLES = 'knowledge_graphs/closed_kg_triples.txt'
AOPWIKI_TRIPLES = './aopwiki/new_aopwiki_triples.txt'
edges = []
node_names = []

new_node_idx = 1

def new_edges_from_triple_dump(triple_dump):
    with open(triple_dump, 'r') as fd:
        for line in pbar(fd.readlines()):
            chunks = line.split(' ')
            if len(chunks) > 1:
                if chunks[0] in node_names:
                    idx1 = node_names.index(chunks[0])
                else:
                    node_names.append(chunks[0])
                    idx1 = new_node_idx
                    new_node_idx += 1
                if chunks[2] in node_names:
                    idx2 = node_names.index(chunks[2])
                else:
                    node_names.append(chunks[2])
                    idx2 = new_node_idx
                    new_node_idx += 1
                edges.append((idx1, idx2))


new_edges_from_triple_dump(KG_TRIPLES)

new_edges_from_triple_dump(AOPWIKI_TRIPLES)

with open('edges_for_node2vec.txt', 'w') as fd:
    for pair in edges:
        fd.write("%s %s\n" % (pair[0], pair[1]))

with open('node_idx_for_node2vec.txt', 'w') as fd:
    node_idx = 1
    for node_name in node_names:
        fd.write("%s\t%s\n" % (node_idx, node_name))
        node_idx += 1

