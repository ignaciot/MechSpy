# -*- coding: utf-8 -*-

import networkx as nx
import pickle
import sys

RDF_EDGE_LIST = './knowledge_graphs/closed_kg_triples.txt'
NON_CLOSED_RDF_EDGE_LIST = './knowledge_graphs/not_closed_kg_triples.txt'
#DIRECTED_NX_GRAPH = './knowledge_graphs/Directed_KG_triples_no-mesh.gpickle'
NON_CLOSED_DIRECTED_NX_GRAPH = './knowledge_graphs/non-closed_Directed_KG_triples_no-mesh.gpickle'


if __name__=='__main__':
    KG = nx.DiGraph()
    line_nr = 0
    #with open(RDF_EDGE_LIST) as rdf_edge_list:
    with open(NON_CLOSED_RDF_EDGE_LIST) as rdf_edge_list:
        for line in rdf_edge_list:
            #if len(line.split(' ')) < 4:
            #    print(line)
            if len(line) > 1:
                subj, pred, obj, _ = line.split(' ')
                KG.add_edge(subj, obj, label=pred)
                line_nr += 1
            if line_nr % 500000 == 0:
                sys.stdout.write('.')
    #nx.write_gpickle(KG, NX_GRAPH)
    #nx.write_gpickle(KG, DIRECTED_NX_GRAPH)
    nx.write_gpickle(KG, NON_CLOSED_DIRECTED_NX_GRAPH)
    sys.stdout.write('\n\n')
