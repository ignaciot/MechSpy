import argparse
from graphviz import Digraph
import importlib
import multiprocessing
from natsort import natsorted, ns
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import textwrap

import mech_utils
import mechanisms_lib


labels_df = pd.read_csv(mech_utils.ONTO_LABELS, sep="\t", na_filter=False,
                        names=['concept', 'label'])

def generate_explanation(experiment_set_results):
    experiment_set = next(iter(experiment_set_results.keys()))
    expl_data = {}
    nr_time_points = len(exp_lib.EXPERIMENTS[experiment_set])
    G = mech_utils.get_non_closed_graph()

    if 'top_3_mech' not in experiment_set_results[experiment_set].keys():
        print("No predictions for %s, skipping." % experiment_set)
        return []

    for top_mech in experiment_set_results[experiment_set]['top_3_mech']:
        gene_paths_per_step = dict()
        direct_gene_paths_per_step = dict()
        all_significant_genes= {0: [], 1: [], 2: [], 3: [], 4: []}

        mech_steps = mechanisms_lib.MECHANISMS[top_mech]
        nr_mech_steps = len(mech_steps)
        found_paths = False
        time_point_index = 0
        for mech_step in mechanisms_lib.MECHANISMS[top_mech]:
            gene_paths_per_step[mech_step] = []
            direct_gene_paths_per_step[mech_step] = []
        print("Processing top mechanism %s for %s" % (top_mech, experiment_set))
        for result_key in natsorted(experiment_set_results[experiment_set].keys()):
            if result_key[-6:] == '_genes':
                # figure out which mechanism steps correspond to this time point (sequential weight of 1.0)
                mech_weight_bins = mech_utils.get_mech_step_weights(nr_mech_steps, 
                                                                    nr_time_points, 
                                                                    time_point_index)
                corresp_mech_step_indices = np.where(mech_weight_bins == 1.0)[0]

                # TODO: good loop to multiprocess:
                gene_nr = 0
                for (gene, gene_change) in list(set(experiment_set_results[experiment_set][result_key])):
                    gene_label = labels_df[(labels_df.concept==gene)].label.values[0]
                    if gene_nr < 5:
                        all_significant_genes[gene_nr].append(gene_label.split(" (")[0])
                    gene_nr += 1
                    for mech_step_idx in corresp_mech_step_indices:
                        mech_step = mechanisms_lib.MECHANISMS[top_mech][mech_step_idx]
                        aggregated_nodes = []
                        #print("Finding paths of significant gene %s to mechanism step %s" % (gene, mech_step))
                        try:
                            for path in nx.all_shortest_paths(G, source=gene, target=mech_step):
                                if len(path) == 2:
                                    direct_gene_paths_per_step[mech_step].append((gene_change, path))
                                if len(path) == 3 and 'geneid' in path[1]:
                                    aggregated_nodes.append(path[1])
                                    found_paths = True
                                elif len(path) <= mech_utils.MAX_NODES or any(['reactome' in x for x in path]):
                                    gene_paths_per_step[mech_step].append((gene_change, path))
                                    found_paths = True
                        except nx.exception.NetworkXNoPath as e:
                            #print("No paths")
                            pass
                        except nx.exception.NodeNotFound as e:
                            #print("No gene node for %s" % gene_node)
                            pass

                        if len(aggregated_nodes) > 1:
                            gene_paths_per_step[mech_step].append((gene_change, [gene, ';'.join(aggregated_nodes), mech_step]))
                        elif len(aggregated_nodes) > 0:
                            gene_paths_per_step[mech_step].append((gene_change, [gene, aggregated_nodes[0], mech_step]))

                time_point_index += 1

        output = ""
        if found_paths:
            output += "\nNarrative for mechanism prediction of %s for %s\n----------------------------------------------------------------------------------" % (top_mech, experiment_set)
        else:
            output += "\nNo particular explanation for the mechanistic prediction of %s for %s based on current domain knowledge.\n" % (top_mech, experiment_set)

        path_nr = 1
        explanation_graph = Digraph(engine="neato")
        explanation_graph.attr(overlap='false')
        explanation_graph.attr(overlap_scaling='-2')
        explanation_graph.attr(splines='curved')
        visited_nodes = []
        visited_edges = []
        graph_row = 0
        for mech_step in gene_paths_per_step.keys():
            prefix = ""
            last_concept_1 = ""
            for (gene_change, path) in gene_paths_per_step[mech_step]:
                destination = 1
                while destination < len(path):
                    try:
                        if '>;<' in path[destination]:
                            relation = "interacts with"
                        elif '>;<' in path[destination - 1]:
                            relation = "participate in"
                        else:
                            relation_concept = G.edges[path[destination - 1], path[destination]]['label']
                            relation = labels_df[(labels_df.concept==relation_concept)].label.values[0]
                    except KeyError as e:
                        print(path[destination - 1], path[destination])
                    except IndexError as e:
                        relation = relation_concept
                    try:
                        if ';' in path[destination - 1]:
                            concept_1 = ', '.join([labels_df[(labels_df.concept==some_node)].label.values[0] for some_node in path[destination - 1].split(';')])
                        else:
                            concept_1 = labels_df[(labels_df.concept==path[destination - 1])].label.values[0]
                    except IndexError as e:
                        concept_1 = path[destination - 1]
                    try:
                        if ';' in path[destination]:
                            concept_2 = ', '.join([labels_df[(labels_df.concept==some_node)].label.values[0] for some_node in path[destination].split(';')])
                        else:
                            concept_2 = labels_df[(labels_df.concept==path[destination])].label.values[0]
                    except IndexError as e:
                        concept_2 = path[destination]

                    concept_1_short = ''
                    if ';' in path[destination - 1]:
                        for current_label in [labels_df[(labels_df.concept==some_node)].label.values[0] for some_node in path[destination - 1].split(';')]:
                            if concept_1_short == '':
                                concept_1_short = current_label.split(" (")[0]
                            else:
                                concept_1_short = concept_1_short + '\n' + current_label.split(" (")[0]
                    else:
                        concept_1_short = concept_1.split(" (")[0]
                        if len(concept_1_short) > 8:
                            concept_1_short = '\n'.join(textwrap.wrap(concept_1_short, 20))

                    concept_2_short = ''
                    if ';' in path[destination]:
                        for current_label in [labels_df[(labels_df.concept==some_node)].label.values[0] for some_node in path[destination].split(';')]:
                            if concept_2_short == '':
                                concept_2_short = current_label.split(" (")[0]
                            else:
                                concept_2_short = concept_2_short + '\n' + current_label.split(" (")[0]
                    else:
                        concept_2_short = concept_2.split(" (")[0]
                        if len(concept_2_short) > 8:
                            concept_2_short = '\n'.join(textwrap.wrap(concept_2_short, 20))

                    if concept_1 == last_concept_1:
                        gene_change = "Also, %s" % gene_change
                    if relation == "is a" and concept_2[0] in ['a', 'e', 'i', 'o', 'u']:
                        relation = "is an"
                    max_len_1 = len(concept_1_short)
                    if max_len_1> 20:
                        max_len_1 = 20
                    max_len_2 = len(concept_2_short)
                    if max_len_2> 20:
                        max_len_2 = 20
                    if destination == 1:
                        output += ".\n%d. %s %s %s %s" % (path_nr, gene_change, concept_1, relation, concept_2.replace("\n", ''))
                        prefix = ", and"
                        if concept_1_short not in visited_nodes:
                            visited_nodes.append(concept_1_short)
                            if destination == 1:
                                explanation_graph.node(concept_1_short, concept_1_short, pos="0,%s!" % graph_row, style="filled", fillcolor="#d3d3d3", fontname="FreeSans", fontsize="12")
                            else:
                                explanation_graph.node(concept_1_short, concept_1_short, pos="0,%s!" % graph_row, style="filled", fillcolor="#f1f1f1", fontname="FreeSans", fontsize="12")
                        if concept_2_short not in visited_nodes:
                            visited_nodes.append(concept_2_short)
                            if path[destination] in mech_steps:
                                explanation_graph.node(concept_2_short, concept_2_short, pos="20,%s!" % graph_row, style="filled", fillcolor="#caaddb", fontname="FreeSans", fontsize="12")
                            else:
                                explanation_graph.node(concept_2_short, concept_2_short, pos="10,%s!" % graph_row, style="filled", fillcolor="#f1f1f1", fontname="FreeSans", fontsize="12")
                    else:
                        if ';' in path[destination - 1]:
                            if path[destination - 1].count(';') == 2:
                                output += "%s they both %s %s" % (prefix, relation, concept_2.replace("\n", ''))
                            else:
                                output += "%s they all %s %s" % (prefix, relation, concept_2.replace("\n", ''))
                        else:
                            output += "%s %s %s %s" % (prefix, concept_1_short, relation, concept_2.replace("\n", ''))
                        if concept_2_short not in visited_nodes:
                            visited_nodes.append(concept_2_short)
                            if path[destination] in mech_steps:
                                explanation_graph.node(concept_2_short, concept_2_short, pos="20,%s!" % graph_row, style="filled", fillcolor="#caaddb", fontname="FreeSans", fontsize="12")
                            else:
                                explanation_graph.node(concept_2_short, concept_2_short, pos="20,%s!" % graph_row, style="filled", fillcolor="#f1f1f1", fontname="FreeSans", fontsize="12")
                    destination += 1
                    last_concept_1 = concept_1
                    if "%s-%s" % (concept_1_short, concept_2_short) not in visited_edges:
                        explanation_graph.edge(concept_1_short, concept_2_short, label=relation, splines='ortho', fontname="FreeSans", fontsize="9")
                        visited_edges.append("%s-%s" % (concept_1_short, concept_2_short))
                graph_row = graph_row - 2 
                prefix = " Then, "
                path_nr += 1
        
        ordered_gene_list = []
        for i in range(5*time_point_index):
            ordered_gene_list = ordered_gene_list + [all_significant_genes[x].pop(0) for x in all_significant_genes.keys() if len(all_significant_genes[x]) > 0]
        output += ".\n\nGene knockout suggestions for further experiments at this dose and time points (in order of relevance): %s\n\n" % ", ".join(ordered_gene_list)
        print(output)
        with open('narratives/%s_%s_explanation.txt' % (experiment_set, top_mech), 'w') as new_mech_narrative_fd:
            new_mech_narrative_fd.write(output)

        explanation_graph.attr(label='Mechanistic explanation for %s of %s' % (top_mech, experiment_set), engine="neato")
        explanation_graph.render('narratives/%s_%s_explanation' % (experiment_set, top_mech), format="png")
        del explanation_graph

        time_point_index += 1        

    return [gene_paths_per_step, direct_gene_paths_per_step]
        


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-time-series', dest='input_ts', help='Name of the .py file where the input time series has been defined.', required=True)
    parser.add_argument('-t', '--nr-threads', dest='nr_threads', help='Numbed of concurrent threads to use for multiprocessing.', default=4, required=False)
    parser.add_argument('-o', '--output-dir', dest='output_dir', help='Output directory where the mechanistic narratives and graphical explanations will be created.', required=True)
    parser.add_argument('-k', '--keyword', dest='keyword', help='Search keyword to only produce an explanation for the given chemical[s] matching this term (e.g. "Doxorubicin").', required=False)
    args = parser.parse_args()
    to_import = args.input_ts.split('.py')[0]
    exp_lib = importlib.import_module(to_import)

    with open('%s_inference_data.pkl' % exp_lib.EXP_SET_ID, 'rb') as handle:
        results = pickle.load(handle)
        
        if args.keyword:
            filtered_results = []
            for result in results:
                if args.keyword in str(result.keys()):
                    filtered_results.append(result)
            results = filtered_results

        pool = multiprocessing.Pool(int(args.nr_threads))
        expl_results = pool.map(generate_explanation, results)
        pool.close()
        pool.join()

