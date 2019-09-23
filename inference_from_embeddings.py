# -*- coding: utf-8 -*-

import argparse
import importlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import random
from scipy.spatial import distance
from sklearn.utils import resample
from sklearn.metrics import confusion_matrix, f1_score, classification_report, roc_curve, auc, accuracy_score
from natsort import natsorted, ns
import pickle

import mech_utils
import mech_labels
import mechanisms_lib

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-time-series', dest='input_ts', help='Name of the .py file where the input time series has been defined.', required=True)
parser.add_argument('-g', '--max-nr-genes', dest='max_genes', help='Maximum number (N) of genes to use for inference (i.e. use no more than N genes from each time point).', required=False, default=100)
args = parser.parse_args()

to_import = args.input_ts.split('.py')[0]
exp_lib = importlib.import_module(to_import)

RESULTS_FOR_CLUSTERING = 'embeddings_for_each_time_point.pkl'
RESULTS_DUMP_FILE = '%s_inference_data.pkl' % exp_lib.EXP_SET_ID

P_VALUE_THRESHOLD = 0.05
BOOTSTRAP_RUNS = 1000


# TODO: CACHE THIS
mech_step_emb = dict()
# Gather embeddings for each mechanism's steps:
# {
#   'M1':   [
#               [ embedding vector ],
#               [ embedding vector ],
#               [ embedding vector ],
#           ],
#   'M2':   { ......
#
for mechanism_label in mechanisms_lib.MECHANISMS.keys():
    mech_step_emb[mechanism_label] = []
    print("Processing mechanism %s" % mechanism_label)
    mech_nodes = mechanisms_lib.MECHANISMS[mechanism_label]
    for node in mech_nodes:
        mech_step_emb[mechanism_label].append(mech_utils.get_node_embedding(node))


# Gather embeddings for each significant gene in each time point in each experiment
entrez_df = pd.read_csv(mech_utils.ENTREZ_MAPPINGS, sep="\t", na_filter=False, \
                        usecols=[0, 2], names=['symbol', 'entrez_id'])
entrez2_df = pd.read_csv(mech_utils.ALIASES, sep="\t", na_filter=False, \
                        usecols=[10, 18], names=['alias', 'entrez_id'])

results = []
results_dump = dict()

for experiment_set in exp_lib.EXPERIMENTS.keys():
    print('Processing %s experiments' % experiment_set)
    
    enrichment_results = dict()
    enrichment_results[experiment_set] = dict()
    for mechanism_label in mechanisms_lib.MECHANISMS.keys():
        enrichment_results[experiment_set][mechanism_label] = dict()


    nr_time_points = len(exp_lib.EXPERIMENTS[experiment_set])
    time_point_index = 0
    # TEST: pseudo-random shuffling of the experimental time points
    #       (shuffle them but ensure they end up out of order)
    #shuffled_time_points = random.sample(exp_lib.EXPERIMENTS[experiment_set].keys(), len(exp_lib.EXPERIMENTS[experiment_set].keys()))
    #while natsorted(exp_lib.EXPERIMENTS[experiment_set].keys()) == shuffled_time_points:
    #    shuffled_time_points = random.sample(exp_lib.EXPERIMENTS[experiment_set].keys(), len(exp_lib.EXPERIMENTS[experiment_set].keys()))
    #for exp_time_point in shuffled_time_points:
    for exp_time_point in natsorted(exp_lib.EXPERIMENTS[experiment_set].keys()):
        gene_emb = []
        log = ""
        log = log + "   Processing experiment time point %s\n" % exp_time_point
        exp_df = pd.read_csv(exp_lib.EXPERIMENTS[experiment_set][exp_time_point],
                             sep="\s+", header=1, na_filter=False, \
                             usecols=[1, 2, 3],
                             names=['symbol', 'fold_change', 'p_value'],
                             dtype={'symbol':'str', 'fold_change':'float',
                                    'p_value':'float'})

        gene_count = 0
        enrichment_results[experiment_set]["%s_genes" % exp_time_point] = []
        for gene in exp_df.itertuples():
            if gene.p_value < P_VALUE_THRESHOLD:
                if gene.symbol != '<NA>':
                    entrez_id = None
                    entrez_row = entrez_df[(entrez_df.symbol=='%s' % gene.symbol)]
                    if len(entrez_row) > 0:
                        entrez_id = entrez_row.entrez_id.values[0]
                    else:
                        entrez2_row = entrez2_df[(entrez2_df.alias=='%s' % gene.symbol)]
                        if len(entrez2_row) > 0:
                            entrez_id = entrez2_row.entrez_id.values[0]
                    if entrez_id:
                        gene_node = "<http://purl.uniprot.org/geneid/%s>" % entrez_id
                    else:
                        continue

                    if gene_node in mech_utils.NODE_INDEX:
                        gene_emb.append(mech_utils.get_node_embedding(gene_node))
                        gene_change = "Upregulated"
                        if gene.fold_change < 0:
                            gene_change = "Downregulated"
                        enrichment_results[experiment_set]["%s_genes" % exp_time_point].append((gene_node, gene_change))

                        gene_count += 1
                        if gene_count == args.max_genes:
                            break
                    #else:
                    #    print("Did not find a node for gene %s (%s)" % (gene_node, gene.symbol))

        print("Gene count for %s: %d" % (exp_time_point, gene_count))
        if gene_count == 0:
            # skip this assay since there were no useful genes
            print("Skipping assay %s since there were no useful, significant gene changes" % exp_time_point)
            continue

        gene_emb = np.array(gene_emb)
        # Average all gene embeddings, calculate the distance to each mechanism step
        gene_emb_avg = np.mean(gene_emb, axis=0)
        
        results_dump[exp_time_point] = gene_emb_avg

        for mechanism_label in mech_step_emb.keys():
            mech_nodes = mech_step_emb[mechanism_label]
            nr_mech_nodes = len(mech_nodes)
            enriched_mech_nodes = np.zeros(nr_mech_nodes)
            mech_weight_bins = mech_utils.get_mech_step_weights(nr_mech_nodes, nr_time_points, time_point_index)
            i = 0
            for node in mech_nodes:
                mech_step_score = distance.cosine(gene_emb_avg, node)
                enriched_mech_nodes[i] = (1-mech_step_score) * mech_weight_bins[i]
                i += 1

            enrichment_results[experiment_set][mechanism_label][exp_time_point] = enriched_mech_nodes
        print(log)
        time_point_index += 1
    results.append(enrichment_results)

pickle.dump(results_dump, open(RESULTS_FOR_CLUSTERING, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

chem_cluster = []
top_choice = []
top_choice_labels = []
chem_cluster_labels = []
print("Mechanistic predictions for %s:" % exp_lib.EXP_SET_ID)
G = mech_utils.get_graph()
all_genes = [x for x in G.nodes if "<http://purl.uniprot.org/geneid/" in x]
hgu133plus2_genes = pickle.load(open("./affy_genes/hgu133plus2_all_URIs.pkl", 'rb'))
predictions = []
i = 0
known_labels_count = 0
top_predictions = 0
close_predictions = 0
somewhat_close_predictions = 0

# stratification: stats for chemicals with one possible label
known_labels_count_strat_1 = 0
top_predictions_strat_1 = 0
close_predictions_strat_1 = 0
somewhat_close_predictions_strat_1 = 0
y_test = []
y_pred = []

# stratification: stats for chemicals with two possible labels
known_labels_count_strat_2 = 0
top_predictions_strat_2 = 0
close_predictions_strat_2 = 0
somewhat_close_predictions_strat_2 = 0

for result in results:
    experiment_set = next(iter(result.keys()))
    print("\n%s" % experiment_set)
    scores = []
    labels = []
    p_vals = []
    z_scores = []
    chem_features = np.array([])
    for mechanism_label in result[experiment_set].keys():
        if mechanism_label[-6:] != '_genes':
            mechanism_score = np.zeros(len(mechanisms_lib.MECHANISMS[mechanism_label]))
            mech_gene_counts = dict()
            for exp_time_point in result[experiment_set][mechanism_label].keys():
                mechanism_score = np.vstack((mechanism_score, result[experiment_set][mechanism_label][exp_time_point]))
                mech_gene_counts[exp_time_point] = len(result[experiment_set]["%s_genes" % exp_time_point])
            mechanism_score = mechanism_score.max(axis=0)
            final_mechanism_score = np.median(mechanism_score)
            # Calculate the significance of this score
            empirical_score_values = []
            mech_nodes = mech_step_emb[mechanism_label]
            nr_mech_nodes = len(mech_nodes)
            nr_time_points = len(mech_gene_counts.keys())
            for _ in range(BOOTSTRAP_RUNS):
                # draw the same number of genes for each time point at random, and calculate the
                # enrichment score for this mechanism from the random genes
                time_point_index = 0
                random_mechanism_score = np.zeros(nr_mech_nodes)
                for exp_time_point in natsorted(mech_gene_counts.keys()):
                    random_genes = resample(hgu133plus2_genes, n_samples=mech_gene_counts[exp_time_point], replace=True) 
                    gene_emb = []
                    for gene_node in random_genes:
                        gene_emb.append(mech_utils.get_node_embedding(gene_node))
                    gene_emb = np.array(gene_emb)
                    gene_emb_avg = np.mean(gene_emb, axis=0)
                    enriched_mech_nodes = np.zeros(nr_mech_nodes)
                    mech_weight_bins = mech_utils.get_mech_step_weights(nr_mech_nodes, nr_time_points, time_point_index)
                    j = 0
                    for node in mech_nodes:
                        mech_step_score = distance.cosine(gene_emb_avg, node)
                        enriched_mech_nodes[j] = (1-mech_step_score) * mech_weight_bins[j]
                        j += 1

                    random_mechanism_score = np.vstack((random_mechanism_score, enriched_mech_nodes))
                    time_point_index += 1
                random_mechanism_score = random_mechanism_score.max(axis=0)
                empirical_score_values.append(np.average(random_mechanism_score))
            # We now have an empirical distribution of enrichment scores for this mechanism. How unlikely
            # is it to get the score we obtained?
            #print("mean: %.6f  -  stdev: %.6f  -  z-score=%.3f" % (np.mean(empirical_score_values), np.std(empirical_score_values), (final_mechanism_score - np.mean(empirical_score_values))/np.std(empirical_score_values)))
            #plt.clf()
            #plt.title(mechanism_label)
            #plt.hist(empirical_score_values, color='b')
            #plt.axvline(final_mechanism_score, color='r')
            #plt.show()
            median_random_score = np.median(empirical_score_values)
            # In each case, the null hypothesis is that our score is not significantly different than
            # the median score on our random bootstrap simulations
            if final_mechanism_score < median_random_score:
                p_val = 1.0
            else:
                p_val = float(sum(empirical_score_values >= final_mechanism_score) + 1) / (BOOTSTRAP_RUNS + 1)
            z_score = (final_mechanism_score - np.mean(empirical_score_values))/np.std(empirical_score_values)
            print("%s score: %.3f (p-val=%.2E , z=%.2f)\t%s" % (mechanism_label, final_mechanism_score, p_val, z_score, str(mechanism_score)))
            if p_val <= 0.05:
                labels.append(mechanism_label)
                p_vals.append(p_val)
                z_scores.append(z_score)
                scores.append(final_mechanism_score)
            chem_features = np.append(chem_features, final_mechanism_score)

    chemical_name = experiment_set.split(" (")[0]
    if len(labels) == 0:
        # just skip this one if none of the calls is statistically significant
        print("No statistically significant mechanism could be determined.")
#        print("PAPER:%s & %s & %s & %s & N\/A & N\/A & N\/A \\\\" % (
#                                                    experiment_set.split(" (")[0],
#                                                    experiment_set.split(" (")[1].split(')')[0].replace("-", "").replace("uM", '$\mu$M').replace("micromolar", '$\mu$M').replace("millimolar", 'mM'),
#                                                    experiment_set.split(" [")[1].split(']')[0].replace("-", " "),
#                                                    ', '.join(mech_labels.CHEMICAL_LABELS[chemical_name]),
#                                                    ))
    else:
        labels = np.array(labels)
        p_vals = np.array(p_vals)
        chem_cluster.append(chem_features)
        chem_cluster_labels.append(experiment_set)
        predictions.append(labels[scores.index(max(scores))])
        sorted_scores = np.sort(scores)
        sorted_score_indices = np.argsort(scores)[::-1]
        chemical_with_known_label = chemical_name in mech_labels.CHEMICAL_LABELS.keys()
        if chemical_with_known_label and len(mech_labels.CHEMICAL_LABELS[chemical_name]) > 0 and sum(scores) > 0:
            print("Known mechanisms: %s" % ', '.join(mech_labels.CHEMICAL_LABELS[chemical_name]))
            known_labels_count += 1
            if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                known_labels_count_strat_1 += 1
                y_test.append(mech_labels.CHEMICAL_LABELS[chemical_name][0])
            elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                known_labels_count_strat_2 += 1
            if labels[sorted_score_indices][0] in mech_labels.CHEMICAL_LABELS[chemical_name]:
                prediction = 'RIGHT ON'
                top_predictions += 1
                top_choice.append(chem_features)
                top_choice_labels.append(experiment_set)
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    top_predictions_strat_1 += 1
                    y_pred.append(labels[sorted_score_indices][0])
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    top_predictions_strat_2 += 1
            elif len(labels) > 1 and labels[sorted_score_indices][1] in mech_labels.CHEMICAL_LABELS[chemical_name]:
                prediction = 'CLOSE'
                close_predictions += 1
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    close_predictions_strat_1 += 1
                    y_pred.append(labels[sorted_score_indices][1])
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    close_predictions_strat_2 += 1
            elif len(labels) > 2 and labels[sorted_score_indices][2] in mech_labels.CHEMICAL_LABELS[chemical_name]:
                prediction = 'SOMEWHAT CLOSE'
                somewhat_close_predictions += 1
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    somewhat_close_predictions_strat_1 += 1
                    y_pred.append(labels[sorted_score_indices][2])
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    somewhat_close_predictions_strat_2 += 1
            else:
                prediction = 'WRONG'
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    y_pred.append(labels[sorted_score_indices][0])
        else:
            prediction = 'UNKNOWN'
        results[i][experiment_set]['top_3_mech'] = labels[sorted_score_indices][:3]

        print("Most likely mechanism: %s (score=%.3f , p=%.2E) - %s, then %s" % 
                (labels[scores.index(max(scores))], 
                max(scores),
                p_vals[scores.index(max(scores))], 
                prediction,
                ', '.join(["\n                       %s (score=%.3f , p=%.2E)" % (labels[x], scores[x], p_vals[x]) for x in sorted_score_indices[1:]])))

        first = labels[sorted_score_indices][0]
        if chemical_with_known_label and first in mech_labels.CHEMICAL_LABELS[chemical_name]:
            first = "\\textbf{%s}" % labels[sorted_score_indices][0]
        first_pval = p_vals[sorted_score_indices][0]
        first = "%s (%.2E)" % (first, first_pval)
        second = "N\/A"
        if len(labels) > 1:
            second = labels[sorted_score_indices][1]
            if chemical_with_known_label and second in mech_labels.CHEMICAL_LABELS[chemical_name]:
                second = "\\textbf{%s}" % labels[sorted_score_indices][1]
            second_pval = p_vals[sorted_score_indices][1]
            second = "%s (%.2E)" % (second, second_pval)
        third = "N\/A"
        if len(labels) > 2:
            third = labels[sorted_score_indices][2]
            if chemical_with_known_label and third in mech_labels.CHEMICAL_LABELS[chemical_name]:
                third = "\\textbf{%s}" % labels[sorted_score_indices][2]
            third_pval = p_vals[sorted_score_indices][2]
            third = "%s (%.2E)" % (third, third_pval)

        known_label = "UNKNOWN"
        if chemical_with_known_label:
            known_label = ', '.join(mech_labels.CHEMICAL_LABELS[chemical_name])
#        print("PAPER:%s & %s & %s & %s & %s & %s & %s \\\\" % (
#                                                    experiment_set.split(" (")[0],
#                                                    experiment_set.split(" (")[1].split(')')[0].replace("-", "").replace("uM", '$\mu$M').replace("micromolar", '$\mu$M').replace("nanomolar", 'nM').replace("millimolar", 'mM'),
#                                                    experiment_set.split(" [")[1].split(']')[0].replace("-", " "),
#                                                    known_label,
#                                                    first,
#                                                    second,
#                                                    third,
#                                                    ))


    i += 1

pickle.dump(results, open(RESULTS_DUMP_FILE, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

if known_labels_count > 0:
    print("-------------------------------------------------------------------")
    print("Stratification: all assays")
    print("Accuracy in top mechanism (all samples with 1+ prediction): %.3f" % (float(top_predictions)/known_labels_count))
    print("\nAccuracy in top 2 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions) + float(close_predictions))/known_labels_count))
    print("\nAccuracy in top 3 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions) + float(close_predictions) + float(somewhat_close_predictions))/known_labels_count))
    print("\n-------------------------------------------------------------------")
    print("Stratification: chemicals with only one known mechanism of toxicity.")
    print("\nAccuracy in top mechanism (all samples with 1+ prediction): %.3f" % (float(top_predictions_strat_1)/known_labels_count_strat_1))
    print("\nAccuracy in top 2 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions_strat_1) + float(close_predictions_strat_1))/known_labels_count_strat_1))
    print("\nAccuracy in top 3 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions_strat_1) + float(close_predictions_strat_1) + float(somewhat_close_predictions_strat_1))/known_labels_count_strat_1))
    print("\n-------------------------------------------------------------------")
    print("Stratification: chemicals with only two known mechanisms of toxicity.")
    print("\nAccuracy in top mechanism (all samples with 1+ prediction): %.3f" % (float(top_predictions_strat_2)/known_labels_count_strat_2))
    print("\nAccuracy in top 2 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions_strat_2) + float(close_predictions_strat_2))/known_labels_count_strat_2))
    print("\nAccuracy in top 3 mechanisms (all samples with 1+ prediction): %.3f" % ((float(top_predictions_strat_2) + float(close_predictions_strat_2) + float(somewhat_close_predictions_strat_2))/known_labels_count_strat_2))
    print("\n-------------------------------------------------------------------")

np.save("%s_chem_cluster_data.npy" % exp_lib.EXP_SET_ID, np.array(chem_cluster))
np.save("%s_chem_cluster_labels.npy" % exp_lib.EXP_SET_ID, np.array(chem_cluster_labels))
np.save("%s_chem_cluster_top-choice_data.npy" % exp_lib.EXP_SET_ID, np.array(top_choice))
np.save("%s_chem_cluster_top-choice_labels.npy" % exp_lib.EXP_SET_ID, np.array(top_choice_labels))

print("\nUsing mechanisms described as:")
for mech in mechanisms_lib.MECHANISMS.keys():
    print("%s:" % mech)
    for step in mechanisms_lib.MECHANISMS[mech]:
        print(step)


