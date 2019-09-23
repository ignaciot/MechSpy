import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)

import mech_labels

sns.set(font_scale=0.6)
#
#prefix = "textbook-Open-TG-Gates"
#prefix2 = "rest-Open-TG-Gates"
#
#data = np.load('%s_chem_cluster_data.npy' % prefix)
#data2 = np.load('%s_chem_cluster_data.npy' % prefix2)
#all_data = np.vstack((data, data2))
#labels = np.load('%s_chem_cluster_labels.npy' % prefix)
#labels2 = np.load('%s_chem_cluster_labels.npy' % prefix2)
#all_data = data

all_data = None
#for filename in glob.glob("*_chem_cluster_data.npy"):

for filename in ["all-experiments_chem_cluster_data.npy"]:
#for filename in ["textbook-Open-TG-Gates_chem_cluster_data.npy"]:
    if all_data is None:
        all_data = np.load(filename)
        all_labels = np.load(filename.replace("data", "labels"))
    else:
        all_data = np.vstack((all_data, np.load(filename)))
        all_labels = np.concatenate((all_labels, np.load(filename.replace("data", "labels"))))

all_data = all_data[~np.all(all_data == 0, axis=1)] # remove entries with all zeros
canonical_mech = []
#all_labels = np.concatenate((labels, labels2))
#all_labels = labels
short_labels = []
for label in all_labels:
    canonical_mech.append(label[label.find("(")+1:label.find(")")])
    chem = label.split(' (')[0]
    short_labels.append("%s (%s)" % (chem, ','.join(mech_labels.CHEMICAL_LABELS[chem])))

# all:
df = pd.DataFrame(all_data, columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M9', 'M10', 'M11', 'M12'])

# only close ones:
#df = pd.DataFrame(all_data[[1,2,5,6,7,8,9,10,11,12,13,14,15,16,19,23,24,25,26,30],:], columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M9', 'M10', 'M11', 'M12'])

# only mechs M1-M5:
#df = pd.DataFrame(all_data[:,:5], columns=['M1', 'M2', 'M3', 'M4', 'M5'])

lut = dict(zip((list(set(canonical_mech))), "rbgmky"))
canonical_mech = pd.Series(canonical_mech)
row_colors = canonical_mech.map(lut)
#for metric in ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "dice", "euclidean", "hamming", "jaccard", "kulsinski", "matching", "minkowski", "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean"]:
for metric in ["braycurtis"]:
    print(metric)
    #g = sns.clustermap(df, yticklabels=all_labels, row_colors=row_colors, figsize=(12,5), cmap="mako", 
    g = sns.clustermap(df, figsize=(15,5), yticklabels=short_labels, cmap="Purples", 
    #g = sns.clustermap(df, figsize=(12,5), yticklabels=[], cmap="Purples", 
    #                metric=metric)
                    metric=metric, z_score=0)
#                    col_cluster=False, metric=metric, z_score=0)
    #                   z_score=0)
                    #col_cluster=False, z_score=0)
    #plt.show()
    plt.tight_layout()
    plt.savefig('%s_dendrogram.png' % metric, dpi=300)
