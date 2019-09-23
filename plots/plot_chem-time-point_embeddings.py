import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection, preprocessing)
from adjustText import adjust_text
from sklearn.decomposition import PCA
import mechanisms_lib

#CHEM_MECH_SCORES = "chem_cluster_data.npy"
CHEM_MECH_SCORES = "textbook-Open-TG-Gates_chem_cluster_top-choice_data.npy"
#CHEM_MECH_LABELS = "chem_cluster_labels.npy"
CHEM_MECH_LABELS = "textbook-Open-TG-Gates_chem_cluster_top-choice_labels.npy"

RESULTS_DUMP_FILE = 'embeddings_for_each_time_point.pkl'

NODE_INDEX_FILENAME = 'embeddings/node_idx_for_node2vec.txt'
EMBEDDINGS_FILENAME = 'embeddings/node2vec_embeddings_no-mesh_q3_d32.txt'
EMBEDDINGS_NUMPY = "%s.npy" % EMBEDDINGS_FILENAME

X = np.load(CHEM_MECH_SCORES)
labels = np.load(CHEM_MECH_LABELS)
#X_norm = preprocessing.normalize(X)
sc = preprocessing.StandardScaler()
X_norm = sc.fit_transform(X)

#pca = PCA(n_components=2)
#pca.fit(X_norm)
#print(pca.explained_variance_ratio_)
#print(pca.singular_values_)
#plt.figure()
#ax = plt.subplot()
#ax.scatter(pca.components_[:, 0], pca.components_[:, 1])
#for i in range(len(labels)):
#    plt.text(pca.components_[0, i], pca.components_[1, i], str(labels[i]), fontdict={'size': 16})
#plt.show()


#plt.figure()
#perp = 100
#tsne = manifold.TSNE(n_components=2, init='random',random_state=0, perplexity=perp)
#X_tsne = tsne.fit_transform(X_norm)
#ax = plt.subplot()
#ax.scatter(X_tsne[:, 0], X_tsne[:, 1])
#texts = []
#for i in range(len(labels)):
#    plt.text(X_tsne[i, 0], X_tsne[i, 1], str(labels[i]), fontdict={'size': 16})
#    #texts.append(ax.text(X_tsne[i, 0], X_tsne[i, 1], str(y[i]), fontsize=18, color=label_color))
#    #texts.append(ax.text(X_tsne[i, 0], X_tsne[i, 1], str(labels[i]), fontsize=10))
##adjust_text(texts, force_points=1, expand_points=(2,2), expand_text=(2,2), arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))
##plt.title('t-SNE (perplexity=%d)\n %s' % (perp, filename))
##plt.savefig('barcode_cluster_%s_tsne.png' % filename, dpi=1200)
#plt.show()





results = pickle.load(open(RESULTS_DUMP_FILE, 'rb'))


y = []
embeddings = []
last_label = ''
current_embeddings = []
for exp in results.keys():
    current_label = ' '.join(exp.split('_')[:2])
    if last_label != '' and last_label != current_label:
        current_embeddings = np.array(current_embeddings)
        embeddings.append(np.mean(current_embeddings, axis=0))
        #y.append(exp.replace('_', ' '))
        y.append(current_label)
        current_embeddings = []
    else:
        current_embeddings.append(results[exp])
    last_label = current_label


node_index = dict()
with open(NODE_INDEX_FILENAME, 'r') as fd:
    for line in fd.readlines():
        chunks = line.split('\t')
        node_index[chunks[1][:-1]] = int(chunks[0])

emb = np.load(EMBEDDINGS_NUMPY)

mech_step_emb = dict()
for mechanism_label in mechanisms_lib.MECHANISMS.keys():
    mech_step_emb[mechanism_label] = []
    print("Processing mechanism %s" % mechanism_label)
    mech_nodes = mechanisms_lib.MECHANISMS[mechanism_label]
    mech_step_nr = 0
    current_embeddings = []
    for node in mech_nodes:
        mech_step_nr += 1
        node_idx = node_index[node]
        current_embeddings.append(emb[node_idx])
    current_embeddings = np.array(current_embeddings)
    embeddings.append(np.mean(current_embeddings, axis=0))
    y.append(mechanism_label)


plt.figure()
perp = 10
tsne = manifold.TSNE(n_components=2, init='random',random_state=0, perplexity=perp)
X_tsne = tsne.fit_transform(embeddings)
ax = plt.subplot()
ax.scatter(X_tsne[:, 0], X_tsne[:, 1])
texts = []
for i in range(len(y)):
    if y[i] in mechanisms_lib.MECHANISMS.keys():
        plt.text(X_tsne[i, 0], X_tsne[i, 1], str(y[i]), fontdict={'size': 16}, color='red')
    else:
        plt.text(X_tsne[i, 0], X_tsne[i, 1], str(y[i]), fontdict={'size': 16}, color='black')
    #texts.append(ax.text(X_tsne[i, 0], X_tsne[i, 1], str(y[i]), fontsize=18, color=label_color))
    #texts.append(ax.text(X_tsne[i, 0], X_tsne[i, 1], str(y[i]), fontsize=10))
#adjust_text(texts, force_points=1, expand_points=(2,2), expand_text=(2,2), arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))
#adjust_text(texts, arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))
#plt.title('t-SNE (perplexity=%d)\n %s' % (perp, filename))
#plt.savefig('barcode_cluster_%s_tsne.png' % filename, dpi=1200)
plt.show()

