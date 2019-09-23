import mech_labels
import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


possible_answers = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11']

top_acc = []
top2_acc = []
top3_acc = []
strat1_top_acc = []
strat1_top2_acc = []
strat1_top3_acc = []
strat2_top_acc = []
strat2_top2_acc = []
strat2_top3_acc = []
mixed_data = []

for i in range(1000):
    top_hits = 0
    top2_hits = 0
    top3_hits = 0

    top_predictions = 0
    close_predictions = 0
    somewhat_close_predictions = 0
    known_labels_count = 0

    # stratification: stats for chemicals with one possible label
    known_labels_count_strat_1 = 0
    top_predictions_strat_1 = 0
    close_predictions_strat_1 = 0
    somewhat_close_predictions_strat_1 = 0

    # stratification: stats for chemicals with two possible labels
    known_labels_count_strat_2 = 0
    top_predictions_strat_2 = 0
    close_predictions_strat_2 = 0
    somewhat_close_predictions_strat_2 = 0

    for chemical_name in mech_labels.CHEMICAL_LABELS.keys():
        random.shuffle(possible_answers)
        top = possible_answers[0]
        top2 = possible_answers[1]
        top3 = possible_answers[2]

        if len(mech_labels.CHEMICAL_LABELS[chemical_name]) > 0:
            known_labels_count += 1
            if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                known_labels_count_strat_1 += 1
            elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                known_labels_count_strat_2 += 1

            if top in mech_labels.CHEMICAL_LABELS[chemical_name]:
                top_predictions += 1
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    top_predictions_strat_1 += 1
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    top_predictions_strat_2 += 1
            elif top2 in mech_labels.CHEMICAL_LABELS[chemical_name]:
                close_predictions += 1
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    close_predictions_strat_1 += 1
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    close_predictions_strat_2 += 1
            elif top3 in mech_labels.CHEMICAL_LABELS[chemical_name]:
                somewhat_close_predictions += 1
                if len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 1:
                    somewhat_close_predictions_strat_1 += 1
                elif len(mech_labels.CHEMICAL_LABELS[chemical_name]) == 2:
                    somewhat_close_predictions_strat_2 += 1

    top_acc.append(float(top_predictions)/known_labels_count)
    top2_acc.append((float(top_predictions) + float(close_predictions))/known_labels_count)
    top3_acc.append((float(top_predictions) + float(close_predictions) + float(somewhat_close_predictions))/known_labels_count)
    mixed_data.append(['None\ntop', float(top_predictions)/known_labels_count])
    mixed_data.append(['None\ntop-2', (float(top_predictions) + float(close_predictions))/known_labels_count])
    mixed_data.append(['None\ntop-3', (float(top_predictions) + float(close_predictions) + float(somewhat_close_predictions))/known_labels_count])

    strat1_top_acc.append(float(top_predictions_strat_1)/known_labels_count_strat_1)
    strat1_top2_acc.append((float(top_predictions_strat_1) + float(close_predictions_strat_1))/known_labels_count_strat_1)
    strat1_top3_acc.append((float(top_predictions_strat_1) + float(close_predictions_strat_1) + float(somewhat_close_predictions_strat_1))/known_labels_count_strat_1)
    mixed_data.append(['1-mech\ntop', float(top_predictions_strat_1)/known_labels_count_strat_1])
    mixed_data.append(['1-mech\ntop-2', (float(top_predictions_strat_1) + float(close_predictions_strat_1))/known_labels_count_strat_1])
    mixed_data.append(['1-mech\ntop-3', (float(top_predictions_strat_1) + float(close_predictions_strat_1) + float(somewhat_close_predictions_strat_1))/known_labels_count_strat_1])

    strat2_top_acc.append(float(top_predictions_strat_2)/known_labels_count_strat_2)
    strat2_top2_acc.append((float(top_predictions_strat_2) + float(close_predictions_strat_2))/known_labels_count_strat_2)
    strat2_top3_acc.append((float(top_predictions_strat_2) + float(close_predictions_strat_2) + float(somewhat_close_predictions_strat_2))/known_labels_count_strat_2)
    mixed_data.append(['2-mechs\ntop', float(top_predictions_strat_1)/known_labels_count_strat_1])
    mixed_data.append(['2-mechs\ntop-2', (float(top_predictions_strat_1) + float(close_predictions_strat_1))/known_labels_count_strat_1])
    mixed_data.append(['2-mechs\ntop-3', (float(top_predictions_strat_1) + float(close_predictions_strat_1) + float(somewhat_close_predictions_strat_1))/known_labels_count_strat_1])

print("No stratification (all chemicals):")
print("Mean top-prediction precision: %.3f (std = %.3f)" % (np.mean(top_acc), np.std(top_acc)))
print("Mean top-2-predictions precision: %.3f (std = %.3f)" % (np.mean(top2_acc), np.std(top2_acc)))
print("Mean top-3-predictions precision: %.3f (std = %.3f)\n" % (np.mean(top3_acc), np.std(top3_acc)))

print("Stratification: chemicals with only one known mechanism:")
print("Mean top-prediction precision: %.3f (std = %.3f)" % (np.mean(strat1_top_acc), np.std(strat1_top_acc)))
print("Mean top-2-predictions precision: %.3f (std = %.3f)" % (np.mean(strat1_top2_acc), np.std(strat1_top2_acc)))
print("Mean top-3-predictions precision: %.3f (std = %.3f)\n" % (np.mean(strat1_top3_acc), np.std(strat1_top3_acc)))

print("Stratification: chemicals with only two known mechanisms:")
print("Mean top-prediction precision: %.3f (std = %.3f)" % (np.mean(strat2_top_acc), np.std(strat2_top_acc)))
print("Mean top-2-predictions precision: %.3f (std = %.3f)" % (np.mean(strat2_top2_acc), np.std(strat2_top2_acc)))
print("Mean top-3-predictions precision: %.3f (std = %.3f)\n" % (np.mean(strat2_top3_acc), np.std(strat2_top3_acc)))


data = pd.DataFrame(data=mixed_data, index=list(range(9000)), columns=['Stratification', 'Precision'])

rdata = [   ['All assays, top',     0.430],
            ['All assays, top-2',   0.647],
            ['All assays, top-3',   0.747],
            ['One mech, top',   0.376],
            ['One mech, top-2',   0.598],
            ['One mech, top-3',   0.632],
            ['Two mechs, top',   0.432],
            ['Two mechs, top-2',   0.648],
            ['Two mechs, top-3',   0.852],
        ]
real_data = pd.DataFrame(data=rdata, index=list(range(9)), columns=['Stratification', 'Precision'])

my_dpi = 300
plt.clf()
fig = plt.figure(figsize=(4500/my_dpi, 2500/my_dpi), dpi=my_dpi)
sns.set(font_scale=2, style="whitegrid")
sns.stripplot(x="Stratification", y="Precision", data=real_data, color='black', edgecolor='gray', size=15)
sns.violinplot(x="Stratification", y="Precision", data=data, palette="Purples")
#plt.axhline(y=0.430, xmin=0, xmax=0.15, color='r', lw=2)
#plt.axhline(y=0.647, xmin=0.15, xmax=0.25, color='r', lw=4)
#plt.axhline(y=0.747, xmin=0.25, xmax=0.35, color='r', lw=4)
#plt.axhline(y=0.376, xmin=0.35, xmax=0.45, color='r', lw=4)
#plt.axhline(y=0.598, xmin=0.45, xmax=0.55, color='r', lw=4)
#plt.axhline(y=0.632, xmin=0.55, xmax=0.65, color='r', lw=4)
#plt.axhline(y=0.432, xmin=0.65, xmax=0.75, color='r', lw=4)
#plt.axhline(y=0.648, xmin=0.75, xmax=0.85, color='r', lw=4)
#plt.axhline(y=0.852, xmin=0.85, xmax=0.95, color='r', lw=4)
plt.tight_layout()
plt.savefig('figures/mech_inf_violin_baseline.png')

