import glob
import os

EXP_SET_ID = "diXa-004-lung"
EXP_CELLS = "lung-epithelial"

EXP_ROOT = "./microarray"

EXPERIMENTS = dict()
EXP_LABELS = []

chem_mech = {
            }

chemicals = set()
for filename in glob.glob("%s/*_24hour_lung-FC_top_10k_genes.txt" % EXP_ROOT):
    chunks = filename.split('/')[-1].split('_')
    chemical = chunks[0]
    chemicals |= set([chemical])

for chemical in list(chemicals):
    concentrations = []
    for filename in glob.glob("%s/%s_*_24hour_lung-FC_top_10k_genes.txt" % (EXP_ROOT, chemical)):
        concentration = filename.split('/')[-1].split('_')[1]
        concentrations.append(concentration)
    for concentration in concentrations:
        time_point_count = 0
        new_time_series = dict()
        for filename in glob.glob("%s/%s_%s_*hour_lung-FC_top_10k_genes.txt" % (EXP_ROOT, chemical, concentration)):
            time_point = filename.split('/')[-1].split('_')[2]
            new_time_series["%s_%s_%s" % (chemical, concentration, time_point)] = filename
            if os.path.getsize(filename) > 1:
                time_point_count += 1
        if time_point_count > 1:    # We want at the very least two time points with enough genes
            EXPERIMENTS["%s (%s) diXa-004 lung [%s]" % (chemical, concentration, EXP_CELLS)] = new_time_series
            EXP_LABELS.append('?')
        else:
            print("Skipped %s (%s) for not reaching the minimum of two time points with significant gene changes." % (chemical, concentration))

print("done")
