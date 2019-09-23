import glob
import os

EXP_SET_ID = "rest-Open-TG-Gates"
EXP_CELLS = "hepatocytes"

EXP_ROOT = "./microarray"

EXPERIMENTS = dict()
EXP_LABELS = []

chem_mech = {
                'adapin': '?',
                'allopurinol': '?',
                'allyl-alcohol': 'M4',
                'amiodarone': 'M5',
                'aspirin': 'M4',
                'azathioprine': '?',
                'benzbromarone': '?',
                'bromobenzene': '?',
                'carbamazepine': '?',
                'chlorpromazine': '?',
                'cimetidine': '?',
                'clofibrate': '?',
                'cyclophosphamide': '?',
                'diazepam': '?',
                'ethionine': '?',
                'fluphenazine': '?',
                'flutamide': '?',
                'gemfibrozil': '?',
                'glibenclamide': '?',
                'griseofulvin': '?',
                'haloperidol': '?',
                'hexachlorobenzene': '?',
                'hydroxyzine': 'M6',
                'ibuprofen': '?',
                'imipramine': '?',
                'indomethacin': '?',
                'interleukin-1-beta,-human': '?',
                'interleukin-6,-human': '?',
                'ketoconazole': '?',
                'labetalol': '?',
                'lomustine': '?',
                'methapyrilene': '?',
                'methyltestosterone': '?',
                'naphthyl-isothiocyanate': '?',
                'nitrofurantoin': '?',
                'omeprazole': '?',
                'perhexiline': '?',
                'phenylbutazone': '?',
                'phenytoin': '?',
                'propylthiouracil': '?',
                'rifampicin': '?',
                'sulfasalazine': '?',
                'tetracycline': '?',
                'thioacetamide': '?',
                'thioridazine': '?',
                'transforming-growth-factor-beta-1': '?',
                'WY-14643': '?',
            }

chemicals = set()
for filename in glob.glob("%s/*-uM_2h_top_10k_genes.txt" % EXP_ROOT):
    chunks = filename.split('/')[-1].split('_')
    chemical = chunks[0]
    chemicals |= set([chemical])

for chemical in list(chemicals):
    concentrations = []
    for filename in glob.glob("%s/%s_*-uM_2h_top_10k_genes.txt" % (EXP_ROOT, chemical)):
        concentration = filename.split('/')[-1].split('_')[1]
        concentrations.append(concentration)
    for concentration in concentrations:
        time_point_count = 0
        new_time_series = dict()
        for filename in glob.glob("%s/%s_%s_*h_top_10k_genes.txt" % (EXP_ROOT, chemical, concentration)):
            time_point = filename.split('/')[-1].split('_')[2]
            new_time_series["%s_%s_%s" % (chemical, concentration, time_point)] = filename
            if os.path.getsize(filename) > 1:
                time_point_count += 1
        if time_point_count > 1:    # We want at the very least two time points with enough genes
            EXPERIMENTS["%s (%s) Open TG-Gates [%s]" % (chemical, concentration, EXP_CELLS)] = new_time_series
            EXP_LABELS.append(chem_mech[chemical])
        else:
            print("Skipped %s (%s) for not reaching the minimum of two time points with significant gene changes." % (chemical, concentration))

print("done")
