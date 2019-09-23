import pandas as pd
import numpy as np
import re

LUNG_FC_FILENAME = "./public_datasets/lung_log2_fold-change_ratios.txt"
CHIP_ANNOTATION_FILENAME = "./public_datasets/probes_annotation_HG-U133_Plus_2.txt"
SAMPLE_INFO = "./public_datasets/s_Lung.txt"
OUTPUT_DIR = "./microarray"


fc_df = pd.read_csv(LUNG_FC_FILENAME, sep='\t', na_filter=False)

sample_df = pd.read_csv(SAMPLE_INFO, sep='\t', na_filter=False)

annotations = dict()
with open(CHIP_ANNOTATION_FILENAME, 'r') as fd:
    for line in fd.readlines():
        chunks = line.split('\t')
        annotations[chunks[0]] = chunks[4]

for sample in fc_df.columns[1:].values:
    sample_info = sample_df[sample_df['Source Name'] == sample]
    if len(sample_info) == 0:
        print("Sample %s not found, skipping..." % sample)
        continue

    symbols = []
    fc = []
    for row in fc_df[['Row Names', sample]].itertuples():
        fc_value = 0
        if row[2] != 'NaN':
            fc_value = abs(float(row[2]))
        fc.append(fc_value)
        symbol = annotations[row[1]]
        if symbol == '':
            symbol = 'NA'
        symbols.append(symbol.split(" : ")[0])
    fc = np.array(fc)

    chem_name = sample_info['Factor Value[Compound]'].values[0]
    chem_name = re.sub("[^0-9a-zA-Z]", "-", chem_name)
    chem_dose = sample_info['Factor Value[Dose]'].values[0]
    chem_dose_unit = sample_info['Characteristics[DoseUnit]'].values[0]
    time_point = sample_info['Factor Value[Dose Duration]'].values[0]
    time_point_unit = sample_info['Characteristics[Dose DurationUnit]'].values[0]
    new_filename = "%s_%s-%s_%s%s_lung-FC_top_10k_genes.txt" % (chem_name,
                                                                chem_dose,
                                                                chem_dose_unit,
                                                                time_point,
                                                                time_point_unit)
    with open("%s/%s" % (OUTPUT_DIR, new_filename), 'w') as new_fd:
        gene_count = 0
        new_fd.write("\t\tdiff_exp.Symbol\tdiff_exp.logFC\tdiff_exp.adj.P.Val\n")
        # iterate over indices of fold change, sorted in descending order
        for idx in fc.argsort()[::-1]:
            gene_count += 1
            if gene_count > 10000:
                break
            # use a fake "p-value" for these to use the existing codebase as is (since we
            # only have one replicate of each condition we can't really use limma or any 
            # statistical analysis other than just fold change)
            new_fd.write("%s\t%s\t%s\t%s\n" % (gene_count, symbols[idx], fc[idx], 1/(1000*fc[idx])))

