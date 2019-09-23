
# Some unique ID for this set of time series assays; alphanumeric and "-" only for automatic formatting.
EXP_SET_ID = "some-id-for-this-set-1234"

# The cell line or primary cell type used
EXP_CELLS = "HepG2"

# The directory that contains the files with the genes displaying significant changes in expression,
# one file per time point, in the following tab-separated format:
#
# Gene number <TAB> Gene symbol <TAB> log(fold change) <TAB> Adjusted p-value
EXP_ROOT = "./example_assay_results"

EXPERIMENTS = dict()

# Define one entry like the one below for each time series. Use the example format to specify the key:
#           "Chemical-name (concentration) ID-of-assays-set [cell-type]"
# ... as well as the example format to specify each time point's key:
#       "Chemical-name_concentration_time-point" : "/path/to/genes/file/chemical-name_concentration_time-point_top_genes.txt"
# This will help MechSpy pick up this information automatically in the generated output. Please use
# only alphanumeric or "-" characters only for any of these variables.
EXPERIMENTS["SomeChemical (50uM) %s [%s]" % (EXP_SET_ID, EXP_CELLS)] = {
        "SomeChemical_50uM_12h": "%s/SomeChemical_50uM_12h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_50uM_24h": "%s/SomeChemical_50uM_24h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_50uM_48h": "%s/SomeChemical_50uM_48h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_50uM_72h": "%s/SomeChemical_50uM_72h_top_genes.txt" % EXP_ROOT,
        }

EXPERIMENTS["SomeChemical (100uM) %s [%s]" % (EXP_SET_ID, EXP_CELLS)] = {
        "SomeChemical_100uM_12h": "%s/SomeChemical_100uM_12h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_100uM_24h": "%s/SomeChemical_100uM_24h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_100uM_48h": "%s/SomeChemical_100uM_48h_top_genes.txt" % EXP_ROOT,
        "SomeChemical_100uM_72h": "%s/SomeChemical_100uM_72h_top_genes.txt" % EXP_ROOT,
        }

