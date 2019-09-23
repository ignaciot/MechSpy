# MechSpy

## Mechanistic inference for toxicology

This repository provide the code for the publication "Applying Knowledge-Driven Mechanistic Inference to Toxiogenomics" (DOI: xxxxxxxxxxxx) and contains the collection of scripts to reproduce the results obtained. Besides the code, this repository provides:

* Vector embeddings of the knowledge graph used for the inference step, generated using node2vec.
* Knowledge graph utilized for this study (extended from the KG created for [PheKnowLator](https://github.com/callahantiff/PheKnowLator/wiki)), exported as a NetworkX object.
* Differential expression analysis output for all time series tested

The following versions of Python v3.7.2 libraries were used when running this code:

    glob2==0.6
    graphviz==0.10.1
    matplotlib==3.1.0
    natsort==6.0.0
    networkx==2.3
    numpy==1.16.4
    Owlready2==0.18
    pandas==0.24.2
    seaborn==0.9.0
    scikit-learn==0.20.3
    scipy==1.3.0

To clone this repo:

    $ git clone git@github.com:UCDenver-ccp/MechSpy.git
    $ cd MechSpy

This will eventually become an easier to use, single-tool with easier to compose configuration files. In the meantime, these are the steps to generate predictions of the most likely mechanisms of toxicity from your own time series of gene expression data:

1. Identify the most significant changes in gene expression at each time point. You can use your favorite tool for this, as long as the output is formatted as the following tab-separated columns: `Gene number <TAB> Gene symbol <TAB> log(fold change) <TAB> Adjusted p-value`, one file per time point. If your time series was conducted using microarrays, any of the provided R scripts under the `differential_analysis_example_scripts` directory can be modified to match your number of replicates and time points and use Limma to generate this list of differential genes. Be sure to run `install_requirements.R` and install any additional pre-requisites (e.g. system libraries) beforehand.

2. Create a file denoting an experiment set to be used, with a `.py` extension, using as a template the provided `example_experiment_set.py` file. The comments in that file will guide you to complete the necessary details.

3. Call MechSpy's main prediction code!

    $ python3 inference_from_embeddings.py -i example_experiment_set.py

This will output the resulting mechanism enrichment scores and p-values, followed by a sorted list of likely mechanisms of toxicity. It will also generate a pickled file containing all the necessary data to produce an explanation, named `[YOUR_EXPERIMENT_SET_ID]_inference_data.pkl`.

3. Generate the mechanistic narrative and diagram for each of the top-3 predictions. This will output the mechanistic narrative for each time series provided, and create an image with the diagram representation of this mechanistic explanation for each. You must call this other tool with the same input file, and specify an output directory to create the narrative text and the graphical explanation figure:

    $ python3 generate_explanation.py -i my_experiment.py -o ./narratives

Note that this may take a while depending on the size of the knowledge graph used and how dense it is. You can also specify a [case sensitive] search keyword to only generate an explanation for a particular chemical, or concentration:

    $ python3 generate_explanation.py -i my_experiment.py -o ./narratives -k "SomeChemical"
    or
    $ python3 generate_explanation.py -i my_experiment.py -o ./narratives -k "50uM"

For further options on either command, you can call them using:

    $ python3 inference_from_embeddings.py --help

    $ python3 generate_explanation.py --help

You can add new mechanisms of your own to MechSpy by editing `mechanisms.py`. A mock, commented-out "M12" mechanism was added to illustrate where would new mechanisms need to go.

-------------------------------

If you would like to reproduce any of the predictions generated in the MechSpy research article, you can just unpack the processed gene results under the `microarray` directory, and then follow the instructions above with any of the datasets, for example:

    $ python3 inference_from_embeddings.py -i open_tg_gates_canonical_mechs.py
    $ python3 generate_explanation.py -i open_tg_gates_canonical_mechs.py -o ./narratives

These are all the experiment sets you can test:

    open_tg_gates_canonical_mechs.py
    open_tg_gates_all_other_chemicals.py
    dixa_heparg_experiments.py
    dixa_hepg2_experiments.py
    dixa_kidney_experiments.py
    dixa_lung_experiments.py
    dixa078_HepG2_experiments.py
    tobacco_nasal_experiments.py
    tobacco_buccal_experiments.py
    tobacco_bronchial_experiments.py
    all_experiments.py   (runs all of the above, note this will take several hours on a laptop computer)

All known mechanisms for each chemical are listed in `mech_labels.py`.
 
-------------------------------

If you have any questions, please (preferably) open an issue on this GitHub repo via the Issue Tracker tab, or email us at ignacio.tripodi (at) colorado.edu and we'll be happy to help!

If you have used MechSpy or any of the processed data shared in this repository for your research, please cite the following article:

Tripodi, I. J. et al, "Mechanistic inference from knowledge representation: a toxicogenomics case study" [UPDATE URL AND DOI]
    


