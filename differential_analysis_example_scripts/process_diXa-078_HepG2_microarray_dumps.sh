#!/bin/bash

# Usage:
# $ bash process_diXa-078_HepG2_microarray_dumps.sh /DATA_DUMP_DIR ./microarray


# arguments:
# path to the directory containing the .zip files from TG-Gates
DATADIR=${1}
# path to the output directory to store the lists of relevant genes per concentration/time point
OUTDIR=${2}

# R script to call out differential genes
GENESRSCRIPT='./process_diXa_microarray_4-time-points_3-replicates_1-dose.R'

cd $DATADIR

for X in `ls *.tsv`;
do
  echo "Processing $X";
  Rscript $GENESRSCRIPT $DATADIR $OUTDIR $X;
done

