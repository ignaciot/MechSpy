#!/bin/bash

# Usage:
# $ bash process_open_tg_gates_microarray_dumps.sh /DIRECTORY_HOLDING_COMPRESSED_DATASETS  ./microarray


# arguments:
# path to the directory containing the .zip files from TG-Gates
ZIPDIR=${1}
# path to the output directory to store the lists of relevant genes per concentration/time point
OUTDIR=${2}

# R script to call out differential genes
GENESRSCRIPT='./process_Open-TG-Gates_microarray_3-time-points_3-doses.R'

cd $ZIPDIR

for X in `ls *.zip`;
do
  echo "Processing $X";
  NEWDIR=${X%.*}
  if [ ! -d $NEWDIR ]; then
    unzip $X;
    Rscript $GENESRSCRIPT $ZIPDIR/$NEWDIR $OUTDIR;
  else
    echo "Already processed";
  fi
done

