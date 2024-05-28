#!/bin/bash

## EDIT THESE VARIABLES ONLY.
#
# Control: the control ChIPseq file.
# Experimental: relative or absolute paths to each experimental SAM file.
# Genome: Genome to use with HOMER.

DATA_DIRECTORY="results/aligned/raw"
GENOME="hg38"

## DO NOT EDIT.
#
# Finds peaks and annotates peaks with HOMER. Makes a BED file for later used
# with BEDtools.

CONTROL_FILE=
EXPERIMENTAL_FILES=
tmp=($(ls $DATA_DIRECTORY))
for nm in ${tmp[@]}; do
    if [[ ${nm^^} == CONTROL* ]]; then
        CONTROL_FILE=${DATA_DIRECTORY:-.}/$nm
    else
        EXPERIMENTAL_FILES=(${EXPERIMENTAL_FILES[@]} ${DATA_DIRECTORY:-.}/$nm)
    fi
done

echo "Searching directory: $DATA_DIRECTORY"
echo "Found control file: $CONTROL_FILE"
echo "Found ${#EXPERIMENTAL_FILES[@]} experimental files: ${EXPERIMENTAL_FILES[@]}"

manual_dir=./_manual/tagdirs
mkdir -p $manual_dir

ctl_base=${CONTROL_FILE##*/}
ctl_ext=${ctl_base##*\.}
ctl_base=${ctl_base%.*}

makeTagDirectory $manual_dir/$ctl_base $CONTROL_FILE -genome $GENOME

for expr_file in ${EXPERIMENTAL_FILES[@]}; do
    base=${expr_file##*/}
    ext=${base##*\.}
    base=${base%.*}

    makeTagDirectory $manual_dir/$base $expr_file -genome $GENOME
    findPeaks $manual_dir/$base -style factor -o auto -i $manual_dir/$ctl_base
    annotatePeaks.pl $manual_dir/$base/peaks.txt $GENOME -annStats \
        $manual_dir/$base/peaks-annStats.txt > $manual_dir/$base/peaks-ann.txt
    pos2bed.pl $manual_dir/$base/peaks.txt > $manual_dir/$base/peaks.ped
done