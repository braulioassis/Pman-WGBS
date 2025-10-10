#!/bin/bash
### ---------------------------------------------------
### Fastp preprocessing of paired-end sequencing data
### - Adapter trimming, PolyG trimming, quality filtering
### ---------------------------------------------------

# Define the directory containing raw FASTQ files
DATA_DIR=dir

# Activate conda environment with fastp installed
# (source the miniconda activate script first)
source ~/miniconda3/bin/activate
conda activate fastp

# Change to the data directory
cd $DATA_DIR

# Loop through all files matching the pattern L*_1.fq.gz (R1 reads)
ls L*_1.fq.gz | while read file; do

        # Extract a reference prefix from the file name
        ref=$(echo "${file}" | cut -d "_" -f 1-7) 

        # Extract just the sample name (last part after the last /)
        name=$(echo "${file}" | cut -d "_" -f 1-7 | cut -d "/" -f 7)

        # Print the sample name to track progress
        echo "${name}"

        # Check if output file already exists; skip if so
        if [ -f "${name}.R1.fastp.fq.gz" ]; then
                echo "exists; skipping..."
        else
               ### fastp options:
               ### -g : enable PolyG tail trimming (common in NovaSeq data)
               ### -Q : disable quality filtering
               ### -L : disable length filtering
               ### -A : disable adapter trimming
               
               # Run fastp on paired-end reads
                fastp -g -Q -w 8 \
                	-h "${name}.fastp.html" \  # HTML report
                	-i "${file}" \              # Input R1
                	-I "${ref}_2.fq.gz" \      # Input R2 (paired read)
                	-o "${name}.R1.fastp.fq.gz" \ # Output R1
                	-O "${name}.R2.fastp.fq.gz"   # Output R2
        fi
done

echo "Done!"
