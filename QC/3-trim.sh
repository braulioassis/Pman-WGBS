#!/bin/bash
### ---------------------------------------------------
### Adapter and quality trimming using Trim Galore
### - Processes paired-end BBduk-trimmed FASTQ files
### - Performs quality trimming and optional filtering
### ---------------------------------------------------

# File containing the list of sample base names (one per line)
samples="samples.txt"

# Loop through each line in the samples file
while read line; do

    # Run Trim Galore on paired-end reads
    # "${line}"*.1.bb*  : matches the R1 file for this sample
    # "${line}"*.2.bb*  : matches the R2 file for this sample
    trim_galore --paired \   # Treat files as paired-end
                 -q 30 \     # Trim low-quality bases with Phred score < 30
                 -j 14 \     # Use 14 threads for parallel processing
                 --length 0 \ # Keep reads of any length (no minimum)
                 "${line}"*.1.bb* \
                 "${line}"*.2.bb* 

# Read lines from the samples file
done < "${samples}"
