#!/bin/bash
### ---------------------------------------------------
### Adapter trimming using BBduk
### - Processes paired-end FASTQ files previously processed by fastp
### - Trims adapter sequences and optional quality trimming
### ---------------------------------------------------

# Loop through all R1 fastp-processed FASTQ files in subdirectories
ls */*.R1.fastp.fq.gz | while read file; do

        # Extract the base name without file extension
        name=$(echo "${file}" | cut -d "." -f 1)

        # Print sample name
        echo "${name}"

        # Skip processing if output already exists
        if [ -f "${name}.bb.fp.fq.gz" ]; then
                echo "exists; skipping..."
        else 
                # Run BBduk for adapter trimming
                bbduk.sh -Xmx1g \           # Allocate 1 GB memory
                in1="${file}" \              # Input R1
                in2="${name}.R2.fastp.fq.gz" \ # Input R2 (paired read)
                out1="${name}.1.bb.fp.fq.gz" \ # Output R1 after BBduk trimming
                out2="${name}.2.bb.fp.fq.gz" \ # Output R2 after BBduk trimming
                ref=adapters.fa \            # Reference file with adapter sequences
                ktrim=r \                    # Trim adapters from the right end of reads
                k=23 \                       # K-mer size for adapter matching
                mink=11 \                     # Minimum k-mer size for partial matches
                hdist=1 \                     # Hamming distance allowed for k-mer matching
                tpe \                         # Trim paired-end reads to same length
                tbo                           # Trim adapters based on overlap detection
        fi
done

echo "Done!"
