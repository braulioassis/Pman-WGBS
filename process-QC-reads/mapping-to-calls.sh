# DNMtools pipeline
# Index reference
dnmtools abismalidx GCF_003704035.1_HU_Pman_2.1.3_genomic.fna ref-index2.1.3

# Map
samples="samples.txt"
while read line; do
        dnmtools abismal -i ref-index2.1.3 -v -t 28 -B -o "${line}".bam "${line}"*.1.*val* "${line}"*.2.*val*
done < "${samples}"

# Format
while read line; do
	dnmtools format -f abismal -v -t 28 -B "${line}".bam "${line}"Format.bam
	done < "${samples}"

# Sort
while read line; do
	samtools sort -o "${line}"Format_Sort.bam "${line}"Format.bam
	done < "${samples}"

# Remove duplicates
while read line; do
	dnmtools uniq -B -t 28 -S "${line}"duplicate-removal-stats.txt "${line}"Format_Sort.bam "${line}"Format_Sort_Uniq.bam
	done < "${samples}"

# Get bisulfite conversion rates
while read line; do
	dnmtools bsrate -t 28 -v -c GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -o "${line}".bsrate "${line}"Format_Sort_Uniq.bam
	done < "${samples}"

# Compute methylation
while read line; do
	dnmtools counts -t 28 -c GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -o "${line}".meth "${line}"Format_Sort_Uniq.bam
	done < "${samples}"

# Collapse symmetric CpGs
while read line; do
	dnmtools sym -o "${line}"CpG.meth "${line}".meth
	done < "${samples}"

# Get methylation statistics
while read line; do
	dnmtools levels -relaxed -o "${line}".levels "${line}"CpG.meth
	done < "${samples}"
