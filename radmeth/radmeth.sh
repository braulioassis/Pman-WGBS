### For population-level differences
# Radmeth proportion table
dnmtools merge -t -radmeth /W/*CpG*.meth > proportion-table.txt
dnmtools merge -t -radmeth /LZ/*CpG*.meth > proportion-table.txt
dnmtools merge -t -radmeth /JZ/*CpG*.meth > proportion-table.txt

# Radmeth for population, with hypoxia treatment as a factor
dnmtools radmeth -v -factor Pop design.w.txt /W/proportion-table.txt > radmeth.w.pop.trt.bed
dnmtools radmeth -v -factor Pop design.lz.txt /LZ/proportion-table.txt > radmeth.lz.pop.trt.bed
dnmtools radmeth -v -factor Pop design.jz.txt /JZ/proportion-table.txt > radmeth.jz.pop.trt.bed

### Filter out C-to-T substitution sites in ME population
# Index reference using hisat2
hisat2-build -p 28 GCF_003704035.1_HU_Pman_2.1.3_genomic.fna ht2

# Trim raw ME WGS
trim_galore --paired -q 0 -j 14 --length 20 *1.fq.gz *2.fq.gz

# Map with hisat2
hisat2 -x ht2 \
-1 *val_1.fq.gz -2 *val_2.fq.gz \
-S ME-2.1.3

# Obtain C-to-T substitutions in ME
# Run identify_c_to_t.py

### Build DMR's 
# Adjust p-values based on neighboring loci
dnmtools radadjust -bins 1:200:1 radmeth.w.pop.trt.bed > radmeth.w.pop.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.lz.pop.trt.bed > radmeth.lz.pop.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.jz.pop.trt.bed > radmeth.jz.pop.trt.adj.bed

# Merge methylation sites into DMR's
dnmtools radmerge -p 0.01 radmeth.w.pop.trt.adj.bed > radmeth.w.pop.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.lz.pop.trt.adj.bed > radmeth.lz.pop.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.jz.pop.trt.adj.bed > radmeth.jz.pop.trt.adj.dmr.bed
