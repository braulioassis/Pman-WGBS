# Radmeth for normoxia vs. hypoxia within each population and tissue type
dnmtools radmeth -v -factor Trt design.w.bw.txt proportion-table.w.bw.txt > radmeth.w.bw.trt.bed
dnmtools radmeth -v -factor Trt design.w.me.txt proportion-table.w.me.txt > radmeth.w.me.trt.bed
dnmtools radmeth -v -factor Trt design.lz.bw.txt proportion-table.lz.bw.txt > radmeth.lz.bw.trt.bed
dnmtools radmeth -v -factor Trt design.lz.me.txt proportion-table.lz.me.txt > radmeth.lz.me.trt.bed
dnmtools radmeth -v -factor Trt design.jz.bw.txt proportion-table.jz.bw.txt > radmeth.jz.bw.trt.bed
dnmtools radmeth -v -factor Trt design.jz.me.txt proportion-table.jz.me.txt > radmeth.jz.me.trt.bed

# Adjust p-value based on neighboring sites
dnmtools radadjust -bins 1:200:1 radmeth.w.bw.trt.bed > radmeth.w.bw.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.w.me.trt.bed > radmeth.w.me.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.lz.bw.trt.bed > radmeth.lz.bw.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.lz.me.trt.bed > radmeth.lz.me.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.jz.bw.trt.bed > radmeth.jz.bw.trt.adj.bed
dnmtools radadjust -bins 1:200:1 radmeth.jz.me.trt.bed > radmeth.jz.me.trt.adj.bed

# Merge into DMR's
dnmtools radmerge -p 0.01 radmeth.w.bw.trt.adj.bed > radmeth.w.bw.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.w.me.trt.adj.bed > radmeth.w.me.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.lz.bw.trt.adj.bed > radmeth.lz.bw.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.lz.me.trt.adj.bed > radmeth.lz.me.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.jz.bw.trt.adj.bed > radmeth.jz.bw.trt.adj.dmr.bed
dnmtools radmerge -p 0.01 radmeth.jz.me.trt.adj.bed > radmeth.jz.me.trt.adj.dmr.bed
