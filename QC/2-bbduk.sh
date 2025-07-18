ls */*.R1.fastp.fq.gz | while read file; do
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
        if [ -f "${name}.bb.fp.fq.gz" ]; then
                echo "exists; skipping..."
        else 
                bbduk.sh -Xmx1g \
                in1="${file}" \
                in2="${name}.R2.fastp.fq.gz" \
                out1="${name}.1.bb.fp.fq.gz" \
                out2="${name}.2.bb.fp.fq.gz" \
                ref=adapters.fa \
                ktrim=r k=23 mink=11 hdist=1 tpe tbo
        fi
done
echo "Done!"
