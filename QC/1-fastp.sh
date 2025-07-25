DATA_DIR=dir

source ~/miniconda3/bin/activate
conda activate fastp

cd $DATA_DIR

ls L*_1.fq.gz  | while read file; do
        ref=$(echo "${file}" | cut -d "_" -f 1-7) 
        name=$(echo "${file}" | cut -d "_" -f 1-7 | cut -d "/" -f 7)
        echo "${name}"
        if [ -f "${name}.R1.fastp.fq.gz" ]; then
                echo "exists; skipping..."
        else
               ### -g (PolyG)
               ### -Q (disable quality filtering)
               ### -L (disable length filtering)
               ### -A (disable adaptor trimming)
                fastp -g -Q -w 8\
                	-h "${name}.fastp.html" \
                	-i "${file}" \
                	-I "${ref}_2.fq.gz" \
                	-o "${name}.R1.fastp.fq.gz" \
                	-O "${name}.R2.fastp.fq.gz"
            
        fi
done
echo "Done!"
