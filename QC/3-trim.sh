samples="samples.txt"
while read line; do
	trim_galore --paired -q 30 -j 14 --length 0 "${line}"*.1.bb* "${line}"*.2.bb*;
	done < "${samples}"
