# W BW proportion table
base_dir="/"
ext="_CpG.meth"

file_names=(
  W17
W18
W19
W20
W239
W241
W242
W243
W244
W245
W246
W247
W248
W249
W250
W256C
W257
W258
W259
W261
W262
W264
W265
W266A
W267A
W268A
W269
W27A
W27B
W28A
W28B
W29A
W300A
W301A
W302A
W303A
W304
W305
W306
)

files=()
for name in "${file_names[@]}"; do
  files+=("$base_dir/$name$ext")
done

dnmtools merge -t -radmeth "${files[@]}" > "$base_dir/proportion-table-W-BW.txt"

# W ME proportion table
base_dir="/SCRATCH/bassis/wgbs/W/counts/CpG"
ext="_CpG.meth"

file_names=(
W126B
W127B
W128B
W149C
W150C
W152C
W153
W154
W155
W172
W173
W174
W175
W179
W182
W183
W184
W188
W189
W193
W194
W196
W198
W202
W205
W206
W209
W210
W218
W219
W223
W224
W228
W229
W230
W231
W233
W234
W235
)

files=()
for name in "${file_names[@]}"; do
  files+=("$base_dir/$name$ext")
done

dnmtools merge -t -radmeth "${files[@]}" > "$base_dir/proportion-table-W-ME.txt"

# LZ BW proportion table
base_dir="/"
ext="_CpG.meth"

file_names=(
  LZ104
LZ106
LZ108
LZ109
LZ110
LZ111
LZ112
LZ114
LZ116
LZ117
LZ121
LZ122
LZ123
LZ39
LZ40
LZ41
LZ47
LZ48
LZ49
LZ52
LZ63
LZ64
LZ71
LZ87
LZ88
LZ92
LZ93
LZ94
LZ95
LZ96
LZ97
LZ99
LZ101
)

files=()
for name in "${file_names[@]}"; do
  files+=("$base_dir/$name$ext")
done

dnmtools merge -t -radmeth "${files[@]}" > "$base_dir/proportion-table-LZ-BW.txt"

# LZ ME proportion table
base_dir="/SCRATCH/bassis/wgbs/counts/CpG"
ext="_CpG.meth"

file_names=(
LZ10
LZ129
LZ13
LZ130
LZ131
LZ132
LZ135
LZ138
LZ140
LZ146
LZ15
LZ166
LZ168
LZ169
LZ2
LZ21
LZ211
LZ214
LZ23
LZ25
LZ26
LZ3
LZ30
LZ32
LZ6
LZ60
LZ61
LZ66
LZ7
LZ72
LZ74
LZ75
LZ78
LZ81
LZ9
)

files=()
for name in "${file_names[@]}"; do
  files+=("$base_dir/$name$ext")
done

dnmtools merge -t -radmeth "${files[@]}" > "$base_dir/proportion-table-LZ-ME.txt"
