import pysam

def identify_c_to_t_and_g_to_a_substitutions(sam_file_path, output_bed_path):
    """
    Identify C-to-T and G-to-A substitutions from a SAM file and export positions to a BED file.

    Args:
    - sam_file_path: Path to the input SAM file.
    - output_bed_path: Path to the output BED file.
    """
    try:
        with pysam.AlignmentFile(sam_file_path, "r") as samfile, open(output_bed_path, "w") as bedfile:
            substitution_count = 0  # Counter for identified substitutions
            processed_reads = 0  # Counter for total processed reads

            for read in samfile:
                processed_reads += 1

                if read.is_unmapped:
                    continue

                if not read.query_sequence:
                    print(f"Skipping read {read.query_name}: missing sequence data.")
                    continue

                ref_start = read.reference_start
                ref_name = samfile.get_reference_name(read.reference_id)
                md_tag = dict(read.tags).get("MD")

                if not md_tag:
                    print(f"Skipping read {read.query_name}: missing MD tag.")
                    continue

                # Parse MD tag into list of numbers and mismatched bases
                ref_pos = ref_start
                seq_pos = 0
                md_parts = []
                num = ''
                for char in md_tag:
                    if char.isdigit():
                        num += char
                    else:
                        if num:
                            md_parts.append(int(num))
                            num = ''
                        md_parts.append(char)
                if num:
                    md_parts.append(int(num))

                # Identify C→T and G→A substitutions
                for part in md_parts:
                    if isinstance(part, int):
                        ref_pos += part
                        seq_pos += part
                    elif isinstance(part, str):
                        if seq_pos < len(read.query_sequence):
                            read_base = read.query_sequence[seq_pos]
                            ref_base = part

                             # Check for C→T or G→A substitution
                            if (ref_base == "C" and read_base == "T") or (ref_base == "G" and read_base == "A"):
                                bedfile.write(f"{ref_name}\t{ref_pos}\t{ref_pos + 1}\n")
                                substitution_count += 1

                        ref_pos += 1
                        seq_pos += 1

            print(f"Processed {processed_reads} reads. Identified {substitution_count} substitutions.")

    except Exception as e:
        print(f"Error: {e}")

# Define input SAM file and output BED file paths
input_sam_file = "ME2.1.3.sam"
output_bed_file = "CT_and_GA_positions2.1.3.bed"

# Generate the BED file
identify_c_to_t_and_g_to_a_substitutions(input_sam_file, output_bed_file)
print(f"C-to-T substitution positions exported to {output_bed_file}")
