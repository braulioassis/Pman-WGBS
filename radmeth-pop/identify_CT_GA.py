import pysam

def identify_c_to_t_and_g_to_a_substitutions(sam_file_path, output_bed_path):
    """
    Identify C-to-T and G-to-A substitutions from a SAM file and export their positions to a BED file.

    Args:
    - sam_file_path (str): Path to the input SAM file containing aligned sequencing reads.
    - output_bed_path (str): Path to the output BED file where substitution positions will be saved.
    """

    try:
        # Open the SAM file for reading, and the BED file for writing
        with pysam.AlignmentFile(sam_file_path, "r") as samfile, open(output_bed_path, "w") as bedfile:
            substitution_count = 0  # Keeps track of the number of detected substitutions
            processed_reads = 0     # Keeps track of total processed reads for reporting

            # Iterate through all reads in the SAM file
            for read in samfile:
                processed_reads += 1  # Increment total reads counter

                # Skip reads that are not aligned to the reference genome
                if read.is_unmapped:
                    continue

                # Skip reads missing a query sequence (e.g., incomplete or malformed entries)
                if not read.query_sequence:
                    print(f"Skipping read {read.query_name}: missing sequence data.")
                    continue

                # Get the starting reference position and the reference chromosome name
                ref_start = read.reference_start
                ref_name = samfile.get_reference_name(read.reference_id)

                # Extract the MD tag — this encodes mismatches between the read and reference
                # Example MD tag: "10A5^C3T2" (means 10 matches, mismatch A, 5 matches, deletion C, etc.)
                md_tag = dict(read.tags).get("MD")

                # Skip reads that don’t contain an MD tag (required to identify mismatches)
                if not md_tag:
                    print(f"Skipping read {read.query_name}: missing MD tag.")
                    continue

                # --- Parse the MD tag into a structured list of match counts and mismatch bases ---
                ref_pos = ref_start  # Tracks the current position on the reference
                seq_pos = 0          # Tracks the current position in the read sequence
                md_parts = []        # Holds parsed MD tag components (integers and mismatched bases)
                num = ''             # Temporary holder for digits representing match counts

                # Parse MD tag into components (numbers = matches, letters = mismatches)
                for char in md_tag:
                    if char.isdigit():
                        num += char
                    else:
                        if num:
                            md_parts.append(int(num))
                            num = ''
                        md_parts.append(char)
                if num:
                    md_parts.append(int(num))  # Append trailing number, if any

                # --- Identify C→T and G→A substitutions based on mismatches ---
                for part in md_parts:
                    if isinstance(part, int):
                        # A number represents matching bases; advance both reference and read positions
                        ref_pos += part
                        seq_pos += part
                    elif isinstance(part, str):
                        # A letter represents a mismatch at this position
                        if seq_pos < len(read.query_sequence):
                            read_base = read.query_sequence[seq_pos]  # Base from the read
                            ref_base = part                            # Base from the reference (via MD tag)

                            # Check for bisulfite-relevant substitutions:
                            #  - C→T on the read strand
                            #  - G→A on the opposite strand
                            if (ref_base == "C" and read_base == "T") or (ref_base == "G" and read_base == "A"):
                                # Write the substitution position to the BED file:
                                # Format: chrom, start, end (1-based half-open)
                                bedfile.write(f"{ref_name}\t{ref_pos}\t{ref_pos + 1}\n")
                                substitution_count += 1  # Increment counter

                        # Move both reference and read positions forward by one after each mismatch
                        ref_pos += 1
                        seq_pos += 1

            # After all reads are processed, report summary
            print(f"Processed {processed_reads} reads. Identified {substitution_count} substitutions.")

    except Exception as e:
        # Catch and display errors
        print(f"Error: {e}")

# --- Define input/output paths ---
input_sam_file = "ME2.1.3.sam"  # Input SAM file (aligned sequencing reads)
output_bed_file = "CT_and_GA_positions2.1.3.bed"  # Output BED file for substitution positions

# --- Run the function ---
identify_c_to_t_and_g_to_a_substitutions(input_sam_file, output_bed_file)

# Confirmation message
print(f"C-to-T and G-to-A substitution positions exported to {output_bed_file}")
