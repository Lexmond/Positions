import pandas as pd
from Bio import SeqIO

#
# --------------------------   INPUT   ------------------------
#

# Translated amino acid sequence alignment from Nextclade
aa_fasta = str(snakemake.input.aa_fasta)

# Set 'positions' to look for per lineage and segment using the input Excel file
positions_table = str(snakemake.input.positions_table)

# Use Snakemake wildcards to set the lineage filter for the input Excel file.
filter_lineage = snakemake.wildcards.lineage  # snakemake.input.filter_lineage

# Use Snakemake wildcards to set the segment filter for the input Excel file.
filter_segment = snakemake.wildcards.segment  # snakemake.input.filter_segment

#
# --------------------------   OUTPUT   ------------------------
#

# Table showing what amino acid is found on specified positions for given sequences.
csv_output = snakemake.output.positions_csv

#########################################################################################
#                                                                                       #
# FUNCTION: Get headers and amino acid sequences from Nextaligns output files.          #
#                                                                                       #
#########################################################################################

#


def get_sequence_names_and_amino_acid_sequences(fasta_file):
    """
    Function to get header names and amino acid sequences
    from Nextaligns output files.
    """
    # Read the FASTA file and extract the nucleotide sequences.
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # If there are no sequences in the fasta file, abort the function.
    if len(records) < 1:
        print("No sequences found in the FASTA file.")
        return

    # Get a list of sequence names from the headers.
    sequence_names = []
    for record in records:
        sequence_names.append(record.id)

    # Get a list of corresponding amino acid sequences
    amino_acid_sequences = []
    for record in records:
        amino_acid_sequences.append(record.seq)

    return sequence_names, amino_acid_sequences


#########################################################################################
#                                                                                       #
# FUNCTION: Get a list of positions of interest for lineage and segment of interest.    #
#                                                                                       #
#########################################################################################

def get_positions_from_excel(excel_file, lineage, segment):
    """
    Get a list of positions of interest for sequences from
    lineage and segment of interest.
    """
    # Import pandas
    import pandas as pd

    # Use pandas to read columns from the Excel file
    position_by_lineage_and_segment = pd.read_excel(excel_file)

    # Make a 'selection' by filtering data to match positions for 'ha' segments belonging to the 'vic'lineage.
    selection = position_by_lineage_and_segment.query(
        f'lineage=="{lineage}" & segment=="{segment}"')

    # Get the selected values form the 'positions' columns as a list.
    positions = selection['position'].to_list()

    return positions


#########################################################################################
#                                                                                       #
# FUNCTION: Make an Excel file listing headers with amino acids at given positions.     #
#                                                                                       #
#########################################################################################

def generate_excel_output(output, sequence_names, amino_acid_sequences, positions):
    """
    Make an Excel file listing headers with amino acids at given positions.
    """
    data = []
    for sequence_name, amino_acid_sequence in zip(sequence_names, amino_acid_sequences):
        row = [sequence_name]
        row.extend(amino_acid_sequence[pos - 1] for pos in positions)
        data.append(row)

    column_names = ["seqName"] + [f"Pos {pos}" for pos in positions]

    df = pd.DataFrame(data, columns=column_names)

    df.to_excel(output, index=False)

#########################################################################################
#                                                                                       #
# FUNCTION: Make an CSV file listing headers with amino acids at given positions.     #
#                                                                                       #
#########################################################################################


def generate_csv_output(output, sequence_names, amino_acid_sequences, positions):
    """
    Make an CSV file listing headers with amino acids at given positions.
    """
    data = []
    for sequence_name, amino_acid_sequence in zip(sequence_names, amino_acid_sequences):
        row = [sequence_name]
        row.extend(amino_acid_sequence[pos - 1] for pos in positions)
        data.append(row)

    column_names = ["seqName"] + [f"Pos {pos}" for pos in positions]

    df = pd.DataFrame(data, columns=column_names)

    df.to_csv(output, index=False)

#########################################################################################
#                                                                                       #
#  Run the functions above...                                                           #
#                                                                                       #
#########################################################################################


# Get 'aa_names' and 'aa_sequences' from input fasta file:
aa_names, aa_sequences = get_sequence_names_and_amino_acid_sequences(
    fasta_file=aa_fasta)

# Get the input positions from the input excel_file that are filtered by lineage and segment:
filtered_positions = get_positions_from_excel(
    excel_file=positions_table,
    lineage=filter_lineage,
    segment=filter_segment
)

# Create an Excel file as final output:
if aa_names and aa_sequences:
    generate_csv_output(
        output=csv_output,
        sequence_names=aa_names,
        amino_acid_sequences=aa_sequences,
        positions=filtered_positions
    )

#
# -------------------------    CREATE OUTPUT EXCEL FILE   -----------------------------
#
# Set the name of the final output file (SNAKEFILE):
# nextstrain.output.excel_output
# '../output/aa_at_positions_for_vic_PA.xlsx'
# excel_output = snakemake.output.positions_excel

# # Create an Excel file as final output:
# if aa_names and aa_sequences:
#     generate_excel_output(
#         output=excel_output,
#         sequence_names=aa_names,
#         amino_acid_sequences=aa_sequences,
#         positions=filtered_positions
#     )
