from Bio import SeqIO
import pandas as pd
import re

# Input from Snakefile
# Note that the input from aggregate_translations needs to be converted to a string for SeqIO to work.
ha_sequences = str(snakemake.input.aa_fasta)

# Output file for Snakefile
glycosylation_csv = snakemake.output.glycosylation_csv

# Define the regular expression pattern for glycosylation sites
glyco_pattern = re.compile(r'N[^P][ST][^P]')

# Read the sequences from the FASTA file
sequences = list(SeqIO.parse(ha_sequences, "fasta"))

# Create a dictionary to store the results
results = {}

# Loop through each sequence and find the glycosylation sites
for seq_record in sequences:
    seq = str(seq_record.seq)
    name = seq_record.id
    matches = glyco_pattern.finditer(seq)

    # Store the found patterns with their positions
    results[name] = {match.start() + 1: match.group() for match in matches}

# Convert the results dictionary to a DataFrame
df = pd.DataFrame(results).T

# Sort the columns (positions) numerically
df = df.reindex(sorted(df.columns, key=int), axis=1)

# Set the name of the index column
df.index.name = "seqName"

# Rename the columns by adding "Pos " to the existing values
df.columns = [f"Pos {col}" if col != 'seqName' else col for col in df.columns]

# Display the DataFrame
print(df)

# Save the results to a CSV file
df.to_csv(glycosylation_csv, index=True)
