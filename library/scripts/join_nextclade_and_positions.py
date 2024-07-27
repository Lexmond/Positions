import pandas as pd

#### JOINING THE POSITIONS TABLE WITH NEXTCLADE TABLE ####

#
# --------------------------   INPUT   ------------------------
#

# Import generated Excel file from rule 'get_positions'
positions_csv = snakemake.input.positions_csv

# Import generated csv file from 'nextclade'
nextclade_csv = snakemake.input.nextclade_csv


#
# --------------------------   OUTPUT   ------------------------
#

# Set the path and name of the final exported CSV output file with added clades.
positions_clades_csv = snakemake.output.positions_clades_csv

# Set the path and name of the final exported CSV output file with added clades
# with only variant positions.
positions_clades_variants_csv = snakemake.output.positions_clades_variants_csv

#
# --------------------------   Read inputs   ------------------------
#

# Read positions_csv with pandas and store as Dataframe in 'df_positions'
df_positions = pd.read_csv(positions_csv)

# Read nextclade csv with pandas and store as Dataframe in 'df_nextclade'
df_nextclade = pd.read_csv(nextclade_csv, sep=";")

# Strip all whitespaces from columns names
df_nextclade.columns = df_nextclade.columns.str.strip()

#
# --------------------------   Join tables   ------------------------
#

# Get 'seqName' and 'clade' columns from df_nextclade
csv_name_clade_columns = df_nextclade[['seqName', 'clade']]
print(csv_name_clade_columns)

# Join two df's by sequence name
joined = df_positions.set_index('seqName').join(
    csv_name_clade_columns.set_index('seqName'))
# joined = pd.merge(df_positions, csv_name_clade_columns, on='Sequence')

# Move the 'clade' column so it is next to the sequence name
clade_column = joined.pop('clade')
joined.insert(0, 'clade', clade_column)

# The 'Sequence' column is returned as an index, here a normal dataframe structure is restored.
joined = joined.reset_index()

#
# --------------------------   Add 'motif' column   ------------------------
#

# Concatenate all amino acids per row for each column that has the string "Pos" in the header.
joined['motif'] = joined[[c for c in joined.columns if 'Pos' in c]].sum(axis=1)

# Sort "joined" by  "clade", "motif", "name".
joined = joined.sort_values(
    ['motif', 'clade', "seqName"], ascending=[True, True, True])


#
# --------------SHOW ONLY AA WITH VARIATIONS IN COLUMNS---------------
#


def process_dataframe(df):
    """
    Function to drop position columns that don't show any variation.
    The remaining rows are sorted by similarity afterwards for grouping.
    """
    # Drop columns with only one unique value and column names starting with 'pos'
    columns_to_drop = []
    for col in df.columns:
        if col.startswith('Pos'):
            unique_values = df[col].nunique()
            if unique_values <= 1:
                columns_to_drop.append(col)

    df = df.drop(columns=columns_to_drop)

    # Group rows by sorting similarity of values within 'Pos' columns
    pos_columns = [col for col in df.columns if col.startswith('Pos')]
    df['similarity'] = df[pos_columns].apply(
        lambda row: ''.join(sorted(row)), axis=1)
    df = df.sort_values(by='similarity').drop(
        columns='similarity').reset_index(drop=True)

    return df


# Get a dataframe that excluded identical columns, keeping only columns that show variation.
processed_joined = process_dataframe(joined)

#
# --------------------------   Saving all output to csv files ------------------------
#

# Save the results with clades showing all culumns and motif to CSV.
joined.to_csv(positions_clades_csv, index=True)

# Save the results with clades showing only culumns with variance and motif CSV.
processed_joined.to_csv(positions_clades_variants_csv, index=True)
