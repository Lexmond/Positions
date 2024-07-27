import pandas as pd

#### JOINING THE POSITIONS TABLE WITH NEXTCLADE OUTPUT ####

# Import generated Excel file from rule 'get_positions'
positions_csv = snakemake.input.positions_csv

# Import generated csv file from 'nextclade'
nextclade_csv = snakemake.input.nextclade_csv

# Set the path and name of the final exported Excel output file.
positions_clades_xlsx = snakemake.output.positions_clades_xlsx

# Set the path and name of the final exported CSV output file with added clades.
positions_clades_csv = snakemake.output.positions_clades_csv

# Set the path and name of the final exported CSV output file with added clades
# with only variant positions.
positions_clades_variants_csv = snakemake.output.positions_clades_variants_csv

# positions_clades_variants_counts_csv = snakemake.output.positions_clades_variants_counts_csv


# Read positions_csv with pandas and store Dataframe in 'df_positions'
df_positions = pd.read_csv(positions_csv)

# Read csv sheet with pandas and store Dataframe in 'df_nextclade'
df_nextclade = pd.read_csv(nextclade_csv, sep=";")

# Strip whitespaces from columns names
df_nextclade.columns = df_nextclade.columns.str.strip()

# Get 'seqName' and 'clade' columns from df_nextclade
csv_name_clade_columns = df_nextclade[['seqName', 'clade']]
print(csv_name_clade_columns)

# csv_name_clade_columns.rename(columns={'seqName': 'Sequence'}, inplace=True)
# print(csv_name_clade_columns)

# Join two df's by sequence name
joined = df_positions.set_index('seqName').join(
    csv_name_clade_columns.set_index('seqName'))
# joined = pd.merge(df_positions, csv_name_clade_columns, on='Sequence')


# Move the 'clade' column so it is next to the sequence name
clade_column = joined.pop('clade')
joined.insert(0, 'clade', clade_column)

# The 'Sequence' column is returned as an index, here a normal dataframe structure is restored.
joined = joined.reset_index()

# Save the results to a CSV file
joined.to_csv(positions_clades_csv, index=True)

#### PROCESS THE DATAFRAME TO ONLY HIGHLIGHT VARIATIONS IN COLUMNS ####

"""
Function to drop position columns that don't show any variation.
The remaining rows are sorted by similarity afterwards for grouping.
"""


def process_dataframe(df):
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

# Save the results to a CSV file
processed_joined.to_csv(positions_clades_variants_csv, index=True)

#### COUNTING VARIATIONS ####

# # Select the 'position' columns (third to last columns) for counting
# pos_columns = processed_joined.columns[2:]

# # Count the variation per row in the above defined columns.
# counts = processed_joined[pos_columns].value_counts()
# print(counts)

# # Save the results to a CSV file
# counts.to_csv(positions_clades_variants_counts_csv, index=True)

# Write dataframes 'joined' and 'processed_joined' per sheet to one excel file
with pd.ExcelWriter(positions_clades_xlsx) as writer:
    joined.to_excel(
        writer, sheet_name="All positions", index=False)
    processed_joined.to_excel(
        writer, sheet_name="Variant positions", index=False)
    # counts.to_excel(
    #     writer, sheet_name="Counts of variants")

# Finally, export the new Excel file with the added clade column.
# processed_joined.to_excel(excel_with_clades, sheet_name="positions_processed")
