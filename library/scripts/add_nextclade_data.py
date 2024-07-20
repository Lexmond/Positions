import pandas as pd

#### JOINING THE POSITIONS TABLE WITH NEXTCLADE OUTPUT ####

# Import generated Excel file from rule 'get_positions'
excel_output = snakemake.input.excel_output

# Import generated csv file from 'nextclade'
nextclade_csv = snakemake.input.nextclade_csv

# Set the path and name of the final exported Excel output file.
excel_with_clades = snakemake.output.excel_with_clades


# Read Excel sheet with pandas and store Dataframe in 'xlsx'
xlsx = pd.read_excel(excel_output, sheet_name="Sheet1")

# Read csv sheet with pandas and store Dataframe in 'csv'
csv = pd.read_csv(nextclade_csv, sep=";")

# Strip whitespaces from columns names
csv.columns = csv.columns.str.strip()

# Get 'seqName' and 'clade' columns from the csv data
csv_name_clade_columns = csv[['seqName', 'clade']]
print(csv_name_clade_columns)

# csv_name_clade_columns.rename(columns={'seqName': 'Sequence'}, inplace=True)
# print(csv_name_clade_columns)

# Join two df's by sequence name
joined = xlsx.set_index('Sequence').join(
    csv_name_clade_columns.set_index('seqName'))
# joined = pd.merge(xlsx, csv_name_clade_columns, on='Sequence')


# Move the 'clade' column so it is next to the sequence name
clade_column = joined.pop('clade')
joined.insert(0, 'clade', clade_column)

# The 'Sequence' column is returned as an index, here a normal dataframe structure is restored.
joined = joined.reset_index()


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

#### COUNTING VARIATIONS ####

# Select the 'position' columns (third to last columns) for counting
pos_columns = processed_joined.columns[2:]

# Count the variation per row in the above defined columns.
counts = processed_joined[pos_columns].value_counts()
print(counts)

# Write dataframes 'joined' and 'processed_joined' per sheet to one excel file
with pd.ExcelWriter(excel_with_clades) as writer:
    joined.to_excel(
        writer, sheet_name="raw_table", index=False)
    processed_joined.to_excel(
        writer, sheet_name="processed_table", index=False)
    counts.to_excel(
        writer, sheet_name="counts")

# Finally, export the new Excel file with the added clade column.
# processed_joined.to_excel(excel_with_clades, sheet_name="positions_processed")
