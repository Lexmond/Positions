import pandas as pd

#
# --------------------------   INPUT   ------------------------
#

# Set 'positions' to look for per lineage and segment using the input Excel file

nextclade_csv = snakemake.input.nextclade_csv
glycosylation_csv = snakemake.input.glycosylation_csv

#
# --------------------------   INPUT   ------------------------
#

merged_csv = snakemake.output.merged_csv

#
# --------------------------   Prepare Dataframes   ------------------------
#

# H1N1pdm PA data
df_nextclade = pd.read_csv(nextclade_csv, sep=";")
df_glycosylation = pd.read_csv(glycosylation_csv, sep=",")

# Reset the indexes for the dataframes
df_nextclade = df_nextclade.reset_index(drop=True)
df_glycosylation = df_glycosylation.reset_index(drop=True)

# Drop the surplus column "index" from the nextclade dataframe
df_nextclade = df_nextclade.drop(["index"], axis=1)

#
# --------------------------  Merge Dataframes   ------------------------
#

# Merge all three tables into one Dataframe 'df'.
df = df_glycosylation.merge(df_nextclade, how="left", on="seqName")

#
# ------  Remove unwanted columns from nextclade after merge   -----------
#

# Check if "subclade" column exists, this is true for nextclade output for ha
if "subclade" in df.columns:
    # Delete all columns to the right of "subclade", keeping columns "short-clade" and "subclade"
    cols_to_keep = df.columns[:df.columns.get_loc("subclade") + 1]
    df = df[cols_to_keep]
else:
    # Delete all columns to the right of "clade", in case that columns "short-clade" and "subclade" are non-exsisting.
    cols_to_keep = df.columns[:df.columns.get_loc("clade") + 1]
    df = df[cols_to_keep]

#
# ------  Sort and renaming Pos columns and adding 'motif' column  -----------
#

# Select columns with names that start with 'Pos '
pos_columns = [col for col in df.columns if col.startswith('Pos ')]

# Identify pos_columns that show variation in their values across the rows
variable_pos_columns = [col for col in pos_columns if df[col].nunique() > 1]

# Create the 'motif' column by concatenating values of the variable 'Pos ' columns
df['motif'] = df[variable_pos_columns].astype(str).agg(''.join, axis=1)

# Remove 'Pos ' part of the column names
df.rename(columns={col: col.replace('Pos ', '')
          for col in variable_pos_columns}, inplace=True)

# Sort the 'Pos ' columns by the numeric part
sorted_columns = sorted([col.replace('Pos ', '')
                        for col in variable_pos_columns], key=int)

# Determine the final order of columns, considering the presence of optional columns
final_columns = ['seqName', 'clade'] + sorted_columns + ['motif']
if "short-clade" in df.columns:
    final_columns.append('short-clade')
if "subclade" in df.columns:
    final_columns.append('subclade')

if "subclade" in df.columns:
    final_columns = ['seqName', 'clade', 'subclade'] + \
        sorted_columns + ['motif']

# Changin the order of columns
if "short-clade" in df.columns:
    final_columns = ['seqName', 'clade', 'short-clade',
                     'subclade'] + sorted_columns + ['motif']

# Store the reordered DataFrame columns to 'df'.
df = df[final_columns]

print(df)


# HA Save the results to CSV.
df.to_csv(merged_csv, index=True)
