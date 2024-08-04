import pandas as pd

#
# --------------------------   INPUT   ------------------------
#

# "../../output/positions/h1n1pdm/ha/positions_nextclade_h1n1pdm_HA1.csv"
positions_merged_csv = snakemake.input.positions_merged_csv
# "../../output/glycosylation/h1n1pdm/ha/glycosylation_nextclade_h1n1pdm_HA1.csv"
glycosylation_merged_csv = snakemake.input.glycosylation_merged_csv

#
# --------------------------   OUTPUT   ------------------------
#

# "../../output/positions/h1n1pdm/ha/positions_nextclade_glycosylation_h1n1pdm_HA1.csv"
merged_csv = snakemake.output.merged_csv

# "../../output/positions/h1n1pdm/ha/positions_nextclade_glycosylation_counts_h1n1pdm_HA1.csv"
counts_csv = snakemake.output.counts_csv

#
# --------------------------   Prepare Dataframes   ------------------------
#

df_positions_merged = pd.read_csv(positions_merged_csv)
df_glycosylation_merged = pd.read_csv(glycosylation_merged_csv)


df = df_positions_merged.merge(
    df_glycosylation_merged, on="seqName", suffixes=('_pos', '_gly'))

# Remove columns that are not of interest.
if "subclade_gly" in df.columns:
    df = df.drop("subclade_gly", axis=1)
if "short-clade_gly" in df.columns:
    df = df.drop("short-clade_gly", axis=1)
if "clade_gly" in df.columns:
    df = df.drop("clade_gly", axis=1)
if "Unnamed: 0_gly" in df.columns:
    df = df.drop("Unnamed: 0_gly", axis=1)
if "Unnamed: 0_pos" in df.columns:
    df = df.drop("Unnamed: 0_pos", axis=1)


# Renaming columns
if "clade_pos" in df.columns:
    df = df.rename(columns={"clade_pos": "clade"})
if "short-clade_pos" in df.columns:
    df = df.rename(columns={"short-clade_pos": "short-clade"})
if "subclade_pos" in df.columns:
    df = df.rename(columns={"subclade_pos": "subclade"})

print(df.head())


# Counts per seqment
df_count_seqName_per_motif = df.groupby(["motif_pos", "motif_gly", "clade"])[
    ["seqName"]].count().sort_values(by="seqName", ascending=False)

print(df_count_seqName_per_motif)

print(f"Total count of seqName: {df_count_seqName_per_motif["seqName"].sum()}")

#
# --------------------------   Save Dataframes   ------------------------
#

# Save the results to CSV.
df.to_csv(merged_csv, index=True)

df_count_seqName_per_motif.to_csv(counts_csv, index=True)
