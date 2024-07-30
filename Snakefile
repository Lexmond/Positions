from pathlib import Path

# Pth of the 'input' directory
INPUT = Path("input")

# Gene lookup dictionary by seqments
GENES = {
    "pb2": "PB2",
    "pb1": "PB1",
    "pa": "PA",
    "ha": "HA1",
    "np": "NP",
    "na": "NA",
    "mp": "M1",
    "ns": "NS1",
}

# Function to scan the input folder and generate lists of lineages, segments, and genes
def input_func_nextalign():
    l1, l2, l3 = [], [], []
    for path in list(INPUT.glob("*.fasta")):
        _, lineage, segment = path.stem.split('_')  # sequences_h1n1pdm_ha.fasta
        gene = GENES[segment]  # Get gene variable from dict by segment
        l1.append(lineage)
        l2.append(segment)
        l3.append(gene)
    return {
        "lineage": l1,
        "segment": l2,
        "gene": l3,
    }

# Get input data dictionary as a dictionary of lists
input_data_nextalign = input_func_nextalign()


# Create lists of tuples that can be used to extract zipped values from the different lists
list_of_2tuples = list(zip(input_data_nextalign['lineage'], input_data_nextalign['segment']))
list_of_3tuples = list(zip(input_data_nextalign['lineage'], input_data_nextalign['segment'], input_data_nextalign['gene']))


# Create a list of output files used rule "nextclade"
nextclade_output_files = [f"output/nextclade/flu/{lineage}/{segment}/nextclade.csv" \
for lineage, segment in list_of_2tuples]
print(f"All output files for checkpoint 'nextclade': \n {nextclade_output_files} \n")

# Create a list of output files used by rule "get_positions"
get_positions_output_files = [f"output/positions/{lineage}/{segment}/positions_{lineage}_{gene}.csv" \
for lineage, segment, gene in list_of_3tuples]
print(f"All output files for rule 'get_positions': \n {get_positions_output_files} \n")

# Create a list of output files used by rule "get_glycosylation"
get_glycosylation_output_files = [f"output/glycosylation/{lineage}/{segment}/glycosylation_{lineage}_{gene}.csv" \
for lineage, segment, gene in list_of_3tuples]
print(f"All output files for rule 'get_glycosylation': \n {get_glycosylation_output_files} \n")

# Merge the output files from rules "nextclade" and "get_positions"
merge_positions_nextclade_output_files = [f"output/positions/{lineage}/{segment}/positions_nextclade_{lineage}_{gene}.csv" \
for lineage, segment, gene in list_of_3tuples]
print(f"All output files for rule 'merge_positions_nextclade': \n {merge_positions_nextclade_output_files} \n")

# Merge the output files from rules "nextclade" and "get_glycosylation"
merge_glycosylation_nextclade_output_files = [f"output/glycosylation/{lineage}/{segment}/glycosylation_nextclade_{lineage}_{gene}.csv" \
for lineage, segment, gene in list_of_3tuples]
print(f"All output files for rule 'merge_glycosylation_nextclade': \n {merge_glycosylation_nextclade_output_files} \n")

# Merge the output files from rules "merge_positions_nextclade" and "merge_glycosylation_nextclade"
merge_positions_nextclade_glycosylation_output_files = [f"output/positions/{lineage}/{segment}/positions_nextclade_glycosylation_{lineage}_{gene}.csv" \
for lineage, segment, gene in list_of_3tuples]
print(f"All output files for rule 'merge_positions_nextclade_glycosylation': \n {merge_positions_nextclade_glycosylation_output_files} \n")


# Snakemake rules
rule all:
    input:
        nextclade_output_files,
        get_glycosylation_output_files,
        get_positions_output_files,
        merge_positions_nextclade_output_files,
        merge_glycosylation_nextclade_output_files,
        merge_positions_nextclade_glycosylation_output_files


# Run Nextclade datasets on the input sequences to get outputs per subtype and segment 
checkpoint nextclade:
    input:
        fasta = "input/sequences_{lineage}_{segment}.fasta",
        dataset = "library/nextclade/data/flu/{lineage}/{segment}"
    output:
        nextclade_csv = "output/nextclade/flu/{lineage}/{segment}/nextclade.csv",
        translations = directory("output/nextclade/flu/{lineage}/{segment}/")  # maak directory en cat om naar checkpoint...
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-all=output/nextclade/flu/{wildcards.lineage}/{wildcards.segment} \
            {input.fasta}
        """

def aggregate_translations(wildcards):
    """The alignment rule produces multiple outputs that we cannot easily name prior
    to running the rule. The names of the outputs depend on the segment being
    aligned and Snakemake's `expand` function does not provide a way to lookup
    the gene names per segment. Instead, we use Snakemake's checkpoint
    functionality to determine the names of the output files after alignment
    runs. Downstream rules refer to this function to specify the translations
    for a given segment.
    """
    checkpoint_output = checkpoints.nextclade.get(**wildcards).output.translations
    return expand("output/nextclade/flu/{lineage}/{segment}/nextclade.cds_translation.{gene}.fasta",
                  lineage=wildcards.lineage,
                  segment=wildcards.segment,
                  gene=GENES[wildcards.segment])


# Look for glycosylation sites per sequences
# Create a positions table per sequence, save to CSV..
rule get_glycosylation:
    input:
        aa_fasta = aggregate_translations
        # fasta_translated = "output/nextclade/flu/{lineage}/{segment}/nextclade.cds_translation.{gene}.fasta"
    output:
        glycosylation_csv = "output/glycosylation/{lineage}/{segment}/glycosylation_{lineage}_{gene}.csv"
    script:
        "library/scripts/get_glycolysation.py"

# Create a positions table per sequence, save to CSV..
rule get_positions:
    input:
        aa_fasta = aggregate_translations,
        positions_table = "input/aa_positions_of_interest.xlsx"
    output:
        positions_csv = "output/positions/{lineage}/{segment}/positions_{lineage}_{gene}.csv"
    script:
        "library/scripts/get_positions.py"

# Merge clade-data from nextstrain and add it to the positions data.
rule merge_positions_nextclade:
    input:
        nextclade_csv = rules.nextclade.output.nextclade_csv,  # "output/nextclade/flu/{lineage}/{segment}/nextclade.csv",
        positions_csv = rules.get_positions.output.positions_csv  # "output/positions/{lineage}/{segment}/positions_{lineage}_{gene}.csv"
    output:
        merged_csv = "output/positions/{lineage}/{segment}/positions_nextclade_{lineage}_{gene}.csv"
    script:
        "library/scripts/merge_positions_nextclade_csv.py"

# Merge clade-data from nextstrain and add it to the glycosylation data.
rule merge_glycosylation_nextclade:
    input:
        nextclade_csv = rules.nextclade.output.nextclade_csv,  # "output/nextclade/flu/{lineage}/{segment}/nextclade.csv",
        glycosylation_csv = rules.get_glycosylation.output.glycosylation_csv # "output/glycosylation/{lineage}/{segment}/glycosylation_{lineage}_{gene}.csv"
    output:
        merged_csv = "output/glycosylation/{lineage}/{segment}/glycosylation_nextclade_{lineage}_{gene}.csv"
    script:
        "library/scripts/merge_glycosylation_nextclade_csv.py"

rule merge_positions_nextclade_glycosylation:
    input: 
        positions_merged_csv = rules.merge_positions_nextclade.output.merged_csv, # "output/positions/{lineage}/{segment}/positions_nextclade_{lineage}_{gene}.csv"
        glycosylation_merged_csv = rules.merge_glycosylation_nextclade.output.merged_csv # "output/glycosylation/{lineage}/{segment}/glycosylation_nextclade_{lineage}_{gene}.csv"
    output:
        merged_csv = "output/positions/{lineage}/{segment}/positions_nextclade_glycosylation_{lineage}_{gene}.csv"
    script:
        "library/scripts/merge_positions_nextclade_glycosylation.py"