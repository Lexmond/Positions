from pathlib import Path

# Directory paths
input_dir = Path("input")
output_dir = "output/nextclade/flu"
dataset_dir = "library/nextclade/data/flu"

# Gene lookup dictionary
gene_lookup = {
    "pb2": "PB2",
    "pb1": "PB1",
    "pa": "PA",
    "ha": "HA1",
    "np": "NP",
    "na": "NA",
    "ma": "M1",
    "ns": "NS",
}

# Function to scan the input folder and generate lists of lineages, segments, and genes
def input_func_nextalign():
    l1, l2, l3 = [], [], []
    for path in list(input_dir.glob("*.fasta")):
        _, lineage, segment = path.stem.split('_')  # sequences_h1n1pdm_ha.fasta
        gene = gene_lookup[segment]  # Get gene variable from dict
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

# Create a list of tuples that can be used to extract zipped values from the different lists
list_of_tuples = list(zip(input_data_nextalign['lineage'], input_data_nextalign['segment']))

nextclade_output_files = [f"output/nextclade/flu/{lineage}/{segment}/nextclade.csv" for lineage, segment in list_of_tuples]
print(f"All output files for Nextclade: \n {nextclade_output_files} \n")

glycosylation_output_files = [f"output/nextclade/flu/{lineage}/ha/glycosylation_sites_HA1.csv" for lineage in input_data_nextalign['lineage']]
print(f"All output files for Glycosylation sites: \n {glycosylation_output_files} \n")

# Snakemake rules
rule all:
    input:
        nextclade_output_files,
        glycosylation_output_files

# Run Nextclade datasets on the input sequences to get outputs per subtype and segment 
rule run_nextclade:
    input:
        fasta = "input/sequences_{lineage}_{segment}.fasta",
        dataset = "library/nextclade/data/flu/{lineage}/{segment}"
    output:
        csv = "output/nextclade/flu/{lineage}/{segment}/nextclade.csv"
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-all=output/nextclade/flu/{wildcards.lineage}/{wildcards.segment} \
            {input.fasta}
        """

# Create a CSV file that contains all the glycosylation sites in HA1 from the input sequences
rule find_glycosylation_sites:
    input:
        fasta_translated_HA1 = "output/nextclade/flu/{lineage}/ha/nextclade.cds_translation.HA1.fasta"
    output:
        glycosylation_csv = "output/nextclade/flu/{lineage}/ha/glycosylation_sites_HA1.csv"
    script:
        "library/scripts/find_glycolysation_sites.py"

