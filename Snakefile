from pathlib import Path

# Directory 'input' path
INPUT = Path("input")

# Gene lookup dictionary
GENES = {
    "ha":  "HA1",
    "na":  "NA",
    "pa":  "PA",
}

# Function to scan the input folder and generate lists of lineages, segments, and genes
def input_func_nextalign():
    l1, l2, l3 = [], [], []
    for path in list(INPUT.glob("*.fasta")):
        _, lineage, segment = path.stem.split('_')  # sequences_h1n1pdm_ha.fasta
        gene = GENES[segment]  # Get gene variable from dict
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



# Create a list of output files for rule nextclade
nextclade_output_files = [f"output/nextclade/flu/{lineage}/{segment}/nextclade.csv" for lineage, segment in list_of_2tuples]
print(f"All output files for Nextclade: \n {nextclade_output_files} \n")

# Create a list of output files for rule glycosylation_sites
glycosylation_output_files = [f"output/glycosylation/flu/{lineage}/{segment}/glycosylation_sites_{gene}.csv" for lineage, segment, gene in list_of_3tuples]
print(f"All output files for Glycosylation sites: \n {glycosylation_output_files} \n")




# Snakemake rules
rule all:
    input:
        nextclade_output_files,
        glycosylation_output_files
        


# Run Nextclade datasets on the input sequences to get outputs per subtype and segment 
checkpoint nextclade:
    input:
        fasta = "input/sequences_{lineage}_{segment}.fasta",
        dataset = "library/nextclade/data/flu/{lineage}/{segment}"
    output:
        csv = "output/nextclade/flu/{lineage}/{segment}/nextclade.csv",
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




# intermediate rule in checkpoint...
# Create a CSV file that contains all the glycosylation sites in HA1 from the input sequences
rule glycosylation_sites:
    input:
        aggregate_translations
        # fasta_translated = "output/nextclade/flu/{lineage}/{segment}/nextclade.cds_translation.{gene}.fasta"
    output:
        glycosylation_csv = "output/glycosylation/flu/{lineage}/{segment}/glycosylation_sites_{gene}.csv"
    script:
        "library/scripts/find_glycolysation_sites.py"

