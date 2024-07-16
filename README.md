# Positions

### Snakemake file that provides information about AA positions in HA, NA and PA sequences.

Here is a tree view of how the 'input' folder should be organized and how sequence input files should be named:

input  
├── positions_by_lineage_and_segment.xlsx  
├── sequences_h1n1pdm_ha.fasta  
├── sequences_h1n1pdm_na.fasta  
├── sequences_h1n1pdm_pa.fasta  
├── sequences_h3n2_ha.fasta  
├── sequences_h3n2_na.fasta  
├── sequences_h3n2_pa.fasta  
├── sequences_h5n1_ha.fasta  
├── sequences_vic_ha.fasta  
├── sequences_vic_na.fasta  
└── sequences_vic_pa.fasta  


## Snakefile
Initiate by running:

```
Snakemake --dryrun  --cores 'all'
```

