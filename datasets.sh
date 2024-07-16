

# H3N2 Datasets
nextclade dataset get --name 'nextstrain/flu/h3n2/pb2' --output-dir 'library/nextclade/data/flu/h3n2/pb2'
nextclade dataset get --name 'nextstrain/flu/h3n2/pb1' --output-dir 'library/nextclade/data/flu/h3n2/pb1'
nextclade dataset get --name 'nextstrain/flu/h3n2/pa' --output-dir 'library/nextclade/data/flu/h3n2/pa'
nextclade dataset get --name 'nextstrain/flu/h3n2/ha/EPI1857216' --output-dir 'library/nextclade/data/flu/h3n2/ha'
nextclade dataset get --name 'nextstrain/flu/h3n2/np' --output-dir 'library/nextclade/data/flu/h3n2/np'
nextclade dataset get --name 'nextstrain/flu/h3n2/na/EPI1857215' --output-dir 'library/nextclade/data/flu/h3n2/na'
nextclade dataset get --name 'nextstrain/flu/h3n2/mp' --output-dir 'library/nextclade/data/flu/h3n2/mp'
nextclade dataset get --name 'nextstrain/flu/h3n2/ns' --output-dir 'library/nextclade/data/flu/h3n2/ns'

# H1N1pdm Datasets
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pb2' --output-dir 'library/nextclade/data/flu/h1n1pdm/pb2'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pb1' --output-dir 'library/nextclade/data/flu/h1n1pdm/pb1' 
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pa' --output-dir 'library/nextclade/data/flu/h1n1pdm/pa'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/ha/MW626062' --output-dir 'library/nextclade/data/flu/h1n1pdm/ha'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/np' --output-dir 'library/nextclade/data/flu/h1n1pdm/np'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/na/MW626056' --output-dir 'library/nextclade/data/flu/h1n1pdm/na'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/mp' --output-dir 'library/nextclade/data/flu/h1n1pdm/mp'
nextclade dataset get --name 'nextstrain/flu/h1n1pdm/ns' --output-dir 'library/nextclade/data/flu/h1n1pdm/ns'

# VIC Datasets
nextclade dataset get --name 'nextstrain/flu/vic/ha/KX058884' --output-dir 'library/nextclade/data/flu/vic/ha'
nextclade dataset get --name 'nextstrain/flu/vic/na/CY073894' --output-dir 'library/nextclade/data/flu/vic/na'

# H5N1 Datasets
nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/all-clades' --output-dir 'library/nextclade/data/flu/h5n1/ha/all-clades'
nextclade dataset get --name 'community/moncla-lab/iav-h5/ha/2.3.4.4' --output-dir 'library/nextclade/data/flu/h5n1/ha/2.3.4.4'


#Sequences downloaded from Gisaid uploaded by EMC or RIVM in the period: 2020-Oct-01 -> 2024-Jul-14




#####----- A/H1N1pdm -----#####

# A/H1N1pdm HA
nextclade run \
   --input-dataset library/nextclade/data/flu/h1n1pdm/ha \
   --output-all=output/nextclade/flu/h1n1pdm/ha \
   input/sequences_h1n1pdm_ha.fasta
### NOTE: A/Netherlands/10534/2023|A_/_H1N1|unassigned|Original|EPI_ISL_18168180 is not H1N1pdm, but H1N1v...

# A/H1N1pdm NA
nextclade run \
   --input-dataset library/nextclade/data/flu/h1n1pdm/na \
   --output-all=output/nextclade/flu/h1n1pdm/na \
   input/sequences_h1n1pdm_na.fasta

# A/H1N1pdm PA
nextclade run \
   --input-dataset library/nextclade/data/flu/h1n1pdm/pa \
   --output-all=output/nextclade/flu/h1n1pdm/pa \
   input/sequences_h1n1pdm_pa.fasta






#####----- A/H3N2 -----#####

# A/H3N2 HA
nextclade run \
   --input-dataset library/nextclade/data/flu/h3n2/ha \
   --output-all=output/nextclade/flu/h3n2/ha \
   input/sequences_h3n2_ha.fasta

# A/H3N2 NA
nextclade run \
   --input-dataset library/nextclade/data/flu/h3n2/na \
   --output-all=output/nextclade/flu/h3n2/na \
   input/sequences_h3n2_na.fasta

# A/H3N2 PA
nextclade run \
   --input-dataset library/nextclade/data/flu/h3n2/pa \
   --output-all=output/nextclade/flu/h3n2/pa \
   input/sequences_h3n2_pa.fasta


#####----- B/Vic -----#####

# B/Vic HA
nextclade run \
   --input-dataset library/nextclade/data/flu/vic/ha \
   --output-all=output/nextclade/flu/vic/ha \
   input/sequences_vic_ha.fasta

# B/Vic NA
nextclade run \
   --input-dataset library/nextclade/data/flu/vic/na \
   --output-all=output/nextclade/flu/vic/na \
   input/sequences_vic_na.fasta

# B/Vic PA (Custom made)
nextclade run \
   --input-dataset library/nextclade/data/flu/vic/pa \
   --output-all=output/nextclade/flu/vic/pa \
   input/sequences_vic_pa.fasta




#####----- A/H5N1 -----#####

# A/H5N1 (all-clades)
nextclade run \
   --input-dataset library/nextclade/data/flu/h5n1/ha/all-clades \
   --output-all=output/nextclade/flu/h5n1/ha/all-clades \
   input/sequences_h5n1_ha.fasta

# A/H5N1 (clade 2.3.4.4)
nextclade run \
   --input-dataset library/nextclade/data/flu/h5n1/ha/2.3.4.4 \
   --output-all=output/nextclade/flu/h5n1/ha/2.3.4.4 \
   input/sequences_h5n1_ha.fasta