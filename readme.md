# Analysis of _Bradyrhizobium_ type IV pili operons

## Programs used for the bioinformatic analysis
Here we describe the scripts used in the manuscripts for the identification of pili clusters in _Bradyrhizobium_.

### Requirements
* NCBI blast and dataset programs
* python > 3.7
* Python pandas library

The programs should be run in the following order.

### gb_files_to_fasta.py
This script extract fasta sequences from a set of genbank files located in a specific folder and generates a single fasta file containing all the genomes sequences to use as subject with blast_contigous_segments.py script. Only required if the genomes were downloaded in the genbank 'gbff' format (Recommended). 

Usage: ./gb_files_to_fasta.py [folder] [genome_file_ext <gb>] [output <fasta>]


### blast_contigous_segments.py
Given a DNA sequence it makes a blastn search and merges neighbor hits to extend the homology region. The objective is to obtain identify clusters which may contain gene insertions causing discontinuous blast hits for a single cluster.
This script requires a blast database made with makeblastdb program or a fasta file with all the concatenated genomes to use as subject (Can be generated from the gbff files with gb_files_to_fasta.py). Names in the fasta file (and database ids) should contain three fields separated by '|': 1) file name or genus_species_strai, 2) a custom descriptions, 3) replicon/scaffold/contig accession number. e.g: >Bradyrhizobium_diazoefficiens_USDA110|Bradyrhizobium diazoefficiens USDA110 chromosome|BA000040.2

```
Usage: ./blast_contiguous_segments.py [query <fasta>] [subject <fasta>] [output file name] join_gap[10000] hsp_qcov[70] p_ident[70] [-db]

e.g.1: ./blast_contiguous_segments.py query.fasta concatenated_genomes.fasta 2000 70 60

e.g.2: ./blast_contiguous_segments.py query.fasta database 2000 70 60 -db
       where database is the basename of a blast formated database created by makeblastdb. Allows blast to run in multiple threads.

join_gap: máximum number of nucleotides between discontinuous blast hits to merge
hsp_qcov: minimum query hsp coverage, percent
p_ident: minimum percent identity 
db: sets on the option to use a NCBI formated database
```

### make_clinker_plot.py
A python script that generates genbank files for all the identified segments and pipes them to clinker (Requires clinker installed and added to the PATH).

Usage: ./make_clinker_plot.py [table <csv file with cluster file_name record_id start end>] [cluster_index <int>] [gb file folder <complete path>]

## Accessory scripts
## download_gb_from_tsv.py
A script that uses NCBI datasets (required) app to download genomes taking as input a TSV file with the refseq genome accession numbers.

Usage: ./download_gb_from_tsv.py <tsv table file with headers> <accession number field name> <output folder>
