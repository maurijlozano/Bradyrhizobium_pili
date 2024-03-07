# Analysis of _Bradyrhizobium_ type IV pili operons

## Programs used for the bioinformatic analysis
Here we describe the scripts used in the manuscripts for the identification of pili clusters in _Bradyrhizobium_.

### Requirements
* NCBI blast and dataset programs
* python > 3.7
* Python pandas library

The programs should be run in the following order.

### gb_files_to_fasta.py
This script extract fasta sequences from a set of genbank files located in a specific folder and generates a single fasta file containing all the genomes sequences to use as subject with blast_contigous_segments.py script. Only required if the genomes were downloaded in the genbank 'gbff' format. 

Usage: ./gb_files_to_fasta.py [folder] [genome_file_ext <gb>] [output <fasta>]


### blast_contigous_segments.py
Given a DNA sequence it makes a blastn search and merges neighbor hits to extend the homology region.

Usage: ./blast_contigous_segments.py [query <fasta>] [subject <fasta>] [output file name]

### make_clinker_plot.py
A python script that generates genbank files for all the identified segments and pipes them clinker.

Usage: ./make_clinker_plot.py [table <csv file with cluster file_name record_id start end>] [cluster_index <int>] [gb file folder <complete path>]

## Accessory scripts
## download_gb_from_tsv.py
A script that uses NCBI datasets (required) app to download genomes taking as input a TSV file with the refseq genome accession numbers.

Usage: ./download_gb_from_tsv.py <tsv table file with headers> <accession number field name> <output folder>
