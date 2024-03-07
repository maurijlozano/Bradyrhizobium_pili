#!/usr/bin/env python3
import glob, sys, os, re
from Bio import SeqIO

if len(sys.argv) != 4:
    sys.exit(f'Usage: {sys.argv[0]} [folder] [genome_file_ext <gb>] [output <fasta>]')

folder = sys.argv[1]
ext = sys.argv[2]
fasta_file = sys.argv[3]

if not os.path.exists(folder):
    sys.exit(f'Error: {folder} not found..')
files = glob.glob(os.path.join(folder,f'*.{ext}'))
if len(files) == 0:
    sys.exit(f'Error: There are no files with extension {ext} in {folder} ..')

with open(fasta_file, 'w+') as f:
    for file in files:
        name_from_file = os.path.splitext(os.path.basename(file))[0]
        name_from_file = re.sub(' ','',name_from_file)
        for r in SeqIO.parse(file,'gb'):
            id = r.id
            description = re.sub(' ','_',r.description.split(',')[0])
            f.write(f'>{name_from_file}|{description}|{id}\n{r.seq}\n')
