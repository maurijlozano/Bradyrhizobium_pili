#!/usr/bin/env python3
#imports
import os, sys, glob
import pandas as pd
from Bio import SeqIO

#from table with regions to align make a clinker plot
print(f'Running {sys.argv[0]}...')
print(f'Checking arguments...')
if len(sys.argv) == 1:
    sys.exit(f'Usage: {sys.argv[0]} [table <csv file with cluster file_name record_id start end strand(+,-) >] [cluster_index <int>, -1 to process all clusters] [gb file folder <complete path>] expand[ <1/0> 1 Sets expand region ON]')

table_file = sys.argv[1]
cluster_index = int(sys.argv[2])-1

if cluster_index == -2:
    all_clusters = True
elif cluster_index < 0:
    sys.exit('Cluster index is 1 indexed.  The first cluster should be numbered "1" and corresponds to the first unique ID in the first column of the table...')
else:
    all_clusters = False

folder = sys.argv[3]
if folder == '':
    folder = './'
if not os.path.exists(folder):
    sys.exit(f'{folder} not found...')
else:
    gb_files = glob.glob(os.path.join(folder,'*.gb'))
    if len(gb_files) == 0:
        sys.exit('No gb files found...')

expand = int(sys.argv[4])
if not expand in [0,1]:
    sys.exit('Cluster index is 1 indexed.  The first cluster should be numbered "1" and corresponds to the first unique ID in the first column of the table...')

output_folder = 'regions'
if not os.path.exists(output_folder):
    print(f'Creating output folder: regions...')
    os.mkdir(output_folder)

table = pd.read_csv(table_file, header=None)

if len(table.columns) != 6:
    print('Error, incorrect table format...')
    sys.exit(f'Usage: {sys.argv[0]} [table <csv file with cluster file_name record_id start end strand(+,-)>] [cluster_index <int>] [gb file folder <complete path>]')

if not all_clusters:
    cluster_name = table.iloc[:,0].unique()[cluster_index]
    ctable = table[~table[2].isna()]
    ctable = ctable[ctable[0] == cluster_name]
else:
    ctable = table[~table[2].isna()]

#generate table file names. Extension = gb
ctable.loc[:,'files'] = [os.path.join(folder, f+'.gb') for f in ctable[1] ]
bigger_region = max([ abs(r[4]-r[3])  for _,r in ctable.iterrows() ]) * 2

#extract regions for clinker
print('Extracting regions from gb files...')
for _,reg in ctable.iterrows():
    file = reg.files
    name = reg[1]
    srecid = reg[2]
    start = int(reg[3])
    end = int(reg[4])
    strand = reg[5]
    #
    if expand:
        expand_region_len = int((bigger_region - abs(end-start)) /2)
        if expand_region_len < 5000:
            expand_region_len = 5000
        if (start-expand_region_len) > 0:
            cstart = start-expand_region_len
        else:
            cstart = 0
        cend = end + expand_region_len
    else:
        cstart = start
        cend = end
    #
    if not os.path.exists(file):
        print(f'File {file} not found, skipping...')
    else:
        region_file_name = os.path.join(output_folder,f'{name}_{srecid}_{start}_{end}.gb')
        print(f'Generating region file: {region_file_name}')
        if not os.path.exists(region_file_name):
            for r in SeqIO.parse(file,'gb'):
                if r.id == srecid:
                    reg_found = True
                    region_rec = r[cstart:cend]
                    if strand == '+': 
                        SeqIO.write(region_rec,region_file_name,'gb')
                    else:
                        region_rec_rc = region_rec.reverse_complement()
                        region_rec_rc.id = r.id
                        region_rec_rc.annotations = r.annotations
                        region_rec_rc.description = r.description
                        region_rec_rc.name = r.name
                        SeqIO.write(region_rec_rc,region_file_name,'gb')
                    break
            if not reg_found:
                print(f'Please verify table, {srecid} not found...')
        else:
            print(f'{region_file_name} found, skipping...')

#running clinker
print('> Running clinker...')
clinker_html_out = f'clinker_plot_cluster_{cluster_index}.html'
os.system(f'clinker {os.path.join(output_folder,"*.gb")} -p {clinker_html_out}')
print('Task done!')

