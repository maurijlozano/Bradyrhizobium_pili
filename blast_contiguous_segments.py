#!/usr/bin/env python3
import sys, os
import pandas as pd

def get_subject_names(subject):
    subject_ids = []
    with open(subject,'r') as f:
        for line in f:
            if '>' in line:
                line = line.strip()
                subject_ids.append(line[1:])
    return(subject_ids)

def get_subject_names_DB(subject):
    subject_ids = []
    os.system(f'blastdbcmd -db {subject} -entry all -out sid.txt -outfmt "%t"')
    with open("sid.txt",'r') as f:
        for line in f:
            line = line.strip()
            subject_ids.append(line)
    return(subject_ids)

def join_continuous_blast_hits(sseq_table, surrounds):
    contiguous_hits_table = pd.DataFrame(columns=sseq_table.columns)
    sseq_table = sseq_table.sort_values(['ab_start']).reset_index(drop=True)
    ab_start = sseq_table.iloc[0,:]['ab_start']
    ab_end = sseq_table.iloc[0,:]['ab_end']
    contiguous_hits_table.loc[0] = sseq_table.loc[0]
    contiguous_hits_table.loc[0,'evalue_list'] =  sseq_table.loc[0,'evalue']
    contiguous_hits_table.loc[0,'pid_list'] =  sseq_table.loc[0,'pident']
    contiguous_hits_table.loc[0,'segments'] =  1
    j = 0
    for i,row in sseq_table.iterrows():
        included_in_previous_hit = range(ab_start,ab_end)
        #surrounds make dependent on hit len, with max == surrounds
        # hit_len = row['length']
        # if surrounds/hit_len >= 2:
        #     surrounds_u = hit_len*2
        # else:
        #     surrounds_u = surrounds
        # if surrounds_u < 4000:
        #     surrounds_u = 4000
        start_range = range(ab_start-surrounds,ab_end+surrounds)
        if i > 0:
            if (row['ab_start'] in included_in_previous_hit) and (row['ab_end'] in included_in_previous_hit):
                continue
            elif row['ab_start'] in start_range:
                contiguous_hits_table.loc[j,'ab_end'] = row['ab_end']
                contiguous_hits_table.loc[j,'qstart'] = min(contiguous_hits_table.loc[j,'qstart'], row['qstart'])
                contiguous_hits_table.loc[j,'qend'] = max(contiguous_hits_table.loc[j,'qend'], row['qend'])
                contiguous_hits_table.loc[j,'hsp_qcov'] = contiguous_hits_table.loc[j,'hsp_qcov'] + row['hsp_qcov']
                contiguous_hits_table.loc[j,'length'] = contiguous_hits_table.loc[j,'length'] + row['length']
                contiguous_hits_table.loc[j,'sstart'] = contiguous_hits_table.loc[j,'ab_start']
                contiguous_hits_table.loc[j,'send'] = contiguous_hits_table.loc[j,'ab_end']
                contiguous_hits_table.loc[j,'evalue_list'] = f"{contiguous_hits_table.loc[j,'evalue_list']}, {row['evalue']}"
                contiguous_hits_table.loc[j,'pid_list'] = f"{contiguous_hits_table.loc[j,'pid_list']}, {row['pident']}"
                contiguous_hits_table.loc[j,'evalue'] = min(contiguous_hits_table.loc[j,'evalue'], row['evalue'])
                contiguous_hits_table.loc[j,'pident'] = min(contiguous_hits_table.loc[j,'pident'], row['pident'])
                contiguous_hits_table.loc[j,'bitscore'] = contiguous_hits_table.loc[j,'bitscore'] + row['bitscore']
                contiguous_hits_table.loc[j,'segments'] += 1
                ab_start = sseq_table.iloc[i,:]['ab_start']
                ab_end = sseq_table.iloc[i,:]['ab_end']
            else:
                j+=1
                contiguous_hits_table.loc[j] = sseq_table.loc[i]
                contiguous_hits_table.loc[j,'evalue_list'] =  sseq_table.loc[0,'evalue']
                contiguous_hits_table.loc[j,'pid_list'] =  sseq_table.loc[0,'pident']
                contiguous_hits_table.loc[j,'segments'] =  1
                ab_start = sseq_table.iloc[i,:]['ab_start']
                ab_end = sseq_table.iloc[i,:]['ab_end']
    return(contiguous_hits_table)

def removeOverlapping(tab):
    tab = tab.sort_values(['sseqid','srecid','ab_start'])
    for i in range(len(tab)-1):
        if (not tab['missing'].iloc[i]) and (not tab['missing'].iloc[i+1]):
            if (tab['srecid'].iloc[i] == tab['srecid'].iloc[i+1]) & ((tab['ab_start'].iloc[i+1] in range(int(tab['ab_start'].iloc[i]),int(tab['ab_end'].iloc[i]))) | (tab['ab_start'].iloc[i] in range(int(tab['ab_start'].iloc[i+1]),int(tab['ab_end'].iloc[i+1]))) ):
                if tab['ab_end'].iloc[i]-tab['ab_start'].iloc[i] < tab['ab_end'].iloc[i+1]-tab['ab_start'].iloc[i+1]:
                    tab.iloc[i,2:] = ''
                    tab['missing'].iloc[i] = True
                else:
                    tab.iloc[i+1,2:] = ''
                    tab['missing'].iloc[i+1] = True
    tab = tab.sort_values(['sseqid','qseqid'])
    return(tab)

if len(sys.argv) < 7 :
    print('Missing arguments --> join_gap[10000] hsp_qcov[70] p_ident[70] are not optional.')
    sys.exit(f'Usage: {sys.argv[0]} [query <fasta>] [subject <fasta>] [output file name] join_gap[10000] hsp_qcov[70] p_ident[70] [-db]')
elif len(sys.argv) == 8 :
    print('Using blast database. -db option')

useBlastDB = False
query = sys.argv[1]
subject = sys.argv[2]

if len(sys.argv) == 8 :
    db = sys.argv[7]
    if db == '-db':
        print(f'DB: {subject}')
        useBlastDB = True

if (not os.path.exists(query)):
    sys.exit(f'{query} missing...')
if (not os.path.exists(subject)) and (not useBlastDB):
    sys.exit(f'{subject} missing...')

table_file = sys.argv[3]

#params
surrounds = int(sys.argv[4])
hsp_qcov = int(sys.argv[5])
pident = int(sys.argv[6])

#
if (not os.path.exists(table_file)) and (not useBlastDB):
    try:
        os.system(f'blastn -task blastn -subject {subject} -query {query} -out {table_file} -perc_identity {pident} -evalue {0.001} -max_target_seqs 10000 -outfmt "6 qseqid sseqid length qlen qstart qend slen sstart send evalue bitscore pident qcovs"')
    except:
        print('Something went wrong, check your files.')
elif (not os.path.exists(table_file)) and (useBlastDB):
    try:
        os.system(f'blastn -task blastn -db {subject} -query {query} -out {table_file} -perc_identity {pident} -evalue {0.001} -max_target_seqs 10000 -outfmt "6 qseqid sseqid length qlen qstart qend slen sstart send evalue bitscore pident qcovs" -num_threads 16')
    except:
        print('Something went wrong, check your files.')
else:
    print(f'--> {table_file} exists, reading as previously run blast results...')

#read table. Blast trimms names in spaces....
table = pd.read_table(table_file, names=['qseqid','sseqid','length','qlen','qstart','qend','slen','sstart','send','evalue','bitscore','pident','qcovs'])
#table = table[table['evalue'] < 0.001]

#format table
table.loc[:,'hsp_qcov'] = [ ((row[5]-row[4])/row[3])*100 for _,row in table.iterrows()]
table.loc[:,'ab_start'] = [ row[7] if (row[7] < row[8]) else row[8] for _,row in table.iterrows()]
table.loc[:,'ab_end'] = [ row[8] if (row[7] < row[8]) else row[7] for _,row in table.iterrows()]
table.loc[:,'strand'] = [ '+' if (row[7] < row[8]) else '-' for _,row in table.iterrows()]

#
if not useBlastDB:
    subject_ids = get_subject_names(subject)
else:
    subject_ids = get_subject_names_DB(subject)

subject_ids = set([id.split('|')[0] for id in subject_ids])
 
continuous_blast_hits_table = pd.DataFrame()
best_hits_per_query = pd.DataFrame()
clusters = table.qseqid.unique()

for cluster in clusters:
    #
    c_continuous_blast_hits_table = pd.DataFrame()
    c_best_hits_per_subject = pd.DataFrame()
    clust_table = table[table['qseqid']== cluster]
    sseqids = clust_table.sseqid.unique() 
    for sseqid in sseqids:
        sseq_table = clust_table[clust_table['sseqid']== sseqid]
        continuous_blast_hits = join_continuous_blast_hits(sseq_table, surrounds)
        c_continuous_blast_hits_table = pd.concat([c_continuous_blast_hits_table,continuous_blast_hits], ignore_index=True)
        #get the best
        continuous_blast_hits_filtered = continuous_blast_hits[(continuous_blast_hits['hsp_qcov'] > hsp_qcov) & (continuous_blast_hits['pident'] > pident)]
        c_best_hits_per_subject = pd.concat([c_best_hits_per_subject,continuous_blast_hits_filtered], ignore_index=True)
    #only one best for genome
    c_best_hits_per_subject.loc[:,'srecid'] = [ hit.split('|')[2] for hit in c_best_hits_per_subject.sseqid]
    c_best_hits_per_subject['sseqid'] = [ hit.split('|')[0] for hit in c_best_hits_per_subject.sseqid]
    c_best_hits_per_subject_unique = pd.DataFrame(columns= c_best_hits_per_subject.columns)
    best_ids = c_best_hits_per_subject.sseqid.unique()
    k=0 
    for best_id in best_ids:
        sseq_table = c_best_hits_per_subject[c_best_hits_per_subject['sseqid']== best_id]
        best_id_hit = sseq_table.sort_values(['hsp_qcov'], ascending=False).reset_index(drop=True).loc[0,:]
        c_best_hits_per_subject_unique.loc[k,:] = best_id_hit
        k += 1
    #detect missing
    c_continuous_blast_hits_table['sseqid'] = [ hit.split('|')[0] for hit in c_continuous_blast_hits_table.sseqid]
    contiguous_sids = c_continuous_blast_hits_table.sseqid.unique()
    contiguous_sids = set([id.split('|')[0] for id in contiguous_sids])
    best_contiguous_sids = c_best_hits_per_subject_unique.sseqid.unique()
    missing_all_list = [ id for id in subject_ids if id not in contiguous_sids ]
    missing_best_list = [ id for id in subject_ids if id not in best_contiguous_sids ]
    #add missing to table
    c_continuous_blast_hits_table.loc[:,'missing'] = 'False'
    j=0
    for i in range(len(c_continuous_blast_hits_table),(len(c_continuous_blast_hits_table) + len(missing_all_list))):
        c_continuous_blast_hits_table.loc[i,'missing'] = 'True'
        c_continuous_blast_hits_table.loc[i,'qseqid'] = cluster
        c_continuous_blast_hits_table.loc[i,'sseqid'] = missing_all_list[j]
        j += 1
    #best
    if len(c_best_hits_per_subject_unique) > 0:
        c_best_hits_per_subject_unique.loc[:,'missing'] = 'False'
    else:
        c_best_hits_per_subject_unique = pd.DataFrame(columns=list(c_best_hits_per_subject.columns)+['missing'])
    j=0
    for i in range(len(c_best_hits_per_subject_unique),(len(c_best_hits_per_subject_unique) + len(missing_best_list))):
        c_best_hits_per_subject_unique.loc[i,'missing'] = 'True'
        c_best_hits_per_subject_unique.loc[i,'qseqid'] = cluster
        c_best_hits_per_subject_unique.loc[i,'sseqid'] = missing_best_list[j]
        j += 1
    #
    continuous_blast_hits_table = pd.concat([continuous_blast_hits_table,c_continuous_blast_hits_table], ignore_index=True)
    best_hits_per_query = pd.concat([best_hits_per_query,c_best_hits_per_subject_unique], ignore_index=True)
#

unfiltered_continuous_blast_hits = os.path.join(os.path.split(table_file)[0],'unfiltered_continuous_blast_hits.csv')
continuous_blast_hits_table.to_csv(unfiltered_continuous_blast_hits)

#analyze overlapping clusters
best_hits_per_query = removeOverlapping(best_hits_per_query)

filtered_continuous_blast_hits_qc = os.path.join(os.path.split(table_file)[0],f'filtered_continuous_blast_hits_qc{hsp_qcov}_pid{pident}.csv')
best_hits_per_query.to_csv(filtered_continuous_blast_hits_qc)

extract_for_clinker = os.path.join(os.path.split(table_file)[0],f'extract_for_clinker.csv')
best_hits_per_query[['qseqid','sseqid','srecid','ab_start','ab_end','strand']].to_csv(extract_for_clinker, index=False, header=False)
print('Task completed.')
