# UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for fast functional annotation and visualization of (meta)-omics datasets 

This repository was created to store the scripts used in the publication "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for fast functional annotation and visualization of (meta)-omics datasets 
". These scripts run the following steps:
1. [Obtention of datasets](https://github.com/iquasere/functional_annotation_publication#1.-obtention-of-datasets)
2. [Installation of tools](https://github.com/iquasere/functional_annotation_publication#2.-installation-of-tools)
3. [Execution of tools](https://github.com/iquasere/functional_annotation_publication#3.-execution-of-tools)
4. [Results analysis](https://github.com/iquasere/functional_annotation_publication#4.-results-analysis)

## 1. Obtention of datasets

```
git clone https://github.com/iquasere/functional_annotation_publication.git
mkdir ann_paper
wget -i functional_annotation_publication/genomes_links.txt -P ann_paper
gunzip ann_paper/*
awk '{print $1}' ann_paper/* >> ann_paper/genomes.fa
```

## 2. Installation of tools

```
mamba create -c conda-forge -c bioconda upimapi recognizer keggcharter dfast prokka eggnog-mapper -n ann_paper
mkdir ~/anaconda3/envs/ann_paper/lib/python3.9/site-packages/data
download_eggnog_data.py

git clone https://github.com/PedroMTQ/mantis.git
cd mantis
conda env create -f mantis_env.yml
conda activate mantis_env
cd ..
python mantis setup_databases
```

## 3. Execution of tools

```
upimapi.py -i ann_paper/proteomes.fasta -o ann_paper/upimapi_proteomes -rd resources_directory --evalue 0.1 -mts 20 -db uniprot -t 15
recognizer.py -f ann_paper/proteomes.fasta -o ann_paper/recognizer_proteomes -rd resources_directory --evalue 0.1 -mts 20 -t 15 --tax-file ann_paper/upimapi_proteomes/UPIMAPI_results.tsv --tax-col "Taxonomic lineage (SPECIES)" --protein-id-col qseqid
emapper.py -i ann_paper/proteomes.fasta -o ann_paper/eggnog_mapper_proteomes/eggnog_results --cpu 15
python mantis run_mantis -t ann_paper/proteomes.fasta -o ann_paper/mantis_proteomes -c 15
```

## 4. Results analysis

```python
out = 'ann_paper'

# Calculation of best e-values
def is_true_positive(l1, l2, threshold=0.5):
    l1, l2 = set(l1), set(l2)       # must be set because some proteins will have multiple times the same identification attributed because of multiple positions of domain
    return sum([i1 in l2 for i1 in l1]) / len(l1) > threshold


def calculate_quality_metrics(tp, fp, fn):
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    return precision, recall, f1_score


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def parse_mantis_consensus(filename):
    file = open(filename).readlines()
    result = []
    for line in file[1:]:
        line = line.rstrip('\n').split('\t')
        result.append(line[:5] + [';'.join(line[6:])])
    df = pd.DataFrame(result, columns=['Query', 'Ref_Files', 'Ref_Hits', 'Consensus_hits', 'Total_hits', 'Links'])
    df = pd.concat([df, pd.DataFrame.from_records(df['Links'].apply(parse_links))], axis=1)
    del df['Links']
    return df


evalues = [1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-18, 1e-21, 1e-24, 1e-27, 1e-30]

# For UPIMAPI
upimapi_meta = pd.read_csv(f'{out}/upimapi_proteomes/UPIMAPI_results.tsv', sep='\t')[['qseqid', 'sseqid','evalue']]
upimapi_meta['qseqid'] = upimapi_meta['qseqid'].apply(lambda x: x.split('|')[1])
upimapi_meta = upimapi_meta.groupby('qseqid')['sseqid', 'evalue'].first().reset_index()
upimapi_meta = upimapi_meta[upimapi_meta['sseqid'] != '*']
n_proteins = count_on_file('>', f'{out}/proteomes.fasta')
upimapi_metrics = []
for evalue in evalues:
    print(evalue)
    upi_meta = upimapi_meta[upimapi_meta['evalue'] < evalue]
    tp = (upi_meta['qseqid'] == upi_meta['sseqid']).sum()
    fp = (upi_meta['qseqid'] != upi_meta['sseqid']).sum()
    fn = n_proteins - len(upi_meta)      # proteins that could not be identified
    precision, recall, f1_score = calculate_quality_metrics(tp, fp, fn)
    upimapi_metrics.append({'evalue': evalue, 'TPs': tp, 'FPs': fp, 'FNs': fn, 'precision': precision,
                            'recall': recall, 'f1_score': f1_score})
upimapi_metrics = pd.DataFrame(upimapi_metrics)
upimapi_metrics.to_excel(f'{out}/upimapi_metrics.xlsx')

# For reCOGnizer
uniprotinfo_meta = pd.read_csv(f'{out}/upimapi_proteomes_uniprotinfo/uniprotinfo.tsv', sep='\t')

def analyze_recognizer(recognizer_excel_filename, metrics_filename):
    recognizer_meta = {}
    for db in tqdm(['CDD', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM', 'Pfam', 'Smart', 'COG', 'KOG'],
                   desc='Reading DBs DFs'):
        if db == 'KOG':
            recognizer_meta[db] = pd.DataFrame(columns=[
                'qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                'qend', 'sstart', 'send', 'evalue', 'bitscore', 'DB ID', 'DB description', 'sequence', 'cog',
                'COG functional category (letter)', 'Protein description', 'COG general functional category',
                'COG functional category'])        # necessary because KOG had no matches
        else:
            df = pd.read_excel(recognizer_excel_filename, sheet_name=db).drop_duplicates()
            df['qseqid'] = df['qseqid'].apply(lambda x: x.split('|')[1])
            df = df.groupby('qseqid')[df.columns.tolist()[1:]].first().reset_index()
            recognizer_meta[db] = pd.merge(df, uniprotinfo_meta, left_on='qseqid', right_on='Entry', how='left')
    recog_to_sp_dbs = {  # DB: (db in uniprot columns, prefix of smps)
        'CDD': ('Cross-references (CDD)', 'cd'),
        'Pfam': ('Cross-references (Pfam)', 'pfam'),
        'NCBIfam': (None, 'NF'),
        'Protein_Clusters': (None, 'PRK'),
        'TIGRFAM': ('Cross-references (TIGRFAMs)', 'TIGR'),
        'Smart': ('Cross-references (SMART)', 'smart'),
        'COG': ('Cross-references (eggNOG)', 'COG'),
        'KOG': ('Cross-references (eggNOG)', 'KOG')}
    recognizer_metrics = []
    for evalue in evalues:
        print(evalue)
        for db, df in recognizer_meta.items():
            print(db)
            v = recog_to_sp_dbs[db]
            if v[0] is not None and len(df) > 0:
                df = df[(df['evalue'].notnull()) & (df['evalue'] < evalue)]     # some proteins are expelled from analysis because blast metrics could not be obtained
                grouped = df[(df['DB ID'].notnull()) & (df[v[0]].notnull())].groupby('qseqid')['DB ID'].apply(
                    ';'.join).reset_index()
                del df['DB ID']
                df = pd.merge(grouped, df, on='qseqid', how='left')
                if db in ['COG', 'KOG']:
                    df = df[df[v[0]].str.startswith(v[1])]
                elif db == 'Pfam':
                    df['DB ID'] = df['DB ID'].str.replace('pfam', 'PF')
                elif db == 'Smart':
                    df['DB ID'] = df['DB ID'].str.replace('smart', 'SM')
                df[v[0]] = df[v[0]].apply(lambda x: x[:-1].split(';'))
                df['DB ID'] = df['DB ID'].apply(lambda x: x.split(';'))
                tps = sum([is_true_positive(df.iloc[i]['DB ID'], df.iloc[i][v[0]]) for i in range(len(df))])
                fps = len(df) - tps
                fns = uniprotinfo_meta[v[0]].notnull().sum() - tps
                print(f'TPs:{tps} FPs:{fps} FNs:{fns}')
                precision, recall, f1_score = calculate_quality_metrics(tps, fps, fns)
                recognizer_metrics.append(
                    {'evalue': evalue, 'db': db, 'TPs': tps, 'FPs': fps, 'FNs': fns, 'precision': precision,
                     'recall': recall, 'f1_score': f1_score})
    recognizer_metrics = pd.DataFrame(recognizer_metrics)
    recognizer_metrics.sort_values(by=['db', 'evalue'], ascending=False)[
        ['db', 'evalue'] + recognizer_metrics.columns.tolist()[2:]].to_excel('sorted_metrics.xlsx', index=False)
    recognizer_metrics.to_excel(metrics_filename, index=False)
    return recognizer_metrics

recognizer_analysis = analyze_recognizer(
    f'{out}/recognizer_proteomes/reCOGnizer_results.xlsx',
    f'{out}/recognizer_proteomes/recognizer_metrics.xlsx')
recognizer_analysis_without_taxonomy = analyze_recognizer(
    f'{out}/recognizer_proteomes_without_taxonomy/reCOGnizer_results.xlsx',
    f'{out}/recognizer_proteomes_without_taxonomy/recognizer_metrics.xlsx')


'''
eggNOG mapper
'''
egg_basename = f'{out}/eggnog_mapper_proteomes/eggnog_results.emapper'
eggnog_annotations = pd.read_csv(f'{egg_basename}.annotations', sep='\t', skiprows=4, skipfooter=3)


'''
mantis
'''
mantis_basename = f'{out}/mantis_proteomes'
mantis_consensus = parse_mantis_consensus(f'{mantis_basename}/consensus_annotation.tsv')
mantis_integrated = pd.read_csv(f'{mantis_basename}/integrated_annotation.tsv', sep='\t')
mantis_output = pd.read_csv(f'{mantis_basename}/output_annotation.tsv', sep='\t')


"""
Alignment of simulated MT reads to proteomes
"""

def run_command(bashCommand, output='', mode='w', sep=' ', print_message=True, verbose=True):
    if print_message:
        print(f"{bashCommand.replace(sep, ' ')}{' > ' + output if output != '' else ''}")
    if output == '':
        subprocess.run(bashCommand.split(sep), stdout=sys.stdout if verbose else None, check=True)
    else:
        with open(output, mode) as output_file:
            subprocess.run(bashCommand.split(sep), stdout=output_file)


def run_pipe_command(bashCommand, output='', mode='w', sep=' ', print_message=True):
    if print_message:
        print(bashCommand)
    if output == '':
        subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True,
                                stdout=subprocess.PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True, stdout=output_file).communicate()


def check_bowtie2_index(index_prefix):
    files = glob.glob(index_prefix + '*.bt2')
    if len(files) < 6:
        return False
    return True


def generate_mg_index(reference, index_prefix):
    run_command(f'bowtie2-build {reference} {index_prefix}')


def align_reads(reads, index_prefix, sam, report, log=None, threads=6):
    run_command(f'bowtie2 -x {index_prefix} -1 {reads[0]} -2 {reads[1]} -S {sam} -p {threads} 1> {report} 2> {log}')


def perform_alignment(reference, reads, basename, threads=1):
    ext = f".{reference.split('.')[-1]}"

    if not check_bowtie2_index(reference.replace(ext, '_index')):
        print('INDEX files not found. Generating new ones')
        generate_mg_index(reference, reference.replace(ext, '_index'))
    else:
        print(f"INDEX was located at {reference.replace(ext, '_index')}")

    if not os.path.isfile(f'{basename}.log'):
        align_reads(reads, reference.replace(ext, '_index'), f'{basename}.sam',
                    f'{basename}_bowtie2_report.txt', log=f'{basename}.log', threads=threads)
    else:
        print(f'{basename}.log was found!')

    run_pipe_command(
        f"""samtools view -F 260 {basename}.sam | cut -f 3 | sort | uniq -c | """
        f"""awk '{{printf("%s\\t%s\\n", $2, $1)}}'""", output=f'{basename}.readcounts')


'''
MT quantification
'''
for l in ['a', 'b', 'c']:
    for n in [0.01, 1, 100]:
        perform_alignment(f'{out}/transcriptomes.fasta', [
            f'SimMOSCA2/rna/{l}/{n}/mt_{n}{l}_R{fr}.fastq' for fr in ['1', '2']],
            f'{out}/alignments/mt_{n}{l}', threads=15)

readcounts = pd.DataFrame(columns=['name'])
for l in ['a', 'b', 'c']:
    for n in [0.01, 1, 100]:
        readcounts = pd.merge(
            readcounts, pd.read_csv(f'{out}/alignments/mt_{n}{l}.readcounts', sep='\t', header=None,
                                    names=['name', f'mt_{n}{l}']), how='outer', on='name')
mt_cols = [f'mt_{n}{l}' for l in ['a', 'b', 'c'] for n in [0.01, 1, 100]]
readcounts[mt_cols] = readcounts[mt_cols].fillna(value=0.0).astype(int)
genbank2uniprot = pd.read_csv(f'{out}/id_conversion.tsv', sep='\t', usecols=[0, 1], index_col=0)
present = readcounts[readcounts['name'].isin(genbank2uniprot.index)]
present['name'] = present['name'].apply(lambda x: genbank2uniprot.loc[x, 'Entry'])
upimapi_meta = pd.read_csv(f'{out}/upimapi_proteomes/UPIMAPI_results.tsv', sep='\t')
upimapi_meta = upimapi_meta[upimapi_meta['sseqid']!='*']
upimapi_meta['qseqid'] = upimapi_meta['qseqid'].apply(lambda x: x.split('|')[1])
upimapi_meta = upimapi_meta.groupby('qseqid')[upimapi_meta.columns.tolist()[1:]].first().reset_index()
recognizer_proteome = pd.read_csv(f'{out}/recognizer_proteomes/reCOGnizer_results.tsv', sep='\t')
keggcharter_input = pd.merge(upimapi_meta, present, left_on='qseqid', right_on='name', how='left')
keggcharter_input['mg'] = [1] * len(keggcharter_input)
keggcharter_input.to_csv(f'{out}/keggcharter_input.tsv', sep='\t', index=False)
'''
keggcharter.py -o ann_paper/keggcharter_proteome -f ann_paper/keggcharter_input.tsv -rd resources_directory -tc "Taxonomic lineage (SPECIES)" -keggc "Cross-references (KEGG)" -tcol mt_0.01a,mt_1a,mt_100a,mt_0.01b,mt_1b','mt_100b,mt_0.01c,mt_1c,mt_100c -gcol mg
'''


def true_positive_funcs(l1, l2, threshold=0.5):
    return sum([i1 in l2 for i1 in l1]) / len(l1) > threshold


def parse_links(data):
    vals = data.split(';')
    result = {}
    for val in vals:
        pair = val.split(':')
        try:
            if pair[0] in result.keys():
                result[pair[0]] += f';{pair[1]}'
            else:
                result[pair[0]] = pair[1]
        except:
            pass
    return result


def parse_fasta_on_memory(file):
    lines = [line.rstrip('\n') for line in open(file)]
    i = 0
    result = dict()
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:].split()[0]
            result[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                result[name] += lines[i]
                i += 1
    return pd.DataFrame.from_dict(result, orient='index', columns=['sequence'])


upimapi_proteome = pd.read_csv(f'{out}/upimapi_proteomes/UPIMAPI_results.tsv', sep='\t')
recognizer_proteome = pd.read_csv(f'{out}/recognizer_proteomes/reCOGnizer_results.tsv', sep='\t')
recognizer_proteome['qseqid'] = recognizer_proteome['qseqid'].apply(lambda x: x.split('|')[1])
cog_ser = recognizer_proteome[recognizer_proteome['cog'].notnull()].groupby('qseqid')['cog'].apply(
    lambda x: '; '.join(set(x)))
ecs_ser = recognizer_proteome[recognizer_proteome['EC number'].notnull()].groupby('qseqid')['EC number'].apply(
    lambda x: '; '.join(set(x)))
eggnog_annotations['#query'] = eggnog_annotations['#query'].apply(lambda x: x.split('|')[1])
eggnog_annotations['eggNOG_OGs'] = eggnog_annotations['eggNOG_OGs'].apply(lambda x: x.split('@')[0])
mantis_consensus['Query'] = mantis_consensus['Query'].apply(lambda x: x.split('|')[1])
upimapi_meta = upimapi_meta.groupby('Entry')[['EC number', 'Cross-reference (eggNOG)']].first().reset_index()

# Build final analysis DF
res_df = parse_fasta_on_memory(f'{out}/proteomes.fasta').reset_index()
res_df.columns = ['Entry', 'sequence']
res_df['Entry'] = res_df['Entry'].apply(lambda x: x.split('|')[1])
res_df = pd.merge(res_df, uniprotinfo_meta[['Entry', 'EC number', 'Cross-reference (eggNOG)']], on='Entry', how='left')
res_df = pd.merge(res_df, upimapi_meta[['Entry', 'EC number', 'Cross-reference (eggNOG)']], on='Entry', how='left')
res_df = pd.merge(res_df, cog_ser, left_on='Entry', right_index=True, how='left')
res_df = pd.merge(res_df, ecs_ser, left_on='Entry', right_index=True, how='left')
res_df = pd.merge(res_df, eggnog_annotations[['#query', 'EC', 'eggNOG_OGs']], left_on='Entry', right_on='#query', how='left')
res_df = pd.merge(res_df, mantis_consensus[['Query', 'enzyme_ec', 'cog']], left_on='Entry', right_on='Query', how='left')
res_df.rename(columns={
        'EC number_x': 'EC number (UniProt)',
        'EC number_y': 'EC number (UPIMAPI)',
        'EC number': 'EC number (reCOGnizer)',
        'EC': 'EC number (eggNOG mapper)',
        'enzyme_ec': 'EC number (mantis)',
        'Cross-reference (eggNOG)_x': 'COG (UniProt)',
        'Cross-reference (eggNOG)_y': 'COG (UPIMAPI)',
        'cog_x': 'COG (reCOGnizer)',
        'eggNOG_OGs': 'COG (eggNOG mapper)',
        'cog_y': 'COG (mantis)'}, inplace=True)
for col in ['sequence', '#query', 'Query']:
    del res_df[col]
res_df.fillna(value='', inplace=True)
for col in res_df.columns.tolist():
    res_df[col] = res_df[col].str.replace(',', ';').str.replace(' ', '').str.rstrip(';').str.split(';')

res = []
for ide in ['EC number', 'COG']:
    df = res_df[~res_df[f'{ide} (UniProt)'].isin([['']])]
    print(len(df))
    for tool in ['UPIMAPI', 'reCOGnizer', 'eggNOG mapper', 'mantis']:
        print(f"Number of #s : {(df[f'{ide} ({tool})'].isin([['']])).sum()}")
        tps = sum([true_positive_funcs(df.iloc[i][f'{ide} (UniProt)'], df.iloc[i][f'{ide} ({tool})']) for i in range(len(df))])
        fps = (~df[f'{ide} ({tool})'].isin([['']])).sum() - tps
        fns = (df[f'{ide} ({tool})'].isin([['']])).sum()
        print(f'TPs:{tps} FPs:{fps} FNs:{fns}')
        precision, recall, f1_score = calculate_quality_metrics(tps, fps, fns)
        res.append({'Qualifier': ide, 'Tool': tool, 'TPs': tps, 'FPs': fps, 'FNs': fns, 'precision': precision,
                    'recall': recall, 'f1_score': f1_score})
pd.DataFrame(res).to_excel('tools_comparison.xlsx',index=False)

joined_df = res_df
for col in joined_df.columns.tolist():
    joined_df[col] = joined_df[col].apply(lambda x: ','.join(x))

lists = {}
for ide in ['EC number', 'COG']:
    for tool in ['UPIMAPI', 'reCOGnizer', 'eggNOG mapper', 'mantis']:
        print(f'{ide} ({tool}):{(joined_df[f"{ide} ({tool})"]!="").sum()}')
        #lists[f'{ide} ({tool})'] = set(','.join(set(joined_df[f'{ide} ({tool})'].tolist())).split(','))

ide = 'COG'

# 1
len(lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)'])
len(lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)'])
len(lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (mantis)'])
len(lists[f'{ide} (mantis)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (reCOGnizer)'])
# 2
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (UPIMAPI)'])) - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)'])
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (eggNOG mapper)'])) - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (mantis)'])
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)']) - len(lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (mantis)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (reCOGnizer)'])
len(set(list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (eggNOG mapper)'])) - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (mantis)'])
len(set(list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (eggNOG mapper)']) - len(lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (mantis)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (reCOGnizer)'])
len(set(list(lists[f'{ide} (eggNOG mapper)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (UPIMAPI)']) - len(lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (reCOGnizer)'] - lists[f'{ide} (mantis)']) - len(lists[f'{ide} (mantis)'] - lists[f'{ide} (UPIMAPI)'] - lists[f'{ide} (eggNOG mapper)'] - lists[f'{ide} (reCOGnizer)'])
# 3
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (eggNOG mapper)'])) - lists[f'{ide} (mantis)'])
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (eggNOG mapper)'])
len(set(list(lists[f'{ide} (reCOGnizer)']) + list(lists[f'{ide} (eggNOG mapper)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (UPIMAPI)'])
len(set(list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (eggNOG mapper)']) + list(lists[f'{ide} (mantis)'])) - lists[f'{ide} (reCOGnizer)'])
# 4
len(set(list(lists[f'{ide} (UPIMAPI)']) + list(lists[f'{ide} (eggNOG mapper)']) + list(lists[f'{ide} (mantis)']) + list(lists[f'{ide} (reCOGnizer)']))) -
```

## Publication

This paper is still not published.