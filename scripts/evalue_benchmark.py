from paper_utils import calculate_quality_metrics, count_on_file
import pandas as pd

def is_true_positive(l1, l2, threshold=0.5):
    l1, l2 = set(l1), set(l2)       # must be set because some proteins will have multiple times the same identification attributed because of multiple positions of domain
    return sum([i1 in l2 for i1 in l1]) / len(l1) > threshold


out = 'ann_paper'
evalues = [1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-18, 1e-21, 1e-24, 1e-27, 1e-30]

# UPIMAPI
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

# reCOGnizer
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

