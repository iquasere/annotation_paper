from paper_utils import parse_mantis_consensus, calculate_quality_metrics, parse_blast, blast_consensus, parse_gff, is_same
import pandas as pd
import numpy as np
from collections import Counter
import re


def true_positive_funcs(l1, l2, threshold=0.5):
    return sum([i1 in l2 for i1 in l1]) / len(l1) > threshold


def add_to_quality_benchmark(quality_benchmark_list, results_df, ide, tool):
    identifications = (~results_df[f'{ide} ({tool})'].isin([['']])).sum()
    df = results_df[~results_df[f'{ide} (UniProt)'].isin([['']])]
    print(f'{tool} results for determining {ide}')
    tps = sum([is_same(df.iloc[i][f'{ide} (UniProt)'], df.iloc[i][f'{ide} ({tool})']) for i in range(len(df))])
    fps = (~df[f'{ide} ({tool})'].isin([['']])).sum() - tps
    fns = (df[f'{ide} ({tool})'].isin([['']])).sum()
    print(f'TPs:{tps} FPs:{fps} FNs:{fns}')
    precision, recall, f1_score = calculate_quality_metrics(tps, fps, fns)
    quality_benchmark_list.append(
        {'Qualifier': ide, 'Tool': tool, 'TPs': tps, 'FPs': fps, 'FNs': fns, 'precision': precision, 'recall': recall,
         'f1_score': f1_score, '# of total identifications': identifications})
    return quality_benchmark_list


def counts_funcs_series(fdf, fcol):
    def count_funcs(series):
        funcs_list = []
        for funcs in series:
            funcs_list += funcs
        return Counter(funcs_list)
    return pd.Series(count_funcs(fdf[fcol]), name=fcol)


def df_post_processing(df):
    df.fillna(value='', inplace=True)
    for col in df.columns.tolist():
        if 'EC number (' in col or 'COG (' in col:
            df[col] = df[col].str.replace(',', ';').str.replace(' ', '').str.rstrip(';').str.split(';')
    return df


def write_res_df(df, filename):
    for col in df.columns.tolist():
        if 'EC number (' in col or 'COG (' in col:
            df[col] = df[col].apply(lambda x: ';'.join(x))
    df.to_csv(filename, sep='\t', index=False)


def results_analysis(out, prokka_date):
    #readcounts = pd.read_csv(f'{out}/readcounts.tsv', sep='\t')
    upimapi_genomes = pd.read_csv(f'{out}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')
    upimapi_genomes = upimapi_genomes.groupby('qseqid')[['EC number', 'Cross-reference (eggNOG)']].first().reset_index()
    recognizer_genomes = pd.read_csv(f'{out}/recognizer_genomes/reCOGnizer_results.tsv', sep='\t')
    cog_ser = recognizer_genomes[recognizer_genomes['cog'].notnull()].groupby('qseqid')['cog'].apply(
        lambda x: '; '.join(set(x)))
    ecs_ser = recognizer_genomes[recognizer_genomes['EC number'].notnull()].groupby('qseqid')['EC number'].apply(
        lambda x: '; '.join(set(x)))
    mantis_genomes = parse_mantis_consensus(f'{out}/mantis_genomes/consensus_annotation.tsv')
    mantis_genomes['og'] = mantis_genomes['cog'].combine_first(mantis_genomes['arcog'])
    eggnog_genomes = pd.read_csv(
        f'{out}/eggnog_mapper_genomes/eggnog_results.emapper.annotations', sep='\t', skiprows=4, skipfooter=3)
    eggnog_genomes.replace('-', np.nan, inplace=True)
    eggnog_genomes['eggNOG_OGs'] = eggnog_genomes['eggNOG_OGs'].apply(lambda x: x.split('@')[0])
    uniprotinfo = pd.read_csv(f'{out}/uniprotinfo.tsv', sep='\t')
    prokka_df = pd.merge(blast_consensus(f'{out}/prokka_genomes/PROKKA_{prokka_date}.blast'), uniprotinfo, on='Entry', how='left')   # must run: diamond blastp -d ann_paper/proteomes.dmnd -q ann_paper/prokka_genomes/PROKKA_01032022.faa -o ann_paper/prokka_genomes/PROKKA_01032022.blast --very-sensitive -p 15
    prokka_tsv = pd.read_csv(f'{out}/prokka_genomes/PROKKA_{prokka_date}.tsv', sep='\t')
    prokka_df = pd.merge(prokka_df, prokka_tsv[prokka_tsv['ftype'] == 'CDS'][['locus_tag', 'EC_number', 'COG']],
                         left_on='qseqid', right_on='locus_tag', how='left')
    prokka_df.rename(
        columns={'EC number': 'EC number (UniProt)', 'EC_number': 'EC number (Prokka)',
                 'Cross-reference (eggNOG)': 'COG (UniProt)', 'COG': 'COG (Prokka)'}, inplace=True)
    dfast_df = pd.merge(blast_consensus(f'{out}/dfast_genomes/protein.blast'), uniprotinfo, on='Entry', how='left')
    gff = parse_gff(f'{out}/dfast_genomes/genome_trimmed.gff')
    gff['COG'] = gff['note'].apply(lambda x: np.nan if 'COG:' not in x else re.split(':| ', x.split('COG:')[-1])[0])
    dfast_df['qseqid'] = dfast_df['qseqid'].apply(lambda x: x.split('|')[0])
    dfast_df = pd.merge(dfast_df, gff, left_on='qseqid', right_on='ID', how='left')
    dfast_df.rename(
        columns={'EC number': 'EC number (UniProt)', 'EC_number': 'EC number (DFAST)',
                 'Cross-reference (eggNOG)': 'COG (UniProt)', 'COG': 'COG (DFAST)'}, inplace=True)
    # Build final analysis DF
    res_df = blast_consensus(f'{out}/genes.blast')
    #res_df = res_df.groupby('qseqid')[res_df.columns.tolist()[1:]].first().reset_index()
    #res_df['Entry'] = res_df['sseqid'].apply(lambda x: x.split('|')[1])
    res_df = pd.merge(res_df, uniprotinfo, on='Entry', how='left')
    res_df = pd.merge(res_df, upimapi_genomes[['qseqid', 'EC number', 'Cross-reference (eggNOG)']], on='qseqid', how='left')
    res_df = pd.merge(res_df, cog_ser, left_on='qseqid', right_index=True, how='left')
    res_df = pd.merge(res_df, ecs_ser, left_on='qseqid', right_index=True, how='left')
    res_df = pd.merge(res_df, eggnog_genomes[['#query', 'EC', 'eggNOG_OGs']], left_on='qseqid', right_on='#query', how='left')
    res_df = pd.merge(res_df, mantis_genomes[['Query', 'enzyme_ec', 'og']], left_on='qseqid', right_on='Query', how='left')
    res_df.rename(columns={
        'EC number_x': 'EC number (UniProt)',
        'EC number_y': 'EC number (UPIMAPI)',
        'EC number': 'EC number (reCOGnizer)',
        'EC': 'EC number (eggNOG mapper)',
        'enzyme_ec': 'EC number (mantis)',
        'Cross-reference (eggNOG)_x': 'COG (UniProt)',
        'Cross-reference (eggNOG)_y': 'COG (UPIMAPI)',
        'cog': 'COG (reCOGnizer)',
        'eggNOG_OGs': 'COG (eggNOG mapper)',
        'og': 'COG (mantis)'}, inplace=True)
    for col in ['#query', 'Query']:
        del res_df[col]
    res_df['EC number (UPIMAPI + reCOGnizer)'] = res_df['EC number (UPIMAPI)'].combine_first(res_df['EC number (reCOGnizer)'])
    res_df['COG (UPIMAPI + reCOGnizer)'] = res_df['COG (UPIMAPI)'].combine_first(res_df['COG (reCOGnizer)'])
    res_df = pd.merge(res_df, prokka_df[['Entry', 'EC number (Prokka)', 'COG (Prokka)']], on='Entry', how='left')
    res_df = pd.merge(res_df, dfast_df[['Entry', 'EC number (DFAST)', 'COG (DFAST)']], on='Entry', how='left')
    res_df = df_post_processing(res_df)
    quality_benchmark = []
    for fide in ['EC number', 'COG']:
        for ftool in ['eggNOG mapper', 'mantis', 'Prokka', 'DFAST', 'UPIMAPI + reCOGnizer', 'UPIMAPI', 'reCOGnizer']:
            quality_benchmark = add_to_quality_benchmark(quality_benchmark, res_df, fide, ftool)
    pd.DataFrame(quality_benchmark).to_excel(f'{out}/table5.xlsx', index=False)
    write_res_df(res_df, f'{out}/res_df.tsv')
    return
    # not in the paper, not very informative
    for fide in ['EC number', 'COG']:
        heatmap_df = pd.DataFrame()
        for tool in ['UPIMAPI + reCOGnizer', 'eggNOG mapper', 'mantis', 'Prokka', 'DFAST']:
            heatmap_df = pd.merge(
                heatmap_df, counts_funcs_series(res_df, f'{fide} ({tool})'), left_index=True, right_index=True, how='outer')
        heatmap_df.fillna(value=0.0).astype(int).to_csv(f'{out}/heatmap_df_{fide}.tsv', sep='\t')


results_analysis('ann_paper', '01032022')
results_analysis('ann_paper/first_group', '03052022')
results_analysis('ann_paper/second_group', '02192022')
results_analysis('ann_paper/third_group', '02252022')
results_analysis('ann_paper/fourth_group', '03072022')
results_analysis('ann_paper/fifth_group', '03112022')
