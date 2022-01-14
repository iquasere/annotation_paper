from paper_utils import count_on_file
import pandas as pd


def percentage(value, total):
    return round(value / total * 100, 2)


out = 'ann_paper'

upimapi_res = pd.read_csv(f'{out}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')
upimapi_res = upimapi_res[upimapi_res['evalue'] < 0.001]
upimapi_res = upimapi_res.groupby('qseqid')[
    ['sseqid', 'Protein names', 'Cross-reference (KEGG)', 'EC number']].first().reset_index()
n_proteins = count_on_file('>', f'{out}/genes.fasta')
annotated = len(upimapi_res)
different_ids = len(set(upimapi_res.sseqid))
uncharacterized = (upimapi_res["Protein names"] == "Uncharacterized protein").sum()
kegg_ids = upimapi_res['Cross-reference (KEGG)'].notnull().sum()
ec_numbers = upimapi_res['EC number'].notnull().sum()


table1 = pd.DataFrame([
    ['Annotated proteins', f'{annotated} ({percentage(annotated, n_proteins)} %)'],
    ['Different UniProt IDs', f'{different_ids} ({percentage(different_ids, n_proteins)} %)'],
    ['Uncharacterized proteins', f'{uncharacterized} ({percentage(uncharacterized, n_proteins)} %)'],
    ['Proteins with assigned KEGG ID', f'{kegg_ids} ({percentage(kegg_ids, n_proteins)} %)'],
    ['Proteins with assigned EC number', f'{ec_numbers} ({percentage(ec_numbers, n_proteins)} %)']]).set_index(0)
table1.columns = ['Number of genes and respective percentage']
table1.to_excel(f'{out}/table1.xlsx')

recognizer_res = pd.read_csv(f'{out}/recognizer_genomes/reCOGnizer_results.tsv', sep='\t')
recognizer_res = recognizer_res[recognizer_res['evalue'] < 0.001]
cog_res = recognizer_res[recognizer_res['DB ID'].str.startswith('COG')]
table2 = pd.merge(
    upimapi_res[upimapi_res['Protein names'] == 'Uncharacterized protein']['qseqid'], cog_res, on='qseqid', how='left')
no_cog_id = table2.groupby('qseqid')['cog'].first().reset_index()['cog'].isnull().sum()
table2 = table2.groupby(['COG general functional category', 'COG functional category'])['COG functional category'].count()
table2 = table2.append(pd.Series({('No COG ID', '-'): no_cog_id}))
table2.to_excel(f'{out}/table2.xlsx')

table3 = []
for db in ['CDD', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM', 'Pfam', 'Smart', 'COG', 'KOG']:
    df = pd.read_excel(f'{out}/recognizer_genomes/reCOGnizer_results.xlsx', sheet_name=db)
    line = {'Database': db,
            'Annotated proteins': f'{len(df)} ({percentage(len(df), n_proteins)} %)'}
    if db not in ['Smart', 'KOG']:
        line['Proteins with assigned EC number'] = \
            f'{df["EC number"].notnull().sum()} ({percentage(df["EC number"].notnull().sum(), n_proteins)} %)'
    else:
        line['Proteins with assigned EC number'] = '0 (0 %)'
    table3.append(line)
table3.to_excel(f'{out}/table3.xlsx')

# Table 4 is manually set

table5 = pd.read_excel(f'{out}/tools_comparison.xlsx')
table5.to_excel(f'{out}/table5.xlsx', index=False)
