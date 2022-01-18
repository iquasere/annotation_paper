from paper_utils import count_on_file, is_same, parse_mantis_consensus, parse_gff
import pandas as pd


def percentage(value, total):
    return round(value / total * 100, 2)


def read_res_df(filename):
    df = pd.read_csv(filename, sep='\t')
    df.fillna(value='', inplace=True)
    for col in df.columns.tolist():
        if 'EC number (' in col or 'COG (' in col:
            df[col] = df[col].str.split(';')
    return df


out = 'ann_paper'
upimapi_res = pd.read_csv(f'{out}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')
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
pd.DataFrame(table3).to_excel(f'{out}/table3.xlsx')

# Table 4 is manually set
# Table 5 is made in tools_benchmark

# Table 6 has one column with times of analysis manually set
print('Table 6 in progress')
n_prot_metagenome = count_on_file('>', f'{out}/genes_MG.fasta')
n_ann_upimapi_and_recognizer = len(set(
    list(pd.read_csv(f'{out}/upimapi_metagenome/UPIMAPI_results.tsv', sep='\t')['qseqid']) +
    list(pd.read_csv(f'{out}/recognizer_metagenome/reCOGnizer_results.tsv', sep='\t')['qseqid'])))
print('UPIMAPI + reCOGnizer done')
n_ann_mantis = len(set(parse_mantis_consensus(f'{out}/mantis_metagenome/consensus_annotation.tsv')['Query']))
print('Mantis done')
n_ann_eggnog_mapper = len(set(pd.read_csv(
    f'{out}/eggnog_mapper_metagenome/eggnog_results.emapper.annotations', sep='\t', skiprows=4, skipfooter=3)['#query']))
print('eggNOG-mapper done')
n_prot_prokka = count_on_file('>', f'{out}/prokka_metagenome/PROKKA_01132022.faa')
prokka_df = pd.read_csv(f'{out}/prokka_metagenome/PROKKA_01132022.tsv', sep='\t')
n_ann_prokka = (prokka_df['ftype'] == 'CDS').sum() - (prokka_df['product'] == 'hypothetical protein').sum()
print('Prokka done')
n_prot_dfast = count_on_file('>', f'{out}/dfast_metagenome/protein.faa')
dfast_df = parse_gff(f'{out}/dfast_metagenome/genome_trimmed.gff')
n_ann_dfast = (dfast_df['feature'] == 'CDS').sum() - (dfast_df['product'] == 'hypothetical protein').sum()
print('DFAST done')
table6 = pd.DataFrame([
    ['UPIMAPI + reCOGnizer', f'{n_ann_upimapi_and_recognizer} ({percentage(n_ann_upimapi_and_recognizer, n_prot_metagenome)} %)'],
    ['Mantis', f'{n_ann_mantis} ({percentage(n_ann_mantis, n_prot_metagenome)} %)'],
    ['eggNOG-mapper', f'{n_ann_eggnog_mapper} ({percentage(n_ann_eggnog_mapper, n_prot_metagenome)} %)'],
    ['Prokka', f'{n_ann_prokka} ({percentage(n_ann_prokka, n_prot_prokka)} %)'],
    ['DFAST', f'{n_ann_dfast} ({percentage(n_ann_dfast, n_prot_dfast)} %)']])
table6.columns = ['Tool', '# of protein annotated']
table6.to_excel(f'{out}/table_6.xlsx', index=False)

# Tables S1, S2, S3 and S4 were manually composed

# Table S5
upimapi_res.iloc[:50].to_csv(f'{out}/Table S5.tsv', sep='\t', index=False)

# Tables S6 and S7 are the e-value benchmarks, composed in the respective script

# Table S8
pd.read_excel(f'{out}/recognizer_genomes/reCOGnizer_results.xlsx', sheet_name='COG').iloc[:50].to_csv(
    f'{out}/Table S8.tsv', sep='\t', index=False)

res_df = read_res_df(f'{out}/res_df.tsv')
tools = ['UPIMAPI + reCOGnizer', 'eggNOG mapper', 'mantis', 'Prokka', 'DFAST']
table_s_9 = []
for fide in ['COG', 'EC number']:
    print(fide)
    for tool in tools:
        print(tool)
        other_tools = [ftool for ftool in tools if ftool != tool]
        df = res_df[res_df[f'{fide} ({tool})'].apply(lambda x: x != [''])].reset_index()
        unique_matches = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 0]); print(unique_matches)
        matches_with_1 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 1]); print(matches_with_1)
        matches_with_2 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 2]); print(matches_with_2)
        matches_with_3 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 3]); print(matches_with_3)
        matches_with_all = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 4]); print(matches_with_all)
        table_s_9.append(
            {'Idenfitier': fide, 'Tool': tool, 'Unique matches': unique_matches, 'Matches with 1 tool': matches_with_1,
             'Matches with 2 tools': matches_with_2, 'Matches with 3 tools': matches_with_3,
             'Matches with all': matches_with_all, '# of assignments': len(df)})
pd.DataFrame(table_s_9).to_excel(f'{out}/table_s_10.xlsx', index=False)
