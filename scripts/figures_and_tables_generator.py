from paper_utils import count_on_file, is_same, parse_mantis_consensus, parse_gff, parse_blast, get_fasta_keys
import pandas as pd
import numpy as np
import re


def percentage(value, total):
    return f'{value} ({round(value / total * 100, 2)} %)'


def read_res_df(filename):
    df = pd.read_csv(filename, sep='\t')
    df.fillna(value='', inplace=True)
    for col in df.columns.tolist():
        if 'EC number (' in col or 'COG (' in col:
            df[col] = df[col].str.split(';')
    return df


out = 'ann_paper'
print('Table 1 in progress')
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
    ['Annotated proteins', percentage(annotated, n_proteins)],
    ['Different UniProt IDs', percentage(different_ids, n_proteins)],
    ['Uncharacterized proteins', percentage(uncharacterized, n_proteins)],
    ['Proteins with assigned KEGG ID', percentage(kegg_ids, n_proteins)],
    ['Proteins with assigned EC number', percentage(ec_numbers, n_proteins)]]).set_index(0)
table1.columns = ['Number of genes and respective percentage']
table1.to_excel(f'{out}/table1.xlsx')

print('Table 2 in progress')
recognizer_res = pd.read_csv(f'{out}/recognizer_genomes/reCOGnizer_results.tsv', sep='\t')
recognizer_res = recognizer_res[recognizer_res['evalue'] < 0.001]
cog_res = recognizer_res[recognizer_res['DB ID'].str.startswith('COG')]
table2 = pd.merge(
    upimapi_res[upimapi_res['Protein names'] == 'Uncharacterized protein']['qseqid'], cog_res, on='qseqid', how='left')
no_cog_id = table2.groupby('qseqid')['cog'].first().reset_index()['cog'].isnull().sum()
table2 = table2.groupby(['COG general functional category', 'COG functional category'])['COG functional category'].count()
table2 = table2.append(pd.Series({('No COG ID', '-'): no_cog_id}))
table2.to_excel(f'{out}/table2.xlsx')

print('Table 3 in progress')
table3 = []
for db in ['CDD', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM', 'Pfam', 'Smart', 'COG', 'KOG']:
    print(db)
    df = pd.read_excel(f'{out}/recognizer_genomes/reCOGnizer_results.xlsx', sheet_name=db)
    cols = df.columns.tolist()
    cols.remove('qseqid')
    df = df.groupby('qseqid')[cols].first().reset_index()
    line = {'Database': db,
            'Annotated proteins': percentage(len(df), n_proteins)}
    if db not in ['Smart', 'KOG']:
        line['Proteins with assigned EC number'] = percentage(df["EC number"].notnull().sum(), n_proteins)
    else:
        line['Proteins with assigned EC number'] = '0 (0 %)'
    table3.append(line)
pd.DataFrame(table3).to_excel(f'{out}/table3.xlsx')


# Table 4
print('Table 4 in progress')
iterations = ['', '/first_group', '/second_group', '/third_group', '/fourth_group', '/fifth_group']
prokka_dates = ["01032022", "03052022", "02192022", "02252022", "03072022", "03112022"]
table5 = pd.ExcelWriter(f'{out}/table5_all.xlsx', engine='xlsxwriter')
df = pd.read_excel(f'{out}/table5.xlsx')
all_t5 = pd.DataFrame(index=df.index, columns=df.columns)
all_t5.fillna('', inplace=True)
for col in all_t5.columns:
    if col in ['Qualifier', 'Tool']:
        all_t5[col] = df[col]
    else:
        all_t5[col] = all_t5[col].apply(lambda x: x.split())
for i in range(len(iterations)):
    df = pd.read_excel(f'{out}{iterations[i]}/table5.xlsx')
    df.to_excel(table5, sheet_name=str(i), index=False)
    n_proteins = count_on_file('>', f'{out}{iterations[i]}/proteomes.fasta')
    n_proteins_dfast = count_on_file('>', f'{out}{iterations[i]}/dfast_genomes/protein.faa')
    n_proteins_prokka = count_on_file('>', f'{out}{iterations[i]}/prokka_genomes/PROKKA_{prokka_dates[i]}.faa')
    if iterations[i] != '/second_group':
        for j in df.index:
            for col in df.columns:
                if col not in ['Qualifier', 'Tool']:
                    if col in ['TPs', 'FPs', 'FNs']:
                        all_t5.at[j, col] += [df.iloc[j][col] / (df.iloc[j]['TPs'] + df.iloc[j]['FPs'] + df.iloc[j]['FNs'])]
                    elif col in ['# of total identifications']:
                        n_p = n_proteins_dfast if df.iloc[j]['Tool'] == 'DFAST' else (n_proteins_prokka if df.iloc[j]['Tool'] == 'Prokka' else n_proteins)
                        all_t5.at[j, col] += [df.iloc[j][col] / n_p]
                    else:
                        all_t5.at[j, col] += [df.iloc[j][col]]
table5.save()
final_t5 = pd.DataFrame(index=all_t5.index, columns=all_t5.columns)
for col in all_t5.columns:
    if col in ['Qualifier', 'Tool']:
        final_t5[col] = all_t5[col]
for j in df.index:
    for col in all_t5.columns:
        if col not in ['Qualifier', 'Tool']:
            final_t5.at[j, col] = f'{round(np.mean(all_t5.iloc[j][col]) * 100, 2)} +- {round(np.std(all_t5.iloc[j][col]) * 100, 2)} %'
final_t5.to_excel(f'{out}/final_t5.xlsx', index=False)

# Table 5
print('Table 5 in progress')
n_prot_real = count_on_file('>', f'{out}/genes_MG.fasta')
upimapi_real = pd.read_csv(f'{out}/upimapi_metagenome/UPIMAPI_results.tsv', sep='\t')
upimapi_real = upimapi_real.groupby('qseqid')[['Entry', 'EC number', 'Cross-reference (eggNOG)']].first().reset_index()
recognizer_real = pd.read_csv(f'{out}/recognizer_metagenome/reCOGnizer_results.tsv', sep='\t')
cog_ser = recognizer_real[recognizer_real['cog'].notnull()].groupby('qseqid')['cog'].apply(
    lambda x: '; '.join(set(x)))
ecs_ser = recognizer_real[recognizer_real['EC number'].notnull()].groupby('qseqid')['EC number'].apply(
    lambda x: '; '.join(set(x)))
recognizer_real = recognizer_real.groupby('qseqid')['DB ID'].first().reset_index()
mantis_real = parse_mantis_consensus(f'{out}/mantis_metagenome/consensus_annotation.tsv')
eggnog_real = pd.read_csv(
    f'{out}/eggnog_mapper_metagenome/eggnog_results.emapper.annotations', sep='\t', skiprows=4, skipfooter=3)
eggnog_real.replace('-', np.nan, inplace=True)
eggnog_real['eggNOG_OGs'] = eggnog_real['eggNOG_OGs'].apply(lambda x: x.split('@')[0])
prokka_real = pd.read_csv(f'{out}/prokka_metagenome/PROKKA_01132022.tsv', sep='\t')
prokka_real = prokka_real[(prokka_real['ftype'] == 'CDS') & (prokka_real['product'] != 'hypothetical protein')]
prokka_real.rename(columns={'EC_number': 'EC number (Prokka)', 'COG': 'COG (Prokka)', 'locus_tag': 'ID (Prokka)'}, inplace=True)
prot_prokka_real = count_on_file('>', f'{out}/prokka_metagenome/PROKKA_01132022.faa')
dfast_real = parse_gff(f'{out}/dfast_metagenome/genome_trimmed.gff')
dfast_real = dfast_real[(dfast_real['feature'] == 'CDS')  & (dfast_real['product'] != 'hypothetical protein')]
dfast_real['COG'] = dfast_real['note'].apply(
    lambda x: np.nan if 'COG:' not in x else re.split(':| ', x.split('COG:')[-1])[0])
dfast_real.rename(columns={'ID': 'ID (DFAST)', 'EC_number': 'EC number (DFAST)', 'COG': 'COG (DFAST)'}, inplace=True)
prot_dfast_real = count_on_file('>', f'{out}/dfast_metagenome/protein.faa')
real_df = pd.DataFrame(get_fasta_keys(f'{out}/genes_MG.fasta'), columns=['qseqid']).iloc[:-1]   # get one empty line
real_df['qseqid'] = real_df['qseqid'].apply(lambda x: x.split()[0])
real_df = pd.merge(
    real_df, upimapi_real[['qseqid', 'Entry', 'EC number', 'Cross-reference (eggNOG)']], on='qseqid', how='left')
real_df = pd.merge(real_df, recognizer_real, on='qseqid', how='left')
real_df = pd.merge(real_df, cog_ser, left_on='qseqid', right_index=True, how='left')
real_df = pd.merge(real_df, ecs_ser, left_on='qseqid', right_index=True, how='left')
real_df = pd.merge(real_df, eggnog_real[['#query', 'EC', 'eggNOG_OGs']], left_on='qseqid', right_on='#query', how='left')
real_df = pd.merge(real_df, mantis_real[['Query', 'enzyme_ec', 'cog']], left_on='qseqid', right_on='Query', how='left')
real_df.rename(columns={
    'Entry': 'ID (UPIMAPI)',
    'DB ID': 'ID (reCOGnizer)',
    '#query': 'ID (eggNOG mapper)',
    'Query': 'ID (Mantis)',
    'EC number_x': 'EC number (UPIMAPI)',
    'EC number_y': 'EC number (reCOGnizer)',
    'EC': 'EC number (eggNOG mapper)',
    'enzyme_ec': 'EC number (Mantis)',
    'Cross-reference (eggNOG)': 'COG (UPIMAPI)',
    'cog_x': 'COG (reCOGnizer)',
    'eggNOG_OGs': 'COG (eggNOG mapper)',
    'cog_y': 'COG (Mantis)'}, inplace=True)
real_df['ID (UPIMAPI + reCOGnizer)'] = real_df['ID (UPIMAPI)'].combine_first(real_df['ID (reCOGnizer)'])
real_df['EC number (UPIMAPI + reCOGnizer)'] = real_df['EC number (UPIMAPI)'].combine_first(real_df['EC number (reCOGnizer)'])
real_df['COG (UPIMAPI + reCOGnizer)'] = real_df['COG (UPIMAPI)'].combine_first(real_df['COG (reCOGnizer)'])

table6 = []
for tool in ['UPIMAPI + reCOGnizer', 'UPIMAPI', 'reCOGnizer', 'Mantis', 'eggNOG mapper']:
    table6.append([tool,
                   percentage(real_df[f'ID ({tool})'].notnull().sum(), n_prot_real),
                   percentage(real_df[f'EC number ({tool})'].notnull().sum(), n_prot_real),
                   percentage(real_df[f'COG ({tool})'].notnull().sum(), n_prot_real)])
table6.append(['Prokka',
               percentage(prokka_real[f'ID (Prokka)'].notnull().sum(), prot_prokka_real),
               percentage(prokka_real[f'EC number (Prokka)'].notnull().sum(), prot_prokka_real),
               percentage(prokka_real[f'COG (Prokka)'].notnull().sum(), prot_prokka_real)])
table6.append(['DFAST',
               percentage(dfast_real[f'ID (DFAST)'].notnull().sum(), prot_dfast_real),
               percentage(dfast_real[f'EC number (DFAST)'].notnull().sum(), prot_dfast_real),
               percentage(dfast_real[f'COG (DFAST)'].notnull().sum(), prot_dfast_real)])
pd.DataFrame(table6, columns=['Tool', '# of proteins annotated', '# of proteins with EC number assigned',
                              '# of proteins with COG assigned']).to_excel(f'{out}/table_6.xlsx', index=False)

# Tables S1, S2, S3 and S4 were manually composed

# Table S5
upimapi_res.iloc[:50].to_csv(f'{out}/Table S5.tsv', sep='\t', index=False)

# Tables S6 and S8 are the e-value benchmarks, composed in the respective script

# Table S7
print('Table S7 in progress')
pd.read_excel(f'{out}/recognizer_genomes/reCOGnizer_results.xlsx', sheet_name='COG').iloc[:50].to_csv(
    f'{out}/Table S7.tsv', sep='\t', index=False)

# Table S9
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
            ) == 0])
        matches_with_1 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 1])
        matches_with_2 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 2])
        matches_with_3 = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 3])
        matches_with_all = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 4])
        table_s_9.append(
            {'Idenfitier': fide, 'Tool': tool, 'Unique matches': unique_matches, 'Matches with 1 tool': matches_with_1,
             'Matches with 2 tools': matches_with_2, 'Matches with 3 tools': matches_with_3,
             'Matches with all': matches_with_all, '# of assignments': len(df)})
pd.DataFrame(table_s_9).to_excel(f'{out}/table_s_10.xlsx', index=False)


# Figure S2
tools = ['UPIMAPI + reCOGnizer', 'eggNOG mapper', 'mantis']
for fide in ['COG', 'EC number']:
    print(fide)
    for tool in tools:
        other_tools = [ftool for ftool in tools if ftool != tool]
        print(tool)
        df = res_df[res_df[f'{fide} ({tool})'].apply(lambda x: x != [''])].reset_index()
        unique_matches = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 0])
        print(f'unique: {unique_matches}')
        for other_tool in other_tools:
            other_other_tool = [otool for otool in other_tools if otool != other_tool]
            matches_with_tool = len([
                df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                    df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i])]
                ) == 1 and sum([is_same(
                    df[f'{fide} ({tool})'][i], df[f'{fide} ({other_other_tool[0]})'][i])]) == 0])
            print(f'{tool} with {other_tool}: {matches_with_tool}')
        matches_with_all = len([
            df[f'{fide} ({tool})'][i] for i in range(len(df)) if sum([is_same(
                df[f'{fide} ({tool})'][i], df[f'{fide} ({other_tool})'][i]) for other_tool in other_tools]
            ) == 2])
        print(f'matches with all: {matches_with_all}')

# some final values for text
# % of proteins annotated
n_proteins = [count_on_file('>', f'{out}{iterations[i]}/proteomes.fasta') for i in range(len(iterations))]
n_proteins_prokka = [count_on_file('>', f'{out}{iterations[i]}/prokka_genomes/PROKKA_{prokka_dates[i]}.faa') for i in range(len(iterations))]
n_proteins_dfast = [count_on_file('>', f'{out}{iterations[i]}/dfast_genomes/protein.faa') for i in range(len(iterations))]
n_ids_u = [len(set(pd.read_csv(f'{out}{iterations[i]}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')['qseqid'])) for i in range(len(iterations))]
n_ids_r = [len(set(pd.read_csv(f'{out}{iterations[i]}/recognizer_genomes/reCOGnizer_results.tsv', sep='\t')['qseqid'])) for i in range(len(iterations))]
n_ids_e = [len(set(pd.read_csv(f'{out}{iterations[i]}/eggnog_mapper_genomes/eggnog_results.emapper.annotations', sep='\t', skiprows=4, skipfooter=3)['#query'])) for i in range(len(iterations))]
n_ids_m = [len(set(parse_mantis_consensus(f'{out}{iterations[i]}/mantis_genomes/consensus_annotation.tsv')['Query'])) for i in range(len(iterations))]
n_ids_p = []
for i in range(len(iterations)):
    df = pd.read_csv(f'{out}{iterations[i]}/prokka_genomes/PROKKA_{prokka_dates[i]}.tsv', sep='\t')
    df = df[(df['ftype'] == 'CDS') & (df['product'] == 'hypothetical protein')]
    n_ids_p.append(len(set(df['locus_tag'])))
n_ids_d = []
for i in range(len(iterations)):
    df = parse_gff(f'{out}{iterations[i]}/dfast_genomes/genome_trimmed.gff')
    df = df[(df['feature'] == 'CDS') & (df['product'] == 'hypothetical protein')]
    n_ids_d.append(len(set(df['locus_tag'])))
for t in [(n_ids_u, n_proteins), (n_ids_r, n_proteins), (n_ids_e, n_proteins), (n_ids_m, n_proteins), (n_ids_p, n_proteins_prokka), (n_ids_d, n_proteins_dfast)]:
    print(np.mean([i/j for i, j in zip(t[0], t[1])]))

