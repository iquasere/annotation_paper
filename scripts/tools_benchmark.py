from paper_utils import parse_mantis_consensus, calculate_quality_metrics
import pandas as pd

out = 'ann_paper'
readcounts = pd.read_csv(f'{out}/readcounts.tsv', sep='\t')

#eggNOG mapper
egg_basename = f'{out}/eggnog_mapper_proteomes/eggnog_results.emapper'
eggnog_annotations = pd.read_csv(f'{egg_basename}.annotations', sep='\t', skiprows=4, skipfooter=3)


#mantis
mantis_basename = f'{out}/mantis_proteomes'
mantis_consensus = parse_mantis_consensus(f'{mantis_basename}/consensus_annotation.tsv')
mantis_integrated = pd.read_csv(f'{mantis_basename}/integrated_annotation.tsv', sep='\t')
mantis_output = pd.read_csv(f'{mantis_basename}/output_annotation.tsv', sep='\t')


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
res_df = pd.merge(res_df, upimapi_proteome[['Entry', 'EC number', 'Cross-reference (eggNOG)']], on='Entry', how='left')
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
