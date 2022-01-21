import pandas as pd

out = 'ann_paper'
readcounts = pd.read_csv(f'{out}/readcounts.tsv', sep='\t')
relation = pd.read_csv(f'{out}/simulated/uniprotinfo.tsv', sep='\t')
relation.columns = relation.columns.tolist()[:-1] + ['name']
readcounts = pd.merge(readcounts, relation[['name', 'Entry']], on='name')
upimapi_res = pd.read_csv(f'{out}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')
cols = upimapi_res.columns.tolist()
cols.remove('qseqid')
upimapi_res = upimapi_res.groupby('qseqid')[cols].first().reset_index()
keggcharter_input = pd.merge(upimapi_res, readcounts, on='Entry', how='left')
keggcharter_input.to_csv(f'{out}/keggcharter_input.tsv', sep='\t', index=False)
