import pandas as pd

out = 'ann_paper'
readcounts = pd.read_csv(f'{out}/readcounts.tsv', sep='\t')
upimapi_res = pd.read_csv(f'{out}/upimapi_genomes/UPIMAPI_results.tsv', sep='\t')
cols = upimapi_res.columns.tolist()
cols.remove('qseqid')
upimapi_res = upimapi_res.groupby('qseqid')[cols].first().reset_index()
keggcharter_input = pd.merge(upimapi_res, readcounts, left_on='Entry', right_on='name', how='left')
keggcharter_input.to_csv(f'{out}/keggcharter_input.tsv', sep='\t', index=False)
