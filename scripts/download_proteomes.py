import requests

taxids = ['2203', '2223', '2162', '119484', '35554', '29543', '863']
out = 'ann_paper'
with open(f'{out}/proteomes.fasta', 'w') as f:
    for taxid in taxids:
        print(f'Downloading proteome for taxid: {taxid}')
        res = requests.get(f'https://www.uniprot.org/uniprot/?query=taxonomy:{taxid}&format=fasta')
        f.write(res.content.decode('utf8'))
