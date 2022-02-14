from multiprocessing import Pool, Manager
from subprocess import check_output
from tqdm import tqdm
from paper_utils import run_pipe_command, split
from random import choice, seed
import pandas as pd


# not in use, awk solution is not working
def db_selection(db_file, ids_file, out_file):
    run_pipe_command(f'bash annotation_paper/scripts/db_selection.sh {db_file} {ids_file} {out_file}')


dbs = [
    '01', '15038598', '30077195', '105270180', '165424568', '45115792', '120308777', '180463165', '60154389',
    '135347374', '195501762', '75192986', '150385971', '210540359', '90231583']


def select_uniprot(up_file, ids_file, output):
    ids = open(ids_file).readlines()
    uniprot = open(up_file)
    line = next(uniprot)
    l = int(check_output(f'grep -c ">" {up_file}', shell=True).decode('utf8'))
    pbar = tqdm(total=l)
    print(f'Length: {l}')
    with open(output, 'w') as f:
        while line is not None:
            if line in ids:
                line = next(uniprot)
                while line is not None and not line.startswith('>'):
                    line = next(uniprot, None)
            else:
                f.write(line)
                line = next(uniprot, None)
                while line is not None and not line.startswith('>'):
                    f.write(line)
                    line = next(uniprot, None)
            pbar.update(1)
    pbar.close()


def get_upper_taxids(taxid, tax_df):
    if taxid == '0':
        return []
    taxids = []
    while taxid != '1' and taxid != 'Taxon':
        taxids.append(taxid)
        taxid = tax_df.loc[taxid]['parent_taxid']
    return taxids


def write_lineages(taxids, taxonomy_df, output):
    with open(output, 'w') as f:
        for taxid in tqdm(taxids):
            lineage = get_upper_taxids(taxid, taxonomy_df)
            f.write(f'{taxid}\t{",".join(lineage)}\n')


def get_lineages_multiprocessing(taxids, taxonomy_df, output, threads=15):
    print(f'Listing all parent tax IDs for {len(taxids)} tax IDs (this may take a while, time for coffee?)')
    taxids_groups = list(split(list(taxids), threads))
    Pool(threads).starmap(write_lineages, [(
            taxids_groups[i], taxonomy_df, f'{output}/taxonomy_{i}.tsv') for i in range(len(taxids_groups))])


def reestructure_taxonomy(tax_tsv, output):
    tax_tsv = pd.read_csv(tax_tsv, sep='\t')
    taxids = tax_tsv[tax_tsv['rank'] == 'Species']['taxid']
    tax_tsv.set_index('taxid', inplace=True)
    get_lineages_multiprocessing(taxids, tax_tsv, out)


out = 'ann_paper'

# 1st iteration
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/ids.txt', f'ann_paper/uniprot_{i}.fasta')
        for i in range(15)])


# 2nd and 3rd iterations
for i in range(2):

    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/ids.txt', f'ann_paper/uniprot_{i}.fasta')
        for i in range(15)])
