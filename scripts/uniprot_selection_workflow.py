import pathlib
from multiprocessing import Pool
from subprocess import check_output
from tqdm import tqdm
from paper_utils import run_pipe_command, run_command
from random import choices, seed
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


def download_genome(entry, out_dir):
    link = f'http://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/{"_".join(entry["core_db"].split("_")[:3])}/' \
           f'{entry["species"]}/dna/{entry["species"].capitalize()}.{entry["assembly"].replace(" ", "_").replace("#", "_")}' \
           f'.dna.toplevel.fa.gz'
    run_command(f'wget {link} -O {out_dir}/{"_".join(entry["name"].split()[:-1])}.fna.gz')


out = 'ann_paper'

# Iteration zero
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/ids.txt', f'ann_paper/uniprot_{i}.fasta')
        for i in range(15)])

run_command(f'wget http://ftp.ensemblgenomes.org/pub/bacteria/current/species_EnsemblBacteria.txt -P {out}')
prok_df = pd.read_csv(
    'species_EnsemblBacteria.txt', sep='\t', encoding='latin-1', skiprows=1, index_col=False,
    names=['name', 'species', 'division', 'taxonomy_id', 'assembly', 'assembly_accession', 'genebuild', 'variation',
           'microarray', 'pan_compara', 'peptide_compara', 'genome_alignments', 'other_alignments', 'core_db',
           'species_id'])

# 1st iteration - Comments of using this instead of that report the retrieval of proteomes from UniProt
seed(0)
first_species = choices(prok_df['name'].tolist(), k=8)
first_species.remove('Halomonas sp. Choline-3u-9 (GCA_002836495)')
# Halomonas sp. Choline-3u-9 is not present in UniProt taxonomy, nor anything close, so was rejected
# Enterococcus faecalis EnGen0354 was selected from the "Enterococcus faecalis"
pathlib.Path(f'{out}/first_group').mkdir(exist_ok=True, parents=True)
for species in first_species:
    download_genome(prok_df[prok_df['name'] == species].reset_index().iloc[0], f'{out}/first_group')
run_pipe_command("awk '{print $1}' ann_paper/first_group/*.faa > ann_paper/first_group/proteomes.fasta")
run_pipe_command('grep ">" ann_paper/first_group/proteomes.fasta > ann_paper/first_group/ids.txt')
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/first_group/ids.txt',
         f'ann_paper/first_group/uniprot_{i}.fasta') for i in range(15)])

# 2nd iteration
seed(1)
second_species = choices(prok_df['name'].tolist(), k=7)
# Serratia fonticola AU-P3(3) was used instead of Serratia fonticola str. 5l
# KLEP7 was selected from the "Klebsiella pneumoniae"
# HISS2 was selected from the "Histophilus somni 2336"
# Neisseria gonorrhoeae DGI2 was selected from the Neisseria gonorrhoeae
pathlib.Path(f'{out}/second_group').mkdir(exist_ok=True, parents=True)
for species in second_species:
    download_genome(prok_df[prok_df['name'] == species].reset_index().iloc[0], f'{out}/second_group')
run_pipe_command("awk '{print $1}' ann_paper/second_group/*.faa > ann_paper/second_group/proteomes.fasta")
run_pipe_command('grep ">" ann_paper/second_group/proteomes.fasta > ann_paper/second_group/ids.txt')
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/second_group/ids.txt',
         f'ann_paper/second_group/uniprot_{i}.fasta') for i in range(15)])
# No OGs were present in this species

# 3rd iteration
seed(2)
third_species = choices(prok_df['name'].tolist(), k=7)
# Thermococcus celericrescens was used instead of Thermococcus celericrescens str. DSM 17994
# Tateyamaria omphalii was used instead of Tateyamaria omphalii str. DOK1-4
# Salmonella enterica subsp. indica serovar 6,14,25:z10:1,(2),7 str. 1121 was used instead of Salmonella enterica subsp. enterica serovar 4,[5],12:i:-
# Novosphingobium sp. PP1Y was used instead of Novosphingobium sp. Rr 2-17 str. RR2-17
pathlib.Path(f'{out}/third_group').mkdir(exist_ok=True, parents=True)
for species in third_species:
    download_genome(prok_df[prok_df['name'] == species].reset_index().iloc[0], f'{out}/third_group')
run_pipe_command("awk '{print $1}' ann_paper/third_group/*.faa > ann_paper/third_group/proteomes.fasta")
run_pipe_command('grep ">" ann_paper/third_group/proteomes.fasta > ann_paper/third_group/ids.txt')
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/third_group/ids.txt',
         f'ann_paper/third_group/uniprot_{i}.fasta') for i in range(15)])

# 4th iteration
seed(3)
fourth_species = choices(prok_df['name'].tolist(), k=7)
# Clostridiales bacterium CHKCI001 was used as the Clostridiales bacterium
# Legionella oakridgensis ATCC 33761 = DSM 21215 used instead of Legionella oakridgensis ATCC 33761 = DSM 21215 str. OR-10, ATCC 33761
# FERIS used as the Fervidobacterium islandicum str. AW-1
# Muribaculaceae bacterium Isolate-036 (Harlan) used as the Muribaculaceae bacterium
# Acidobacteria bacterium AB60 used as the Acidobacteria bacterium
pathlib.Path(f'{out}/fourth_group').mkdir(exist_ok=True, parents=True)
for species in fourth_species:
    download_genome(prok_df[prok_df['name'] == species].reset_index().iloc[0], f'{out}/fourth_group')
run_pipe_command("awk '{print $1}' ann_paper/fourth_group/*.faa > ann_paper/fourth_group/proteomes.fasta")
run_pipe_command('grep ">" ann_paper/fourth_group/proteomes.fasta > ann_paper/fourth_group/ids.txt')
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/fourth_group/ids.txt',
         f'ann_paper/fourth_group/uniprot_{i}.fasta') for i in range(15)])

seed(4)
fifth_species = choices(prok_df['name'].tolist(), k=8)
# Clavibacter michiganensis subsp. michiganensis (strain NCPPB 382) used instead of Clavibacter michiganensis str. CFBP8017
# Bacillus pseudomycoides DSM 12442 used as the Bacillus pseudomycoides
# Bordetella genomosp. 6 only contains tiny proteomes in UniProt (in NCBI has 4758, in UniProt only 10 %)
# GEOS8 used as the Geobacter sp.
# STRCO used as the Streptomyces albidoflavus
# Rhodobacteraceae bacterium HLUCCA08 used instead of Rhodobacteraceae bacterium SM1903
pathlib.Path(f'{out}/fifth_group').mkdir(exist_ok=True, parents=True)
for species in fifth_species:
    download_genome(prok_df[prok_df['name'] == species].reset_index().iloc[0], f'{out}/fifth_group')
run_pipe_command("awk '{print $1}' ann_paper/fifth_group/*.faa > ann_paper/fifth_group/proteomes.fasta")
run_pipe_command('grep ">" ann_paper/fifth_group/proteomes.fasta > ann_paper/fifth_group/ids.txt')
with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/fifth_group/ids.txt',
         f'ann_paper/fifth_group/uniprot_{i}.fasta') for i in range(15)])
