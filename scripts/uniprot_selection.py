from multiprocessing import Pool
from paper_utils import run_pipe_command

def db_selection(db_file, ids_file, out_basename):
    run_pipe_command(f'bash annotation_paper/scripts/db_selection.sh {db_file} {ids_file} {out_basename}.fasta')


dbs = [
    '01', '15038598', '30077195', '105270180', '165424568', '45115792', '120308777', '180463165', '60154389',
    '135347374', '195501762', '75192986', '150385971', '210540359', '90231583']

with Pool(processes=1) as p:
    p.starmap(db_selection, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/ids.qtxt', f'ann_paper/uniprot_{i}.fasta')
        for i in range(15)])

