from multiprocessing import Pool
from subprocess import check_output

from tqdm import tqdm
from paper_utils import run_pipe_command

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


with Pool(processes=15) as p:
    p.starmap(select_uniprot, [
        (f'resources_directory/split_uniprot.{dbs[i]}', 'ann_paper/ids.txt', f'ann_paper/uniprot_{i}.fasta')
        for i in range(15)])
