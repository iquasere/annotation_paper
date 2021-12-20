import os

import pandas as pd
import glob
from paper_utils import run_command, run_pipe_command


def check_bowtie2_index(index_prefix):
    files = glob.glob(index_prefix + '*.bt2')
    if len(files) < 6:
        return False
    return True


def generate_mg_index(reference, index_prefix):
    run_command(f'bowtie2-build {reference} {index_prefix}')


def align_reads(reads, index_prefix, sam, report, log=None, threads=6):
    run_command(f'bowtie2 -x {index_prefix} -1 {reads[0]} -2 {reads[1]} -S {sam} -p {threads} 1> {report} 2> {log}')


def perform_alignment(reference, reads, basename, threads=1):
    ext = f".{reference.split('.')[-1]}"

    if not check_bowtie2_index(reference.replace(ext, '_index')):
        print('INDEX files not found. Generating new ones')
        generate_mg_index(reference, reference.replace(ext, '_index'))
    else:
        print(f"INDEX was located at {reference.replace(ext, '_index')}")

    if not os.path.isfile(f'{basename}.log'):
        align_reads(reads, reference.replace(ext, '_index'), f'{basename}.sam',
                    f'{basename}_bowtie2_report.txt', log=f'{basename}.log', threads=threads)
    else:
        print(f'{basename}.log was found!')

    run_pipe_command(
        f"""samtools view -F 260 {basename}.sam | cut -f 3 | sort | uniq -c | """
        f"""awk '{{printf("%s\\t%s\\n", $2, $1)}}'""", output=f'{basename}.readcounts')


out = 'ann_paper'
readcounts = pd.DataFrame(columns=['name'])
for letter in ['a', 'b', 'c']:
    for n in [0.01, 1, 100]:
        perform_alignment(f'{out}/transcriptomes.fasta', [
            f'SimMOSCA2/rna/{letter}/{n}/mt_{n}{letter}_R{fr}.fastq' for fr in ['1', '2']],
            f'{out}/alignments/mt_{n}{letter}', threads=15)
        readcounts = pd.merge(
            readcounts, pd.read_csv(f'{out}/alignments/mt_{n}{letter}.readcounts', sep='\t', header=None,
                                    names=['name', f'mt_{n}{letter}']), how='outer', on='name')
mt_cols = [f'mt_{n}{letter}' for letter in ['a', 'b', 'c'] for n in [0.01, 1, 100]]
readcounts[mt_cols] = readcounts[mt_cols].fillna(value=0.0).astype(int)
readcounts.to_csv(f'{out}/readcounts.tsv', sep='\t', index=False)