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
    index = reference.replace(f".{reference.split('.')[-1]}", '_index')
    if not check_bowtie2_index(index):
        print('INDEX files not found. Generating new ones')
        generate_mg_index(reference, index)
    else:
        print(f"INDEX was located at {index}")
    if not os.path.isfile(f'{basename}.log'):
        align_reads(
            reads, index, f'{basename}.sam', f'{basename}_bowtie2_report.txt', log=f'{basename}.log', threads=threads)
    else:
        print(f'{basename}.log was found!')
    run_pipe_command(
        f"""samtools view -F 260 -S {basename}.sam | cut -f 3 | sort | uniq -c | awk '{{printf("%s\\t%s\\n", $2, $1)}}'""",
        output=f'{basename}.readcounts')


out = 'ann_paper'


def triplicates_quantification():
    readcounts = pd.DataFrame(columns=['name'])
    for letter in ['a', 'b', 'c']:
        for n in [0.01, 1, 100]:
            perform_alignment(
                f'{out}/transcriptomes.fasta',
                [f'{out}/simulated/rna/Preprocess/Trimmomatic/quality_trimmed_rnaseq_{letter}{n}_{fr}_paired.fq'
                 for fr in ['forward', 'reverse']],
                f'{out}/alignments/mt_{n}{letter}', threads=15)
            readcounts = pd.merge(
                readcounts, pd.read_csv(
                    f'{out}/alignments/mt_{n}{letter}.readcounts', sep='\t', header=None,
                    names=['name', f'mt_{n}{letter}']), how='outer', on='name')
    mt_cols = [f'mt_{n}{letter}' for letter in ['a', 'b', 'c'] for n in [0.01, 1, 100]]
    readcounts[mt_cols] = readcounts[mt_cols].fillna(value=0.0).astype(int)
    factors = pd.read_csv('factors.txt', sep='\n', header=None)
    readcounts[mt_cols] = (readcounts[mt_cols] / factors[0].tolist()).astype(int)
    for n in [0.01, 1, 100]:
        readcounts[f'mt_{n}'] = readcounts[[f'mt_{n}{l}' for l in ['a', 'b', 'c']]].mean(axis=1)
    readcounts[[f'mt_{n}' for n in [0.01, 1, 100]]] = readcounts[[f'mt_{n}' for n in [0.01, 1, 100]]].astype(int)
    readcounts[['name'] + [f'mt_{n}' for n in [0.01, 1, 100]]].to_csv(f'{out}/readcounts.tsv', sep='\t', index=False)


def lipids_quantification():
    readcounts = pd.DataFrame(columns=['name'])
    for condition in ['1', '2']:
        perform_alignment(
            f'{out}/transcriptomes.fasta',
            [f'{out}/simulated/lipids/Preprocess/Trimmomatic/quality_trimmed_rnaseq_{condition}_{fr}_paired.fq'
             for fr in ['forward', 'reverse']],
            f'{out}/alignments/mt_{condition}', threads=15)
        readcounts = pd.merge(
            readcounts, pd.read_csv(
                f'{out}/alignments/mt_{condition}.readcounts', sep='\t', header=None,
                names=['name', f'mt_{condition}']), how='outer', on='name')
    mt_cols = [f'mt_{condition}' for condition in ['1', '2']]
    readcounts[mt_cols] = readcounts[mt_cols].fillna(value=0.0).astype(int)
    readcounts.to_csv(f'{out}/readcounts.tsv', sep='\t', index=False)


triplicates_quantification()
lipids_quantification()
