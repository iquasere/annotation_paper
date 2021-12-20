import pandas as pd
from pathlib import Path
import os
from paper_utils import run_command, get_fasta_ids


def set_abundance_mt(out, simulated, profiles, uniprotinfo, output_dir, factor=1, base=0):
    simulated = simulated[simulated['Transcriptome link'].notnull()].reset_index()
    uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t', index_col = 0)
    handler = open(f'{output_dir}/abundance{factor}.config', 'w')
    print(f'Setting expressions for {len(simulated)} organisms.')
    for i in range(len(simulated)):
        print(f'Dealing with {simulated.iloc[i]["Species"]}')
        rna_file = f'{out}/{simulated.iloc[i]["Transcriptome link"].split("/")[-1].split(".gz")[0]}'
        abundance = simulated.iloc[i]['Abundance']
        profile = profiles[simulated.iloc[i]['Profile']]
        for rna in get_fasta_ids(rna_file):
            found = False
            ide = rna.split()[0]
            if ide not in uniprotinfo.index:
                continue
            info = uniprotinfo.loc[ide]
            for path in profile['Pathway'].keys():
                if path in info['Pathway']:
                    print(rna)
                    handler.write(f'{rna}\t{base + factor * abundance}\n')
                    found = True
            if not found:
                if 'ATP synthase' in info['Protein names']:
                    handler.write(f'{rna}\t{base + factor * abundance}\n')
                else:
                    handler.write(f'{rna}\t{base + abundance}\n')


def run_grinder(reference_file, abundance_file, output_dir, coverage_fold = None, seed = None):
    command = f'/usr/bin/perl ~/anaconda3/bin/grinder -fastq_output 1 -qual_levels 30 10 -read_dist 151 ' \
              f'-insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8 -mate_orientation FR -output_dir {output_dir} ' \
              f'-reference_file {reference_file} -abundance_file {abundance_file} -coverage_fold {coverage_fold} ' \
              f'-random_seed {seed}'
    run_command(command)


def divide_fq(file, output1, output2):
    run_command(f'bash functional_annotation_publication/assets/unmerge-paired-reads.sh {file} {output1} {output2}')


out = 'ann_paper'
uniprotinfo = pd.read_csv(f'{out}/simulated/uniprotinfo.tsv', sep='\t')
uniprotinfo.columns = ['Ensembl_IDs'] + uniprotinfo.columns.tolist()[1:]
letters = ['a', 'b', 'c']
factors = ['0.01', '1', '100']
simulated = pd.read_csv('functional_annotation_publication/assets/simulated_taxa.tsv', sep='\t')
sim_dir = f'{out}/simulated/rna'
profiles = {
    'archaea_co2': {
        'Pathway': {
            'One-carbon metabolism; methanogenesis from CO(2)': 1,
            'Cofactor biosynthesis': 1},
        'Protein names': {'ATP synthase': 1}},
    'archaea_acetate': {
        'Pathway': {
            'One-carbon metabolism; methanogenesis from acetate': 1,
            'Cofactor biosynthesis': 1},
        'Protein names': {'ATP synthase': 1}},
    'bacteria': {
        'Pathway': {
            'lipid metabolism': 1,
            'Cofactor biosynthesis': 1,
            'Metabolic intermediate biosynthesis': 1},
        'Protein names': {'ATP synthase': 1}}}

# this might be needed: download grinder from ;  perl Makefile.PL; make; sudo make install; cpan Bio::Perl
for factor in [100, 1, 0.01]:
    set_abundance_mt(out, simulated, profiles, f'{out}/simulated/uniprotinfo.tsv', sim_dir, factor=factor)
    seed = 13
    for letter in ['a', 'b', 'c']:
        Path(f'{sim_dir}/{letter}/{factor}').mkdir(parents=True, exist_ok=True)
        run_grinder(
            f'{out}/transcriptomes.fasta', f'{sim_dir}/abundance{factor}.config',
            f'{sim_dir}/{letter}/{factor}', seed=seed, coverage_fold='25')
        seed += 1
        # grinder produces interleaved reads, here we break them into forward and reverse reads
        divide_fq(
            f'{sim_dir}/{letter}/{factor}/grinder-reads.fastq',
            f'{sim_dir}/{letter}/{factor}/rnaseq_{letter}{factor}reads_R1.fastq',
            f'{sim_dir}/{letter}/{factor}/rnaseq_{letter}{factor}reads_R2.fastq')
        # preprocessing
        run_command(os.path.expanduser(
            'python ~/anaconda3/envs/mosca/share/MOSCA/scripts/preprocess.py '
            '-i paper_pipeline/sim_formicicum/rna/{0}/{1}/rnaseq_{0}{1}reads_R1.fastq,'
            'paper_pipeline/sim_formicicum/rna/{0}/{1}/rnaseq_{0}{1}reads_R2.fastq -t 14 -d mrna '
            '-o paper_pipeline/formicicum_mt/Preprocess/ -adaptdir resources_directory/adapters/ '
            '-rrnadbs resources_directory/rRNA_databases/ --avgqual 20 --minlen 100'.format(letter, factor)))
