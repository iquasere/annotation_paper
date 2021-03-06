import pandas as pd
from pathlib import Path
import os
from paper_utils import run_command, run_pipe_command, get_fasta_ids
from tqdm import tqdm

out = 'ann_paper'

def set_abundance_mt(output_dir, simulated, profiles, uniprotinfo, factor=1, base=0):
    simulated = simulated[simulated['Transcriptome link'].notnull()].reset_index()
    uniprotinfo = pd.read_csv(uniprotinfo, sep='\t', index_col=-1)
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
                if type(info['Pathway']) == str and path in info['Pathway']:
                    handler.write(f'{rna}\t{base + factor * abundance}\n')
                    found = True
            if not found:
                if 'ATP synthase' in info['Protein names']:
                    handler.write(f'{rna}\t{base + factor * abundance}\n')
                else:
                    handler.write(f'{rna}\t{base + abundance}\n')


def get_ec_factor(ecs, condition, short_chain_ecs, long_chain_ecs):
    for ec in ecs:
        if (ec in short_chain_ecs and condition == '1') or (ec in long_chain_ecs and condition == '2'):
            return 100
        elif ec in short_chain_ecs and condition == '2':
            return 50
        elif ec in long_chain_ecs and condition == '1':
            return 0
        else:
            return 1


def set_abundance_lipids(output_dir, simulated, uniprotinfo, condition, base=0):
    short_chain_ecs = ["1.3.3.6", "1.3.8.5", "1.3.8.1", "1.3.8.7", "4.2.1.17", "1.1.1.35", "2.3.1.16", "2.3.1.9"]
    long_chain_ecs = [
        "1.3.3.6", "1.3.8.5", "1.3.8.1", "1.3.8.7", "4.2.1.17", "1.1.1.35", "2.3.1.16", "2.3.1.9", "6.2.1.3", "1.3.8.8",
        "4.2.1.74", "4.2.1.17", "1.1.1.211"]
    simulated = simulated[simulated['Transcriptome link'].notnull()].reset_index()
    uniprotinfo = pd.read_csv(uniprotinfo, sep='\t', index_col=-1)
    uniprotinfo['EC number'] = uniprotinfo['EC number'].fillna(value='').str.split('; ')
    handler = open(f'{output_dir}/abundance{condition}.config', 'w')
    print(f'Setting expressions for {len(simulated)} organisms.')
    for i in range(len(simulated)):
        print(f'Dealing with {simulated.iloc[i]["Species"]}')
        rna_file = f'{out}/{simulated.iloc[i]["Transcriptome link"].split("/")[-1].split(".gz")[0]}'
        abundance = simulated.iloc[i]['Abundance']
        for rna in tqdm(get_fasta_ids(rna_file)):
            ide = rna.split()[0]
            if ide not in uniprotinfo.index:
                continue
            info = uniprotinfo.loc[ide]
            factor = get_ec_factor(info['EC number'], condition, short_chain_ecs, long_chain_ecs)
            handler.write(f'{rna}\t{base + factor * abundance}\n')


def run_grinder(reference_file, abundance_file, output_dir, coverage_fold=None, seed=None):
    command = (
        f'/usr/bin/perl {os.path.expanduser("~")}/anaconda3/bin/grinder -fastq_output 1 -qual_levels 30 10 -read_dist '
        f'151 -insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8 -mate_orientation FR -output_dir {output_dir} '
        f'-reference_file {reference_file} -abundance_file {abundance_file} -coverage_fold {coverage_fold} '
        f'-random_seed {seed}')
    run_command(command)


def divide_fq(file, output1, output2):
    run_command(f'bash annotation_paper/assets/unmerge-paired-reads.sh {file} {output1} {output2}')


def simulate_rnaseq_triplicates():
    info_df = pd.read_csv('annotation_paper/assets/simulated_taxa.tsv', sep='\t')
    sim_dir = f'{out}/simulated/rna'
    profiles_dict = {
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
    s = 13
    for n in [100, 1, 0.01]:
        set_abundance_mt(sim_dir, info_df, profiles_dict, f'{out}/simulated/uniprotinfo.tsv', factor=n)
        for letter in ['a', 'b', 'c']:

            Path(f'{sim_dir}/{letter}/{n}').mkdir(parents=True, exist_ok=True)

            run_grinder(
                f'{out}/transcriptomes.fasta', f'{sim_dir}/abundance{n}.config',
                f'{sim_dir}/{letter}/{n}', seed=s, coverage_fold='25')

            s += 1
            # split interleaved reads into forward and reverse
            divide_fq(
                f'{sim_dir}/{letter}/{n}/grinder-reads.fastq',
                f'{sim_dir}/{letter}/{n}/rnaseq_{letter}{n}_R1.fastq',
                f'{sim_dir}/{letter}/{n}/rnaseq_{letter}{n}_R2.fastq')

            for r in ['R1', 'R2']:
                # convert @1/1 ... to 1 1
                run_pipe_command(
                    f"awk '{{print $1}}' {sim_dir}/{letter}/{n}/rnaseq_{letter}{n}_{r}.fastq | sed 's/\\// /g' > tmp.fastq "
                    f"&& mv tmp.fastq {sim_dir}/{letter}/{n}/rnaseq_{letter}{n}_{r}.fastq")

            # preprocessing
            run_command(
                'python MOSCA/workflow/scripts/preprocess.py '
                '-i {2}/{0}/{1}/rnaseq_{0}{1}_R1.fastq,{2}/{0}/{1}/rnaseq_{0}{1}_R2.fastq -t 14 -d mrna '
                '-o {2}/Preprocess -rd resources_directory --avgqual 20 --minlen 100'.format(letter, n, sim_dir))


def simulate_lipids_degradation():
    info_df = pd.read_csv('annotation_paper/assets/simulated_taxa.tsv', sep='\t')
    sim_dir = f'{out}/simulated/lipids'
    # this might be needed: download grinder from ;  perl Makefile.PL; make; sudo make install; cpan Bio::Perl
    s = 13
    for condition in ['1', '2']:
        set_abundance_lipids(sim_dir, info_df, f'{out}/simulated/uniprotinfo.tsv', condition)
        Path(f'{sim_dir}/{condition}').mkdir(parents=True, exist_ok=True)
        run_grinder(
            f'{out}/transcriptomes.fasta', f'{sim_dir}/abundance{condition}.config',
            f'{sim_dir}/{condition}', seed=s, coverage_fold='25')
        s += 1
        # split interleaved reads into forward and reverse
        divide_fq(
            f'{sim_dir}/{condition}/grinder-reads.fastq',
            f'{sim_dir}/{condition}/rnaseq_{condition}_R1.fastq',
            f'{sim_dir}/{condition}/rnaseq_{condition}_R2.fastq')

        for r in ['R1', 'R2']:
            # convert @1/1 ... to 1 1
            run_pipe_command(
                f"awk '{{print $1}}' {sim_dir}/{condition}/rnaseq_{condition}_{r}.fastq | sed 's/\\// /g' > tmp.fastq "
                f"&& mv tmp.fastq {sim_dir}/{condition}/rnaseq_{condition}_{r}.fastq")

        # preprocessing
        run_command(
            'python MOSCA/workflow/scripts/preprocess.py '
            '-i {1}/{0}/rnaseq_{0}_R1.fastq,{1}/{0}/rnaseq_{0}_R2.fastq -t 14 -d mrna '
            '-o {1}/Preprocess -rd resources_directory --avgqual 20 --minlen 100'.format(condition, sim_dir))


simulate_rnaseq_triplicates()
simulate_lipids_degradation()
