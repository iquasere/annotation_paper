#!/bin/bash

gunzip ann_paper"$1"/*.gz
cat ann_paper"$1"/*.fna > ann_paper"$1"/genomes.fasta
cat ann_paper"$1"/uniprot_*.fasta > resources_directory/uniprot.fasta
awk '{print $1}' ann_paper"$1"/genomes.fasta > ann_paper"$1"/genomes_fixed.fasta
prodigal -i ann_paper"$1"/genomes_fixed.fasta -a ann_paper"$1"/genes.fasta
awk '{print $1}' ann_paper"$1"/genes.fasta | sed 's/*//' > ann_paper"$1"/genes_trimmed.fasta
upimapi.py -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/upimapi_genomes -rd resources_directory -db uniprot -t 15
recognizer.py -f ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/recognizer_genomes -rd resources_directory -t 15
mkdir ann_paper"$1"/eggnog_mapper_genomes
emapper.py -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/eggnog_mapper_genomes/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data
conda activate mantis_env
python mantis run_mantis -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/mantis_genomes -c 15
conda deactivate
dfast -g ann_paper"$1"/genomes.fasta -o ann_paper"$1"/dfast_genomes --cpu 15 --dbroot resources_directory/dfast
tr ' ' '_' ann_paper"$1"/genomes.fasta > ann_paper"$1"/genomes_fixed.fasta
prokka --outdir ann_paper"$1"/prokka_genomes --metagenome --cpus 15 --centre X --compliant ann_paper"$1"/genomes_fixed.fasta
rm ann_paper"$1"/uniprot_*.fasta
