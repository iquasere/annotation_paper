#!/bin/bash

echo "Building selected UniProt"
cat ann_paper"$1"/uniprot_*.fasta > resources_directory/uniprot.fasta
rm resources_directory/uniprot.dmnd

echo "Building genomes.fasta"
cat ann_paper"$1"/*.fna > ann_paper"$1"/genomes.fasta
awk '{print $1}' ann_paper"$1"/genomes.fasta > ann_paper"$1"/genomes_fixed.fasta

echo "Gene calling"
prodigal -i ann_paper"$1"/genomes_fixed.fasta -a ann_paper"$1"/genes.fasta -p meta > ann_paper"$1"/prodigal_out.txt
awk '{print $1}' ann_paper"$1"/genes.fasta | sed 's/*//' > ann_paper"$1"/genes_trimmed.fasta

echo "Running UPIMAPI"
upimapi.py -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/upimapi_genomes -rd resources_directory -db uniprot -t 15

echo "Running reCOGnizer"
recognizer.py -f ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/recognizer_genomes -rd resources_directory -t 15

echo "Running eggNOG-mapper"
mkdir ann_paper"$1"/eggnog_mapper_genomes
emapper.py -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/eggnog_mapper_genomes/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data

echo "Running DFAST"
dfast -g ann_paper"$1"/genomes.fasta -o ann_paper"$1"/dfast_genomes --cpu 15 --dbroot resources_directory/dfast

echo "Running Prokka"
prokka --outdir ann_paper"$1"/prokka_genomes --metagenome --cpus 15 --centre X --compliant ann_paper"$1"/genomes_fixed.fasta

echo "Running mantis"
eval "$(conda shell.bash hook)"
conda activate mantis_env
python mantis run_mantis -i ann_paper"$1"/genes_trimmed.fasta -o ann_paper"$1"/mantis_genomes -c 15
conda deactivate

echo "Removing artifacts"
rm ann_paper"$1"/uniprot_*.fasta
