#!/bin/bash

# Assign IDs to genes identified with Prokka
diamond makedb --in ann_paper/"$1"/proteomes.fasta -d ann_paper/"$1"/proteomes.dmnd
diamond blastp -d ann_paper/"$1"/proteomes.dmnd -q ann_paper/"$1"/genes_trimmed.fasta -o ann_paper/"$1"/genes.blast --very-sensitive -p 15
awk '{print $1}' ann_paper/"$1"/prokka_genomes/PROKKA_"$2".faa > ann_paper/"$1"/prokka_genomes/PROKKA_"$2"_trimmed.faa
diamond blastp -d ann_paper/"$1"/proteomes.dmnd -q ann_paper/"$1"/prokka_genomes/PROKKA_"$2"_trimmed.faa -o ann_paper/"$1"/prokka_genomes/PROKKA_"$2".blast --very-sensitive -p 15
awk '{print $1}' ann_paper/"$1"/dfast_genomes/protein.faa > ann_paper/"$1"/dfast_genomes/protein_trimmed.faa
diamond blastp -d ann_paper/"$1"/proteomes.dmnd -q ann_paper/"$1"/dfast_genomes/protein_trimmed.faa -o ann_paper/"$1"/dfast_genomes/protein.blast --very-sensitive -p 15
# get information for all IDs from the proteomes
grep '>' ann_paper/"$1"/proteomes.fasta | awk '{print substr($0, 2)}' | tr '\n' ',' | sed 's/.$//' | upimapi.py -ot ann_paper/"$1"/uniprotinfo.tsv --no-annotation -rd resources_directory -cols "Entry&Entry name&EC number" -dbs "Conserved Domains Database&Pfam protein domain database&TIGRFAMs; a protein family database&Simple Modular Architecture Research Tool; a protein domain database&evolutionary genealogy of genes: Non-supervised Orthologous Groups"
grep 'ID=' ann_paper/"$1"/dfast_genomes/genome.gff > ann_paper/"$1"/dfast_genomes/genome_trimmed.gff
