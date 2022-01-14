# UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for fast functional annotation and visualization of (meta)-omics datasets 

This repository was created to store the scripts used in the publication "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for fast functional annotation and visualization of (meta)-omics datasets 
". These scripts run the following steps:
1. [Obtention of datasets](https://github.com/iquasere/annotation_paper#obtention-of-datasets)
2. [Installation of tools](https://github.com/iquasere/annotation_paper#installation-of-tools)
3. [Execution of tools](https://github.com/iquasere/annotation_paper#execution-of-tools)
4. [Results analysis](https://github.com/iquasere/annotation_paper#results-analysis)

## Obtention of datasets

Obtain genomes for 7 prokaryotes (archaea and bacteria), and join them in a single file.
Simulate RNA-Seq data for those prokaryotes for three different conditions: over-, normal, and underexpression of a collection of genes.
```
git clone https://github.com/iquasere/annotation_paper.git
mkdir ann_paper
mkdir ann_paper/assets
awk 'BEGIN{FS="\t"}{print $9}' annotation_paper/assets/simulated_taxa.tsv | tail -n +2 > ann_paper/assets/genomes_links.txt 
wget -i ann_paper/assets/genomes_links.txt -P ann_paper
gunzip ann_paper/*.gz
awk 'BEGIN{FS=":"}{if ($0 ~ /^>/) {print ">"$3} else {print $0}}' ann_paper/*.fa >> ann_paper/genomes.fasta
rm ann_paper/*.fa
```

## Installation of tools

Requires **mamba** installed in the current environment. Mamba can be installed with ```conda install -c conda-forge mamba```
```
mamba create -c conda-forge -c bioconda upimapi recognizer keggcharter dfast prokka eggnog-mapper bowtie2 bioconductor-deseq2 r-pheatmap r-rcolorbrewer -n ann_paper -y
mkdir ann_paper/eggnog_data
conda activate ann_paper
download_eggnog_data.py --data_dir ann_paper/eggnog_data -y
git clone https://github.com/PedroMTQ/mantis.git
conda env create -f mantis/mantis_env.yml
conda activate mantis_env
python mantis setup_databases
mv annotation_paper/assets/MANTIS.config mantis/MANTIS.config
git clone https://github.com/iquasere/MOSCA.git
```

## RNA-Seq simulation

```
awk 'BEGIN{FS="\t"}{print $10}' annotation_paper/assets/simulated_taxa.tsv | tail -n +2 > ann_paper/assets/transcriptomes_links.txt 
wget -i ann_paper/assets/transcriptomes_links.txt -P ann_paper
gunzip ann_paper/*.gz
bash annotation_paper/scripts/clean_transcriptomes.bash
cat ann_paper/*.cds.all.fa >> ann_paper/transcriptomes.fasta
mkdir ann_paper/simulated ann_paper/simulated/rna 
grep '>' ann_paper/transcriptomes.fasta | awk '{print substr($0, 2)}' > ann_paper/ids_to_map.txt
```
Then, go to UniProt's ID mapping [web service](https://www.uniprot.org/uploadlists/), pick the ```ids_to_map.txt``` file, select ```Ensembl Genomes Protein``` in the "From" field, and download results as "Tab-separated". Must contain the following columns: ```yourlist:...```, ```Pathway``` and ```Protein names```. ID map results saved as ```ann_paper/simulated/uniprotinfo.tsv```
```
python annotation_paper/scripts/rnaseq_simulation.py
```

## RNA-Seq quantification

Quantification followed preprocessing with MOSCA and read alignment with Bowtie2.
```
mkdir ann_paper/alignments
python annotation_paper/scripts/rnaseq_quantification.py
```

## Run tools for reference genomes

```
conda activate ann_paper
prodigal -i ann_paper/genomes.fasta -a ann_paper/genes.fasta -o ann_paper/prodigal_out.txt
grep '>' ann_paper/genes.fasta | awk '{print substr($1, 2)"\t"$3"\t"$5}' > ann_paper/genes.tsv
awk '{print $1}' ann_paper/genes.fasta | sed 's/*//' > ann_paper/genes_trimmed.fasta
upimapi.py -i ann_paper/genes_trimmed.fasta -o ann_paper/upimapi_genomes -rd resources_directory --evalue 0.1 -mts 20 -db uniprot -t 15
recognizer.py -f ann_paper/genes_trimmed.fasta -o ann_paper/recognizer_genomes -rd resources_directory -mts 20 -t 15 --tax-file ann_paper/upimapi_genomes/UPIMAPI_results.tsv --tax-col "Taxonomic lineage IDs (SPECIES)" --protein-id-col qseqid
emapper.py -i ann_paper/genes_trimmed.fasta -o ann_paper/eggnog_mapper_genomes/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data
mkdir ann_paper/eggnog_mapper_genomes
python mantis run_mantis -i ann_paper/genes_trimmed.fasta -o ann_paper/mantis_genomes -c 15
dfast_file_downloader.py --protein dfast --dbroot resources_directory/dfast
dfast_file_downloader.py --cdd Cog --hmm TIGR --dbroot resources_directory/dfast
dfast -g ann_paper/genomes.fasta -o ann_paper/dfast_genomes --threshold 0,75,75,1e-3 --cpu 15 --force --dbroot resources_directory/dfast
prokka --outdir ann_paper/prokka_genomes --metagenome --evalue 1e-3 --cpus 15 ann_paper/genomes.fasta

keggcharter.py -f paper_pipeline/keggcharter_output_recognizer/KEGGCharter_results.tsv -o paper_pipeline/keggcharter_output_recognizer -rd resources_directory -tcol mt_0.01a,mt_1a,mt_100a,mt_0.01b,mt_1b,mt_100b,mt_0.01c,mt_1c,mt_100c -ecc "EC number_x" -gcol mg -it "meta-omics community"
keggcharter.py -o ann_paper/keggcharter_proteome -f ann_paper/keggcharter_input.tsv -rd resources_directory -tc "Taxonomic lineage (SPECIES)" -keggc "Cross-references (KEGG)" -tcol mt_0.01a,mt_1a,mt_100a,mt_0.01b,mt_1b','mt_100b,mt_0.01c,mt_1c,mt_100c -gcol mg
keggcharter.py -f paper_pipeline/keggcharter_output_recognizer/KEGGCharter_results.tsv -o paper_pipeline/keggcharter_output_recognizer -rd resources_directory -tcol mt_0.01a,mt_1a,mt_100a,mt_0.01b,mt_1b,mt_100b,mt_0.01c,mt_1c,mt_100c -ecc "EC number_x" -gcol mg -it "meta-omics community" --resume  -rd resources_directory -tcol mt_0.01a,mt_1a,mt_100a,mt_0.01b,mt_1b,mt_100b,mt_0.01c,mt_1c,mt_100c -ecc "EC number_x"
```
 
## Run tools for metagenome and metatranscriptome

```
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1613-MiSeqKapa/4478-DNA-S1613-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1613-MiSeqKapa/4478-DNA-S1613-MiSeqKapa_R2.fastq.gz -t 15 -d dna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1616-MiSeqKapa/4478-DNA-S1616-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1616-MiSeqKapa/4478-DNA-S1616-MiSeqKapa_R2.fastq.gz -t 15 -d dna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1618-MiSeqKapa/4478-DNA-S1618-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1618-MiSeqKapa/4478-DNA-S1618-MiSeqKapa_R2.fastq.gz -t 15 -d dna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1611-MiSeqKapa/4478-R1-1-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1611-MiSeqKapa/4478-R1-1-MiSeqKapa_R2.fastq.gz -t 15 -d mrna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1613-MiSeqKapa/4478-R2-1-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1613-MiSeqKapa/4478-R2-1-MiSeqKapa_R2.fastq.gz -t 15 -d mrna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1616-MiSeqKapa/4478-R3-1-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1616-MiSeqKapa/4478-R3-1-MiSeqKapa_R2.fastq.gz -t 15 -d mrna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/preprocess.py -i Datasets/4478-DNA-S1618-MiSeqKapa/4478-R4-1-MiSeqKapa_R1.fastq.gz,Datasets/4478-DNA-S1618-MiSeqKapa/4478-R4-1-MiSeqKapa_R2.fastq.gz -t 15 -d mrna -o ann_paper/metagenome/Preprocess -rd resources_directory --avgqual 20 --minlen 100
```
 
## Results analysis

First, download the proteomes corresponding to the reference genomes. Clean the proteomes, build DMND database from them, and align genes called with Prodigal to that database, obtaining the most likely identification of each called gene.
```
python annotation_paper/scripts/download_proteomes.py
awk '{print $1}' ann_paper/proteomes.fasta > tmp.fasta && mv tmp.fasta ann_paper/proteomes.fasta
# Assign IDs to genes identified with Prokka
diamond makedb --in ann_paper/proteomes.fasta -d ann_paper/proteomes.dmnd
diamond blastp -d ann_paper/proteomes.dmnd -q ann_paper/genes_trimmed.fasta -o ann_paper/genes.blast --very-sensitive -p 15
awk '{print $1}' ann_paper/prokka_genomes/PROKKA_01032022.faa > tmp.faa | mv tmp.faa ann_paper/prokka_genomes/PROKKA_01032022.faa
diamond blastp -d ann_paper/proteomes.dmnd -q ann_paper/prokka_genomes/PROKKA_01032022.faa -o ann_paper/prokka_genomes/PROKKA_01032022.blast --very-sensitive -p 15
awk '{print $1}' ann_paper/dfast_genomes/protein.faa > tmp.faa | mv tmp.faa ann_paper/dfast_genomes/protein.faa
diamond blastp -d ann_paper/proteomes.dmnd -q ann_paper/dfast_genomes/protein.faa -o ann_paper/dfast_genomes/protein.blast --very-sensitive -p 15
# get information for all IDs from the proteomes
grep '>' ann_paper/proteomes.fasta | awk '{print substr($0, 2)}' | tr '\n' ',' | sed 's/.$//' | upimapi.py -ot ann_paper/uniprotinfo.tsv --no-annotation --no-local-mapping -rd resources_directory -cols "Entry&Entry name&EC number" -dbs "Conserved Domains Database&Pfam protein domain database&TIGRFAMs; a protein family database&Simple Modular Architecture Research Tool; a protein domain database&evolutionary genealogy of genes: Non-supervised Orthologous Groups"
grep 'ID=' ann_paper/dfast_genomes/genome.gff > ann_paper/dfast_genomes/genome_trimmed.gff
```
Then, different evalues are tested for each tool.
```
python annotation_paper/scripts/evalue_benchmark.py
```

## Metagenomics analysis

```
python MOSCA/workflow/scripts/preprocess.py -i Datasets/EST6_S1_L001/EST6_S1_L001_R1_001.fastq.gz,Datasets/EST6_S1_L001/EST6_S1_L001_R2_001.fastq.gz -t 15 -d dna -o ann_paper/Preprocess_MG -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/assembly.py -r ann_paper/Preprocess_MG/Trimmomatic/quality_trimmed_EST6_S1_L001_forward_paired.fq,ann_paper/Preprocess_MG/Trimmomatic/quality_trimmed_EST6_S1_L001_reverse_paired.fq -o ann_paper/Assembly_MG -t 15
conda activate ann_paper
prodigal -i ann_paper/Assembly_MG/contigs.fasta -a ann_paper/genes_MG.fasta


upimapi.py -i ann_paper/genes_MG.fasta -o ann_paper/upimapi_metagenome -rd resources_directory -db uniprot -t 15
02d01h01m04s
recognizer.py -f ann_paper/genes_MG.fasta -o ann_paper/recognizer_metagenome -rd resources_directory -mts 20 -t 15

python mantis run_mantis -i ann_paper/genes_MG.fasta -o ann_paper/mantis_metagenome -c 15
7h28m15s
mkdir ann_paper/eggnog_mapper_metagenome
emapper.py -i ann_paper/genes_MG.fasta -o ann_paper/eggnog_mapper_metagenome/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data
8330 s
awk '{if ($0 ~ /^>/){split($0, a, "_"); print a[1]"_"a[2]} else print $0}' ann_paper/Assembly_MG/scaffolds.fasta > ann_paper/Assembly_MG/prokka_scaffolds.fasta  # For Prokka "Contig ID must <= 37 chars long"
prokka --outdir ann_paper/prokka_metagenome --metagenome --cpus 15 ann_paper/Assembly_MG/prokka_scaffolds.fasta
10h49m41s
dfast -g ann_paper/Assembly_MG/scaffolds.fasta -o ann_paper/dfast_metagenome --cpu 15 --force --dbroot resources_directory/dfast


```

## Publication

This paper is still not published.