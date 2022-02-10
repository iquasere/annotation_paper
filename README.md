# UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets 

![Graphical Abstract](assets/graphical_abstract.jpg "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets ")

This repository was created to store the scripts used in the publication "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets". Main chapters are:
1. [Obtention of datasets](https://github.com/iquasere/annotation_paper#obtention-of-datasets)
2. [Installation of tools](https://github.com/iquasere/annotation_paper#installation-of-tools)
3. [Run tools for reference genomes](https://github.com/iquasere/annotation_paper#run-tools-for-reference-genomes)
4. [Results analysis](https://github.com/iquasere/annotation_paper#results-analysis)
5. [RNA-Seq simulation and quantification into readcounts](https://github.com/iquasere/annotation_paper#rna-seq-simulation-and-quantification-into-readcounts)
6. [Metagenomics analysis](https://github.com/iquasere/annotation_paper#metagenomics-analysis)
7. [Run KEGGCharter](https://github.com/iquasere/annotation_paper#run-keggcharter)
8. [Publication](https://github.com/iquasere/annotation_paper#publication)

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
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_trembl.fasta.gz
zcat uniprot_*.fasta.gz > resources_directory/uniprot.fasta
awk '{ print $1 }' resources_directory/uniprot.fasta > tmp.fa | mv tmp.fa resources_directory/uniprot.fasta
# uniprot had 225578953 sequences, so divided by 15 cores
awk -v size=15038597 -v pre=resources_directory/split_uniprot -v pad=2 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
   { print >> fname }
' resources_directory/uniprot.fasta
python annotation_paper/scripts/download_proteomes.py
awk '{print $1}' ann_paper/proteomes.fasta > tmp.fasta && mv tmp.fasta ann_paper/proteomes.fasta
grep '>' ann_paper/proteomes.fasta > ann_paper/ids.txt
python annotation_paper/scripts/uniprot_selection.py
```

## Run tools for reference genomes

```
conda activate ann_paper
prodigal -i ann_paper/genomes.fasta -a ann_paper/genes.fasta -o ann_paper/prodigal_out.txt
grep '>' ann_paper/genes.fasta | awk '{print substr($1, 2)"\t"$3"\t"$5}' > ann_paper/genes.tsv
awk '{print $1}' ann_paper/genes.fasta | sed 's/*//' > ann_paper/genes_trimmed.fasta
upimapi.py -i ann_paper/genes_trimmed.fasta -o ann_paper/upimapi_genomes -rd resources_directory -db uniprot -t 15
recognizer.py -f ann_paper/genes_trimmed.fasta -o ann_paper/recognizer_genomes -rd resources_directory -t 15
mkdir ann_paper/eggnog_mapper_genomes
emapper.py -i ann_paper/genes_trimmed.fasta -o ann_paper/eggnog_mapper_genomes/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data
conda activate mantis_env
python mantis run_mantis -i ann_paper/genes_trimmed.fasta -o ann_paper/mantis_genomes -c 15
conda deactivate
dfast_file_downloader.py --protein dfast --dbroot resources_directory/dfast
dfast_file_downloader.py --cdd Cog --hmm TIGR --dbroot resources_directory/dfast
dfast -g ann_paper/genomes.fasta -o ann_paper/dfast_genomes --cpu 15 --dbroot resources_directory/dfast
prokka --outdir ann_paper/prokka_genomes --metagenome --cpus 15 ann_paper/genomes.fasta
```
 
## Results analysis

First, download the proteomes corresponding to the reference genomes. Clean the proteomes, build DMND database from them, and align genes called with Prodigal to that database, obtaining the most likely identification of each called gene.
```
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
Finally, the main show: benchmarking the tools against each other.
```
python annotation_paper/scripts/tools_benchmark.py
```

## RNA-Seq simulation and quantification into readcounts

```
awk 'BEGIN{FS="\t"}{print $10}' annotation_paper/assets/simulated_taxa.tsv | tail -n +2 > ann_paper/assets/transcriptomes_links.txt 
wget -i ann_paper/assets/transcriptomes_links.txt -P ann_paper
gunzip ann_paper/*.gz
bash annotation_paper/scripts/clean_transcriptomes.sh
cat ann_paper/*.cds.all.fa > ann_paper/transcriptomes.fasta
mkdir ann_paper/simulated ann_paper/simulated/rna 
grep '>' ann_paper/transcriptomes.fasta | awk '{print substr($0, 2)}' > ann_paper/ids_to_map.txt
```
Then, go to UniProt's ID mapping [web service](https://www.uniprot.org/uploadlists/), pick the ```ids_to_map.txt``` file, select ```Ensembl Genomes Protein``` in the "From" field, and download results as "Tab-separated". Must contain the following columns: ```yourlist:...```, ```Pathway``` and ```Protein names```. ID map results saved as ```ann_paper/simulated/uniprotinfo.tsv```
```
python annotation_paper/scripts/rnaseq_simulation.py
```
RNA-Seq quantification followed preprocessing with MOSCA and read alignment with Bowtie2.
```
mkdir ann_paper/alignments
python annotation_paper/scripts/rnaseq_quantification.py
```

## Metagenomics analysis

```
python MOSCA/workflow/scripts/preprocess.py -i Datasets/EST6_S1_L001/EST6_S1_L001_R1_001.fastq.gz,Datasets/EST6_S1_L001/EST6_S1_L001_R2_001.fastq.gz -t 15 -d dna -o ann_paper/Preprocess_MG -rd resources_directory --avgqual 20 --minlen 100
python MOSCA/workflow/scripts/assembly.py -r ann_paper/Preprocess_MG/Trimmomatic/quality_trimmed_EST6_S1_L001_forward_paired.fq,ann_paper/Preprocess_MG/Trimmomatic/quality_trimmed_EST6_S1_L001_reverse_paired.fq -o ann_paper/Assembly_MG -t 15
conda activate ann_paper
prodigal -i ann_paper/Assembly_MG/contigs.fasta -a ann_paper/genes_MG.fasta
upimapi.py -i ann_paper/genes_MG.fasta -o ann_paper/upimapi_metagenome -rd resources_directory -db uniprot -t 15
recognizer.py -f ann_paper/genes_MG.fasta -o ann_paper/recognizer_metagenome -rd resources_directory -t 15
python mantis run_mantis -i ann_paper/genes_MG.fasta -o ann_paper/mantis_metagenome -c 15
mkdir ann_paper/eggnog_mapper_metagenome
emapper.py -i ann_paper/genes_MG.fasta -o ann_paper/eggnog_mapper_metagenome/eggnog_results --cpu 15 --data_dir ann_paper/eggnog_data
awk '{if ($0 ~ /^>/){split($0, a, "_"); print a[1]"_"a[2]} else print $0}' ann_paper/Assembly_MG/scaffolds.fasta > ann_paper/Assembly_MG/prokka_scaffolds.fasta  # For Prokka, "Contig ID must <= 37 chars long"
prokka --outdir ann_paper/prokka_metagenome --metagenome --cpus 15 ann_paper/Assembly_MG/prokka_scaffolds.fasta
dfast -g ann_paper/Assembly_MG/scaffolds.fasta -o ann_paper/dfast_metagenome --cpu 15 --force --dbroot resources_directory/dfast
grep 'ID=' ann_paper/dfast_metagenome/genome.gff > ann_paper/dfast_metagenome/genome_trimmed.gff
```

## Run KEGGCharter

```
python annotation_paper/scripts/prepare4keggcharter.py
keggcharter.py -f ann_paper/keggcharter_input.tsv -o ann_paper/keggcharter_genomes -rd resources_directory -tc "Taxonomic lineage (SPECIES)" -keggc "Cross-reference (KEGG)" -tcol mt_1,mt_100 -iq -not 7
keggcharter.py -f ann_paper/upimapi_metagenome/UPIMAPI_results.tsv -o ann_paper/keggcharter_genomes -rd resources_directory -tc "Taxonomic lineage (SPECIES)" -keggc "Cross-reference (KEGG)" -iq
```

## Submit assembly to webin

```
grep '>' ann_paper/Assembly_MG/scaffolds.fasta | awk 'BEGIN{FS="_";n=0;sum=0}{n+=1;sum+=$6}END{print sum/n}'
mamba install -c conda-forge -c bioconda ena-webin-cli -y
mkdir ann_paper/EBIsubmission
gunzip ann_paper/Assembly_MG/scaffolds.fasta
ena-webin-cli -context genome -userName username -password 'password' -manifest annotation_paper/assets/manifest.xml -outputDir ann_paper/EBIsubmission -inputDir ann_paper/Assembly_MG -submit
```

## Publication

This paper is still not published.