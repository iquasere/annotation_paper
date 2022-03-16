# UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets 

![Graphical Abstract](assets/graphical_abstract.jpg "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets ")

This repository was created to store the scripts used in the publication "UPIMAPI, reCOGnizer and KEGGCharter: bioinformatics tools for functional annotation and visualization of (meta)-omics datasets". Main chapters are:
1. [Installation of tools](https://github.com/iquasere/annotation_paper#installation-of-tools)
2. [Retrieval of random selections of queries and construction of independent UniProts](https://github.com/iquasere/annotation_paper#retrieval-of-random-selections-of-queries-and-construction-of-independent-uniprots)
3. [Run tools](https://github.com/iquasere/annotation_paper#run-tools)
4. [Results analysis](https://github.com/iquasere/annotation_paper#results-analysis)
5. [RNA-Seq simulation and quantification into readcounts](https://github.com/iquasere/annotation_paper#rna-seq-simulation-and-quantification-into-readcounts)
6. [Metagenomics analysis](https://github.com/iquasere/annotation_paper#metagenomics-analysis)
7. [Run KEGGCharter](https://github.com/iquasere/annotation_paper#run-keggcharter)
8. [Submit assembly to webin](https://github.com/iquasere/annotation_paper#submit-assembly-to-webin)
9. [Publication](https://github.com/iquasere/annotation_paper#publication)


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

## Retrieval of random selections of queries and construction of independent UniProts

Download UniProt database and split it for 15 cores.
```
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_trembl.fasta.gz
zcat uniprot_*.fasta.gz > resources_directory/uniprot.fasta
rm uniprot_*.fasta.gz
awk '{ print $1 }' resources_directory/uniprot.fasta > tmp.fa & mv tmp.fa resources_directory/uniprot.fasta
# uniprot had 225578953 sequences, so divided by 15 cores
awk -v size=15038597 -v pre=resources_directory/split_uniprot -v pad=2 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
   { print >> fname }
' resources_directory/uniprot.fasta
```
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
python annotation_paper/scripts/download_proteomes.py
awk '{print $1}' ann_paper/proteomes.fasta > tmp.fasta && mv tmp.fasta ann_paper/proteomes.fasta
grep '>' ann_paper/proteomes.fasta > ann_paper/ids.txt
python annotation_paper/scripts/uniprot_selection.py
```
Prepare databases (UniProts without query sequences) and queries for multiple iterations.
```
python annotation_paper/scripts/uniprot_selection_workflow.py
bash annotation_paper/scripts/preprocess_queries.sh
```
After this, edit all ```genomes_fixed.fasta``` to add "_N" (where N is iteration 1,2,...) to ">Chromosome"s (I did it by hand, use Notepad++). Finally, run the tools:

## Run tools

```
bash annotation_paper/scripts/run_tools.sh ""
bash annotation_paper/scripts/run_tools.sh "/first_group"
bash annotation_paper/scripts/run_tools.sh "/second_group"
bash annotation_paper/scripts/run_tools.sh "/third_group"
bash annotation_paper/scripts/run_tools.sh "/fourth_group"
bash annotation_paper/scripts/run_tools.sh "/fifth_group"
```
 
## Results analysis

First, download the proteomes corresponding to the reference genomes. Clean the proteomes, build DMND database from them, and align genes called with Prodigal to that database, obtaining the most likely identification of each called gene (all this happens in the script ```preprocess_results.sh``.
```
bash annotation_paper/scripts/preprocess_results.sh "" "01032022"
bash annotation_paper/scripts/preprocess_results.sh "/first_group" "03052022"
bash annotation_paper/scripts/preprocess_results.sh "/second_group" "02192022"
bash annotation_paper/scripts/preprocess_results.sh "/third_group" "02252022"
bash annotation_paper/scripts/preprocess_results.sh "/fourth_group" "03072022"
bash annotation_paper/scripts/preprocess_results.sh "/fifth_group" "03112022"
```
Different evalues are tested for UPIMAPI and reCOGnizer.
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