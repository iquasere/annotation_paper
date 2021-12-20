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
awk 'BEGIN{FS="\t"}{print $9}' annotation_paper/assets/simulated_taxa.tsv > ann_paper/genomes_links.txt 
wget -i ann_paper/genomes_links.txt -P ann_paper
gunzip ann_paper/*.gz
awk 'BEGIN{FS=":"}{if ($0 ~ /^>/) {print ">"$3} else {print $0}}' ann_paper/*.fa >> ann_paper/genomes.fasta
rm ann_paper/*.fa
```

## Installation of tools

Requires **mamba** installed in the current environment. Mamba can be installed with ```conda install -c conda-forge mamba```
```
mamba create -c conda-forge -c bioconda upimapi recognizer keggcharter dfast prokka eggnog-mapper -n ann_paper -y
mkdir ann_paper/eggnog_data
conda activate ann_paper
download_eggnog_data.py --data_dir ann_paper/eggnog_data -y
git clone https://github.com/PedroMTQ/mantis.git
conda env create -f mantis/mantis_env.yml
conda activate mantis_env
python mantis setup_databases
```

## Execution of tools

```
conda activate ann_paper
prodigal -i ann_paper/genomes.fasta -a ann_paper/genes.fasta -o ann_paper/prodigal_out.txt
grep '>' ann_paper/genes.fasta | awk '{print $1"\t"$3"\t"$5}' > genes.tsv
awk '{print $1}' ann_paper/genes.fasta | sed 's/*//' - > ann_paper/genes_trimmed.fasta
upimapi.py -i ann_paper/genes_trimmed.fasta -o ann_paper/upimapi_genomes -rd resources_directory --evalue 0.1 -mts 20 -db uniprot -t 15
recognizer.py -f ann_paper/genes_trimmed.fasta -o ann_paper/recognizer_proteomes -rd resources_directory --evalue 0.1 -mts 20 -t 15 --tax-file ann_paper/upimapi_proteomes/UPIMAPI_results.tsv --tax-col "Taxonomic lineage (SPECIES)" --protein-id-col qseqid
emapper.py -i ann_paper/genes_trimmed.fasta -o ann_paper/eggnog_mapper_proteomes/eggnog_results --cpu 15
python mantis run_mantis -t ann_paper/genes_trimmed.fasta -o ann_paper/mantis_proteomes -c 15
```

## RNA-Seq simulation

```
awk 'BEGIN{FS="\t"}{print $10}' annotation_paper/assets/simulated_taxa.tsv > ann_paper/transcriptomes_links.txt 
wget -i ann_paper/transcriptomes_links.txt -P ann_paper
gunzip ann_paper/*.gz
awk '{print $1}' ann_paper/Geobacter_sulfurreducens.ASM96103v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Geobacter_sulfurreducens.ASM96103v1.cds.all.fa
awk '{print $1}' ann_paper/Methanobacterium_formicicum.DSM1535.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Methanobacterium_formicicum.DSM1535.cds.all.fa
awk '{print $1}' ann_paper/Methanosaeta_concilii_gp6.ASM20441v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Methanosaeta_concilii_gp6.ASM20441v1.cds.all.fa
awk '{print $1}' ann_paper/Methanospirillum_hungatei_jf_1.ASM1344v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Methanospirillum_hungatei_jf_1.ASM1344v1.cds.all.fa
awk '{print $1}' ann_paper/Pelobacter_propionicus_dsm_2379.ASM1504v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Pelobacter_propionicus_dsm_2379.ASM1504v1.cds.all.fa
awk '{print $1}' ann_paper/Syntrophobacter_fumaroxidans_mpob.ASM1496v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Syntrophobacter_fumaroxidans_mpob.ASM1496v1.cds.all.fa
awk '{print $1}' ann_paper/Syntrophomonas_wolfei_subsp_wolfei_str_goettingen_g311.ASM1472v1.cds.all.fa > tmp.fa && mv tmp.fa ann_paper/Syntrophomonas_wolfei_subsp_wolfei_str_goettingen_g311.ASM1472v1.cds.all.fa
cat ann_paper/*cds.all.fa >> ann_paper/transcriptomes.fasta
mkdir ann_paper/simulated ann_paper/simulated/rna 
grep '>' ann_paper/transcriptomes.fasta | awk '{print substr($0, 2)}' > ann_paper
```
Then, go to UniProt's ID mapping [web service](https://www.uniprot.org/uploadlists/), pick the ```ids_to_map.txt``` file, select ```Ensembl Genomes Protein``` in the "From" field, and download results as "Tab-separated". Must contain the following columns: ```yourlist:...```, ```Pathway``` and ```Protein names```. ID map results saved as ```ann_paper/simulated/uniprotinfo.tsv```
```
python annotation_paper/scripts/rnaseq_simulation.py
```

## RNA-Seq quantification

Quantification followed preprocessing with MOSCA and read alignment with Bowtie2.
```
python annotation_paper/scripts/rnaseq_quantification.py
```

## Results analysis

```python
```

## Publication

This paper is still not published.