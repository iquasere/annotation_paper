#!/bin/bash

for group in "" "/first_group" "/second_group" "/third_group" "/fourth_group" "/fifth_group"
do
  gunzip ann_paper"$group"/*.gz
  cat ann_paper"$group"/*.fna > ann_paper"$group"/genomes.fasta
  awk '{print $1}' ann_paper"$group"/genomes.fasta > ann_paper"$group"/genomes_fixed.fasta
done
