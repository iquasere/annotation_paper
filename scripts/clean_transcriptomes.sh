#!/bin/bash

for file in ann_paper/*.cds.all.fa; do
    awk '{print $1}' "$file" > tmp.fa && mv tmp.fa "$file"
done