#!/bin/bash

echo -e "Sample_ID\tsR1\tsR2" > 450_samples.tsv

find /prj/DECODE/ea_analysis/data/reads  \( -name "*R1.fq" -o -name "*R2.fq" \) | sort | \
while read -r file; do
    dir=$(dirname "$file")
    sample_id=$(basename "$file" | sed 's/_filtered.*fq//g')

    if [[ $file == *R1.fq ]]; then
        sR1="$file"
    elif [[ $file == *R2.fq ]]; then
        sR2="$file"
        echo -e "$sample_id\t$sR1\t$sR2" >> 450_samples.tsv
    fi
done


# Remove duplicate lines based on sample ID and add sample_ID column
#awk -F'\t' '!seen[$1]++ { OFS="\t"; if ($3 ~ /concatenated_[12]\.fq\.gz$/) print $1, $2, $3; else if (!($1 in samples)) { samples[$1] = 1; print $1, $2, $3 } }' samples.tsv > 450_samples.tsv

