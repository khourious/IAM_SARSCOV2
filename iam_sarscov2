#!/bin/bash

start=$(date +%s.%N)

[ ! -d $HOME/IAM_SARSCOV2/analysis ] && mkdir $HOME/IAM_SARSCOV2/analysis -v

for i in $(find ./ -type f -name "*.fastq.gz"); do cp "$i" $HOME/IAM_SARSCOV2/analysis -v; done

cd $HOME/IAM_SARSCOV2/analysis

for i in $(find ./ -type f -name "*R1*" -exec basename {} \;); do mv "$i" "$(echo "$i" | rev | cut -c 25- | rev)_R1.fastq.gz" -v; done

for i in $(find ./ -type f -name "*R2*" -exec basename {} \;); do mv "$i" "$(echo "$i" | rev | cut -c 25- | rev)_R2.fastq.gz" -v; done

source activate iam_sarscov2

for i in $(find ./ -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do $HOME/IAM_SARSCOV2/sars2_assembly.sh $HOME/IAM_SARSCOV2/nCoV-2019.reference.fasta "$i"_R1.fastq.gz "$i"_R2.fastq.gz "$i" 12 10 75 $HOME/IAM_SARSCOV2/nCoV-2019.primers.fasta; done

end=$(date +%s.%N)

runtime=$(python -c "print(${end} - ${start})")

echo "" && echo "Done. The runtime was $runtime seconds."
