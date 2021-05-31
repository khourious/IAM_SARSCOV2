#!/bin/bash

start=$(date +%s.%N)

DATE="$(date +'%Y-%m-%d')"
ANALYSISDIR="$HOME/IAM_SARSCOV2/ANALYSIS" # analysis path directory
THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"
DEPTH="10" # minimum depth to mask unanssembled regions
MINLEN="75" # minimum length to trimm sequences

[ ! -d "$ANALYSISDIR" ] && mkdir "$ANALYSISDIR" -v

for i in $(find ./ -type f -name "*R1*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d- -f1)"_R1.fastq.gz -v; done
for i in $(find ./ -type f -name "*R2*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d- -f1)"_R2.fastq.gz -v; done
# for i in $(find ./ -type f -name "*R1*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d_ -f1)"_R1.fastq.gz -v; done
# for i in $(find ./ -type f -name "*R2*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d_ -f1)"_R2.fastq.gz -v; done

cd "$ANALYSISDIR"

for i in $(find ./ -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do bash $HOME/IAM_SARSCOV2/sars2_assembly.sh $HOME/IAM_SARSCOV2/nCoV-2019.reference.fasta "$i"_R1.fastq.gz "$i"_R2.fastq.gz "$i" "$THREADS" "$DEPTH" "$MINLEN" $HOME/IAM_SARSCOV2/nCoV-2019.primers.fasta; done

for i in $(find ./ -type f -name "*chromosomes.report"); do cat "$i"; done | (sed -u 1q; sort) | uniq | sed '$d' > "$DATE"_stats.txt
for i in $(find ./ -type f -name "*.depth10.minor.fa"); do cat "$i" >> "$DATE"_consensus.tmp; done
for i in $(find ./ -type f -name "*.depth10.fa"); do cat "$i" >> "$DATE"_consensus.tmp; done
seqkit grep -vip "MN908947.3_minor" "$DATE"_consensus.tmp | seqkit sort -n - > "$DATE"_consensus.fa

source activate plot
for i in $(find ./ -type f -name "*.bam" | while read o; do basename $o | cut -d. -f1; done | sort | uniq); do fastcov -l "$i".results/"$i".sorted.bam -o "$i".results/"$i".coverage.pdf; gs -dSAFER -r3000 -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -sOUTPUTFILE="$DATE"_coverage_depth.pdf "$i".results/*.pdf; done

rm -rf *.fastq.gz *.tmp

end=$(date +%s.%N)

runtime=$(python -c "print(${end} - ${start})")

echo "" && echo "Done. The runtime was $runtime seconds."
