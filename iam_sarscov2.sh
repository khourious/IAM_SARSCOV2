#!/bin/bash

function background {

   ANALYSISDIR="$HOME/IAM_SARSCOV2/ANALYSIS" # analysis path directory

   DEPTH="10" # minimum depth to mask unanssembled regions

   MINLEN="75" # minimum length to trimm sequences

   RAWDIR=$PWD

   THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}' | sed -n '1p')"

   [ ! -d "$ANALYSISDIR" ] && mkdir "$ANALYSISDIR" -v

   chmod 700 -R "$ANALYSISDIR"

   for i in $(find ./ -type f -name "*R1*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d_ -f1)"_R1.fastq.gz -v; done
   # for i in $(find ./ -type f -name "*R1*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d- -f1)"_R1.fastq.gz -v; done

   for i in $(find ./ -type f -name "*R2*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d_ -f1)"_R2.fastq.gz -v; done
   # for i in $(find ./ -type f -name "*R2*"); do cp "$i" "$ANALYSISDIR"/"$(echo "$i" | awk -F'/' '{print $(NF-1)}' | cut -d- -f1)"_R2.fastq.gz -v; done

   for i in $(find "$ANALYSISDIR" -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do bash $HOME/IAM_SARSCOV2/sars2_assembly.sh $HOME/IAM_SARSCOV2/nCoV-2019.reference.fasta "$ANALYSISDIR"/"$i"_R1.fastq.gz "$ANALYSISDIR"/"$i"_R2.fastq.gz "$i" "$THREADS" "$DEPTH" "$MINLEN" $HOME/IAM_SARSCOV2/nCoV-2019.primers.fasta "$ANALYSISDIR"; done

   for i in $(find "$ANALYSISDIR" -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do cat "$ANALYSISDIR"/*.results/chromosomes.report | (sed -u 1q; sort) | uniq | sed '$d' > "$(pwd | awk -F/ '{print $NF}')_stats.txt"; done

   for i in $(find "$ANALYSISDIR" -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do cat "$ANALYSISDIR"/"$i".results/"$i".depth10.minor.fa >> "$(pwd | awk -F/ '{print $NF}')_consensus.tmp"; done

   for i in $(find "$ANALYSISDIR" -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do cat "$ANALYSISDIR"/"$i".results/"$i".depth10.fa >> "$(pwd | awk -F/ '{print $NF}')_consensus.tmp"; done

   source activate iam_sarscov2

   seqkit grep -vip "MN908947.3_minor" "$(pwd | awk -F/ '{print $NF}')_consensus.tmp" | seqkit sort -n > "$(pwd | awk -F/ '{print $NF}')_consensus.fa"

   source activate plot

   for i in $(find "$ANALYSISDIR" -type f -name "*.fastq.gz" | while read o; do basename $o; done | cut -d_ -f1 | sort |uniq); do fastcov -l "$ANALYSISDIR"/"$i".results/"$i".sorted.bam -o "$ANALYSISDIR"/"$i".results/"$i".coverage.pdf; gs -dSAFER -r3000 -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -sOUTPUTFILE="$(pwd | awk -F/ '{print $NF}')_coverage_depth.pdf" "$ANALYSISDIR"/*.results/*.pdf; done

   mv $(pwd | awk -F/ '{print $NF}')* "$ANALYSISDIR"

   rm -rf "$ANALYSISDIR"/*.fastq.gz "$ANALYSISDIR"/*.tmp

}

export -f background

nohup bash -c background > "$(pwd | awk -F/ '{print $NF}')_log.txt" &
