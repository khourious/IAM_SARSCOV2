#!/bin/bash

DATE="$(date +'%Y-%m-%d')"

THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}')"

rm -rf *_all_*

source activate nextclade

mamba update -y -n nextclade --all

nextclade -i *consensus.fa -j "$THREADS" -c nextclade.tmp

cat nextclade.tmp | (sed -u 1q; sort) | sed -e 's/\;/\t/g' > nextclade_all_"$DATE".txt

source activate pangolin

pangolin --update

pangolin *consensus.fa -t "$THREADS" --outfile pangolin.tmp

cat pangolin.tmp | (sed -u 1q; sort) | sed -e 's/\,/\t/g' > pangolin_all_"$DATE".txt

rm -rf *tmp
