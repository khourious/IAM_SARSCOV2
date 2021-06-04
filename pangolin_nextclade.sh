#!/bin/bash

function background

{
    DATE="$(date +'%Y-%m-%d')"

    THREADS="$(lscpu | grep 'CPU(s):' | awk '{print $2}')"

    cat *consensus.fa > pangolin_nextclade.tmp

    source activate pangolin

    pangolin --update

    pangolin pangolin_nextclade.tmp -t "$THREADS" --outfile pangolin.tmp

    cat pangolin.tmp | (sed -u 1q; sort) | sed -e 's/\,/\t/g' > pangolin_all_"$DATE".txt

    source activate nextclade

    mamba update -y -n nextclade --all

    nextclade -i pangolin_nextclade.tmp -j "$THREADS" -c nextclade.tmp

    cat nextclade.tmp | (sed -u 1q; sort) | sed -e 's/\;/\t/g' > nextclade_all_"$DATE".txt

    rm -rf *tmp
}

export -f background

nohup bash -c background > pangolin_nextclade_log_"$(date +'%Y-%m-%d')".txt &
