#!/bin/bash

#================================================================
# THIS SCRIPTS WAS MODIFIED BY <LPMOR22> TO WORK LOCALLY
#================================================================

#================================================================
# HEADER
#================================================================
#% USAGE
#+    $bash sars2_assembly.sh <REFERENCEGENOME> <001.fastq.gz> <002.fastq.gz> <PREFIX> <NUM_THREADS> <DEPTH> <MIN_LEN> <ADAPTERS>
#%
#% DESCRIPTION
#%    This script performs a reference guided genome assembly of SARS-CoV-2. Python scripts were developed based on the wuhan SARS-CoV-2 reference genome NC_045512.2.
#%    The workflow was developed to work with Illumina paired-end reads. Tests with other technologies should be performed.
#%
#% OPTIONS
#%    <REFERENCEGENOME>    -   Fasta file with reference genome
#%    <001.fastq.gz>       -   Fasqt file with positive sense reads (R1)
#%    <002.fastq.gz>       -   Fastq file with negative sense reads (R2)
#%    <PREFIX>             -   Prefix string to store results and to rename consensus genome
#%    <NUM_THREADS>        -   Number of threads
#%    <DEPTH>              -   Minimum depth to mask unanssembled regions
#%    <MIN_LEN>            -   Minimum length to trimm sequences
#%    <ADAPTERS>           -   Fasta file with adapters used in the sequencing analysis
#%
#% EXAMPLES
#%    $bash sars2_assembly.sh reference.fasta code_R1.fastq.gz code_R2.fastq.gz prefix_name 8 5 75 adapters.fa
#%
#% DEPENDENCIES
#%    BWA Version: 0.7.17-r1188
#%    samtools 1.10 Using htslib 1.10.2
#%    fastp 0.20.1
#%    iVar version 1.3.1
#%    bam-readcount version: 0.8.0-unstable-7-625eea2
#%    Python 3.6.10
#%    mafft v7.480 (2021/May/21)
#%    bedtools v2.30.0
#%    bamdst 1.0.9
#%
#================================================================
#- IMPLEMENTATION
#-    version         $sars2_assembly 0.0.2
#-    authors         Filipe Dezordi and Gabriel Wallau
#-    maintainer      Filipe Dezordi (zimmer.filipe@gmail.com)
#-    username        dezordi
#-    license         GPL
#-    information     dezordi.github.io
#================================================================
#  HISTORY
#     2021/04/17 : dezordi : Script creation
#     2021/06/22 : dezordi : Update, improve indels recognition into intrahost variant analysis
#
#================================================================
# ARGUMENTS
#================================================================

FASTA=$1 #reference genome
FASTQ1=$2 #foward reads
FASTQ2=$3 #reverse reads
PREFIXOUT=$4 #prefix for output
THREADS=$5 #number of threads
DEPTH=$6 #minum depth to mask regions
MIN_LEN=$7 #minimum length to trimm reads
ADAPTERS=$8 #fasta file with adapters
ANALYSISDIR=$9 #analysis directory

#================================================================
# WORKFLOW
#================================================================

#CREATING INDEX OF REFERENCE GENOME
python bwa_index.py -in $FASTA

#CREATING DIRECTORY TO STORE RESULTS
cd "$ANALYSISDIR"
mkdir $PREFIXOUT.results/
cd $PREFIXOUT.results

#QUALITY CHECK
echo "FASTP:" > $PREFIXOUT.time.txt
start=$(date +%s%3N)
python fastp.py -f1 $FASTQ1 -f2 $FASTQ2 -pr $PREFIXOUT -mi $MIN_LEN -p $THREADS -a $ADAPTERS
end=$(date +%s%3N)
analysis_in_miliseconds=$(expr $end - $start)
analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
echo $analysis_in_minutes >> $PREFIXOUT.time.txt
#MAPPING
echo "BWA and ivar:" >> $PREFIXOUT.time.txt
start=$(date +%s%3N)
python bwa_mem.py -f $FASTA -pr $PREFIXOUT -p $THREADS
python ivar.py -f $FASTA -pr $PREFIXOUT -dp $DEPTH
end=$(date +%s%3N)
analysis_in_miliseconds=$(expr $end - $start)
analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
echo $analysis_in_minutes >> $PREFIXOUT.time.txt
#GET PUTATIVE MINOR VARIANTS STEP
echo "Minor Variant Analysis:" >> $PREFIXOUT.time.txt
start=$(date +%s%3N)
python get_mvs.py -f $FASTA -pr $PREFIXOUT -dp $DEPTH -p $THREADS
end=$(date +%s%3N)
analysis_in_miliseconds=$(expr $end - $start)
analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
echo $analysis_in_minutes >> $PREFIXOUT.time.txt
#PANGOLIN AND NEXTCLADE
#echo "Pangolin and Nextclade Analysis:" >> $PREFIXOUT.time.txt
#start=$(date +%s%3N)
#python pango_nextclade.py -pr $PREFIXOUT -dp $DEPTH -p $THREADS
#end=$(date +%s%3N)
#analysis_in_miliseconds=$(expr $end - $start)
#analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
#echo $analysis_in_minutes >> $PREFIXOUT.time.txt
#GET ASSEMBLY METRICS
assembly_metrics.py -pr $PREFIXOUT
cd ..
