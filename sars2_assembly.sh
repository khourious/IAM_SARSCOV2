#!/bin/bash

#================================================================
# modified by lpmor22; 2021-05-29
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
#%    mafft v7.475 (2020/Nov/23)
#%    bedtools v2.30.0
#%    bamdst 1.0.9
#%
#================================================================
#- IMPLEMENTATION
#-    version        $sars2_assembly 0.0.1
#-    authors        Filipe Dezordi and Gabriel Wallau
#-    maintainer     Filipe Dezordi (zimmer.filipe@gmail.com)
#-    username       dezordi
#-    license        GPL
#-    information    dezordi.github.io
#================================================================
#  HISTORY
#     2021/04/17 : dezordi : Script creation
#
#================================================================
# END_OF_HEADER
#================================================================

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

#================================================================
# WORKFLOW
#================================================================

source activate iam_sarscov2

#CREATING INDEX OF REFERENCE GENOME
bwa index $FASTA

#CREATING DIRECTORY TO STORE RESULTS
mkdir $PREFIXOUT.results/
cd $PREFIXOUT.results

#IF THE USER WANT TO RE-ASSEMBLY THE GENOME WITH A DIFFERENT DEPTH THRESHOLD, ONLY THE CONSENSUS GENERATION RUN WILL BE PERFORMED
if [ -f "$PREFIXOUT.sorted.bam" ]; then
    samtools mpileup -a -d 50000 --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar variants -p $PREFIXOUT -q 30 -t 0.05
    samtools mpileup -d 50000 -a --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar consensus -p  $PREFIXOUT -q 30 -t 0 -m $DEPTH -n N
    samtools mpileup -d 50000 -a --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar consensus -p  $PREFIXOUT.ivar060 -q 30 -t 0.60 -m $DEPTH -n N
    mv $PREFIXOUT.fa $PREFIXOUT.depth$DEPTH.fa
    sed -i -e 's/>.*/>'$PREFIXOUT'/g' $PREFIXOUT.depth$DEPTH.fa
    sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.fa
    mv $PREFIXOUT.ivar060.fa $PREFIXOUT.depth$DEPTH.amb.fa
    sed -i -e 's/>.*/>'$PREFIXOUT'/g' $PREFIXOUT.depth$DEPTH.amb.fa
    sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.amb.fa
else
    #QUALITY CHECK
    echo "FASTP:" > $PREFIXOUT.time.txt
    start=$(date +%s%3N)
    fastp -i ../$FASTQ1 -I ../$FASTQ2 -o $PREFIXOUT.R1.fq.gz -O $PREFIXOUT.R2.fq.gz --cut_front --cut_tail --qualified_quality_phred 20 -l $MIN_LEN -h $PREFIXOUT.quality.html --thread $THREADS --adapter_fasta $ADAPTERS
    end=$(date +%s%3N)
    analysis_in_miliseconds=$(expr $end - $start)
    analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
    echo $analysis_in_minutes >> $PREFIXOUT.time.txt
    #MAPPING
    echo "BWA and ivar:" >> $PREFIXOUT.time.txt
    start=$(date +%s%3N)
    bwa mem -t $THREADS $FASTA $PREFIXOUT.R1.fq.gz $PREFIXOUT.R2.fq.gz | samtools sort -o $PREFIXOUT.sorted.bam
    samtools index $PREFIXOUT.sorted.bam
    #GENERATING CONSENSUS WITH MAJOR ALLELE FREQUENCIES
    samtools mpileup -a -d 50000 --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar variants -p $PREFIXOUT -q 30 -t 0.05
    samtools mpileup -d 50000 -a --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar consensus -p  $PREFIXOUT -q 30 -t 0 -m $DEPTH -n N
    samtools mpileup -d 50000 -a --reference $FASTA $PREFIXOUT.sorted.bam  -Q 30 | ivar consensus -p  $PREFIXOUT.ivar060 -q 30 -t 0.60 -m $DEPTH -n N
    mv $PREFIXOUT.fa $PREFIXOUT.depth$DEPTH.fa
    sed -i -e 's/>.*/>'$PREFIXOUT'/g' $PREFIXOUT.depth$DEPTH.fa
    sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.fa
    mv $PREFIXOUT.ivar060.fa $PREFIXOUT.depth$DEPTH.amb.fa
    sed -i -e 's/>.*/>'$PREFIXOUT'/g' $PREFIXOUT.depth$DEPTH.amb.fa
    sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.amb.fa
    echo $analysis_in_minutes >> $PREFIXOUT.time.txt
    #GET PUTATIVE MINOR VARIANTS STEP
    echo "Minor Variant Analysis:" >> $PREFIXOUT.time.txt
    bam-readcount -d 50000 -b 30 -q 30 -w 0 -f $FASTA $PREFIXOUT.sorted.bam > $PREFIXOUT.depth$DEPTH.fa.bc
    python $HOME/IAM_SARSCOV2/minor_finder.py -in $PREFIXOUT.depth$DEPTH.fa.bc
    sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv
    python3 $HOME/IAM_SARSCOV2/major_minor.py -in $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv
    #IF THE LIBRARY CONTAIN MINOR VARIANTS
    if [ `wc -l $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv.fmt | awk '{print $1}'` -ge "2" ];then
        mafft --quiet --thread $THREADS --keeplength --add $PREFIXOUT.depth$DEPTH.fa $FASTA > $PREFIXOUT.depth$DEPTH.fa.algn
        python $HOME/IAM_SARSCOV2/put_minor.py -in $PREFIXOUT.depth$DEPTH.fa.algn -mv $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv.fmt
        mv $PREFIXOUT.depth$DEPTH.fa.algn.minor.fa $PREFIXOUT.depth$DEPTH.minor.fa
        cat $PREFIXOUT.depth$DEPTH.fa $PREFIXOUT.depth$DEPTH.minor.fa > $PREFIXOUT.depth$DEPTH.all.fa
        # nextclade -i $PREFIXOUT.depth$DEPTH.all.fa -c $PREFIXOUT.depth$DEPTH.all.fa.nextclade.csv --jobs $THREADS
        # pangolin $PREFIXOUT.depth$DEPTH.all.fa --outfile $PREFIXOUT.depth$DEPTH.all.fa.pango.csv
    else
        echo
		# nextclade -i $PREFIXOUT.depth$DEPTH.fa -c $PREFIXOUT.depth$DEPTH.nextclade.csv --jobs $THREADS
        # pangolin $PREFIXOUT.depth$DEPTH.fa --outfile $PREFIXOUT.depth$DEPTH.fa.pango.csv
    fi
    #GET ASSEMBLY METRICS
    bedtools bamtobed -i $PREFIXOUT.sorted.bam > $PREFIXOUT.sorted.bed
    samtools view $PREFIXOUT.sorted.bam -u | bamdst -p $PREFIXOUT.sorted.bed -o .
    gunzip region.tsv.gz
    gunzip depth.tsv.gz
    sed -i -e 's/MN908947\.3/'$PREFIXOUT'/g' chromosomes.report
    sed -i -e 's/#Chromosome/Sample ID/g' chromosomes.report
    sed -i -e 's/__/\//g' -e 's/--/|/g' chromosomes.report
    end=$(date +%s%3N)
    analysis_in_miliseconds=$(expr $end - $start)
    analysis_in_minutes="$(($analysis_in_miliseconds / 60000)).$(($analysis_in_miliseconds % 60000))"
    echo $analysis_in_minutes >> $PREFIXOUT.time.txt
fi
cd ..
