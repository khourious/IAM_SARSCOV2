#rule fastṕ
    #fastp -i ../$FASTQ1 -I ../$FASTQ2 -o $PREFIXOUT.R1.fq.gz -O $PREFIXOUT.R2.fq.gz --cut_front --cut_tail --qualified_quality_phred 20 -l $MIN_LEN -h $PREFIXOUT.quality.html --thread $THREADS --adapter_fasta ../$ADAPTERS
#rule bwa_index
    #bwa index $FASTA
#rule bwa_mem
    #bwa mem -t $THREADS ../$FASTA $PREFIXOUT.R1.fq.gz $PREFIXOUT.R2.fq.gz | samtools sort -o $PREFIXOUT.sorted.bam
    #samtools index $PREFIXOUT.sorted.bam
#rule ivar_consensus
    #samtools mpileup -d 50000 -A --reference ../$FASTA -Q 30 -q 30 $PREFIXOUT.sorted.bam | ivar consensus -p  $PREFIXOUT -q 30 -t 0 -m $DEPTH -n N
#rule ivar_variants
    #samtools mpileup -aa -A -d 50000 --reference ../$FASTA -Q 30 -q 30 $PREFIXOUT.sorted.bam | ivar variants -p $PREFIXOUT -q 30 -t 0.05
#rule minor_variants
    #bam-readcount -d 50000 -b 30 -q 30 -w 0 -f ../$FASTA $PREFIXOUT.sorted.bam > $PREFIXOUT.depth$DEPTH.fa.bc
    #python ../minor_finder.py -in $PREFIXOUT.depth$DEPTH.fa.bc
    #sed -i -e 's/__/\//g' -e 's/--/|/g' $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv
    #python ../major_minor.py -in $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv
    #if [ `wc -l $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv.fmt | awk '{print $1}'` -ge "2" ];then
        #mafft --thread $THREADS --keeplength --add $PREFIXOUT.depth$DEPTH.fa ../$FASTA > $PREFIXOUT.depth$DEPTH.fa.algn
        #python ../put_minor.py -in $PREFIXOUT.depth$DEPTH.fa.algn -mv $PREFIXOUT.depth$DEPTH.fa.bc.fmt.minors.tsv.fmt
#rule pangolin
    #pangolin $PREFIXOUT.depth$DEPTH.fa -t $THREADS --outfile $PREFIXOUT.depth$DEPTH.fa.pango.csv
#rule nextclade
    #nextclade -i $PREFIXOUT.depth$DEPTH.fa -c $PREFIXOUT.depth$DEPTH.nextclade.csv --jobs $THREADS   
#rule assembly_metrics
    #bedtools bamtobed -i $PREFIXOUT.sorted.bam > $PREFIXOUT.sorted.bed
    #samtools view $PREFIXOUT.sorted.bam -u | bamdst -p $PREFIXOUT.sorted.bed -o .
    #gunzip ./region.tsv.gz
    #gunzip ./depth.tsv.gz