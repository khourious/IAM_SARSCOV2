import argparse, os

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)

#Storing argument on variables
args = parser.parse_args()
prefixout = args.prefixout

try:
    os.system(f"bedtools bamtobed -i {prefixout}.sorted.bam > {prefixout}.sorted.bed \
                    && samtools view {prefixout}.sorted.bam -u | bamdst -p {prefixout}.sorted.bed --cutoffdepth 1000 -o . \
                    && gunzip ./region.tsv.gz \
                    && gunzip ./depth.tsv.gz \
                    && sed -i -e 's/NC_045512\.2/'{prefixout}'/g' -e 's/MN908947\.3/'{prefixout}'/g' -e 's/#Chromosome/Sample_ID/g' -e 's/__/\//g' -e 's/--/|/g' -e 's/Avg depth/Avg_depth/g' -e 's/Cov /Cov_/g' -e 's/%//g' chromosomes.report \
                    && awk '{$1=$1} 1' OFS=, chromosomes.report > chromosomes.report.tmp1 \
                    && awk -F, '{print $1,$3,$7,$9,$10}' chromosomes.report.tmp1 > chromosomes.report.tmp2 \
                    && cat chromosomes.report.tmp2 | tr ' ' '\t' > chromosomes.report")
except:
    print("Error in assembly metrics step")