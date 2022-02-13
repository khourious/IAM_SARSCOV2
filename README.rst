************
IGM_SARSCOV2
************

This repository contains a forked from **dezordi/ViralFlow** of a set of scripts to performs a reference guided genome assembly of SARS-CoV-2 created by Filipe Dezordi and Gabriel Wallau (FIOCRUZ-IAM).

.. code-block:: text

    IGM_SARSCOV2
    ------------
    Reference guided genome assembly of SARS-CoV-2 (Using ViralFlow v0.0.5 - dezordi/ViralFlow)

    Assembly + summary stats:
      igm_sarscov2 -w 1 -p <primer scheme> -i <input path> -d <depth> -t <threads>

    Assembly:
      igm_sarscov2 -w 2 -p <primer scheme> -i <input path> -d <depth> -t <threads>

    Summary stats:
      igm_sarscov2 -w 3 -p <primer scheme> -i <input path> -t <threads>

    Command line parameters:
      -d arg    Minimum depth to variant calling (default: 10).
      -i arg    Path containing the Illumina fastQ sequencing data.
      -p arg    Primer scheme name used for generate tiling amplicons (i.e. ARTIC_V4.1).
      -t arg    Max number of threads (default: all cores minus 2).
      -u        Update conda dependencies.
      -v        Display version information and check for update.
      -w arg    Set the workflow option (1- assembly+summary; 2- assembly; 3- summary).

-----------------------
Setting up the pipeline
-----------------------

Download and install the pipeline from the github repo:

.. code:: bash

    git clone --recursive https://github.com/khourious/IGM_SARSCOV2.git; cd IGM_SARSCOV2
    chmod 700 -R INSTALL
    bash INSTALL

--------------------------------
How to use the complete workflow
--------------------------------

* For use, requires:
    * Set the workflow option number 1.
    * Primer scheme name used for generate tiling amplicons (i.e.: ARTIC_V4, ARTIC_V4-1 or FIOCRUZ-IOC_V2).
    * Path containing the Illumina fastQ sequencing data.

* Optional:
    * Max number of threads (default: all cores minus 2).
    * Minimum depth to variant calling (default: 10).

.. code:: bash

    igm_sarscov2 -w 1 -p ARTIC_V4-1 -i /home/user/BaseSpaceFolder/LIBRARRY_NAME -d 10 -t 12

--------------------------------
How to use the assembly workflow
--------------------------------

* For use, requires:
    * Set the workflow option number 2.
    * Primer scheme name used for generate tiling amplicons (i.e.: ARTIC_V4, ARTIC_V4-1 or FIOCRUZ-IOC_V2).
    * Path containing the Illumina fastQ sequencing data.

* Optional:
    * Max number of threads (default: all cores minus 2).
    * Minimum depth to variant calling (default: 10).

.. code:: bash

    igm_sarscov2 -w 2 -p ARTIC_V4-1 -i /home/user/BaseSpaceFolder/LIBRARRY_NAME -d 10 -t 12

-------------------------------------
How to use the summary stats workflow
-------------------------------------

* For use, requires:
    * Set the workflow option number 3.
    * Primer scheme name used for generate tiling amplicons (i.e.: ARTIC_V4, ARTIC_V4-1 or FIOCRUZ-IOC_V2).
    * Path containing the Illumina fastQ sequencing data.

* Optional:
    * Max number of threads (default: all cores minus 2).

.. code:: bash

    igm_sarscov2 -w 3 -p ARTIC_V4-1 -i /home/user/BaseSpaceFolder/LIBRARRY_NAME -t 12

----------
Files info
----------

.. code-block:: text

    IGM_SARSCOV2/
     ├── INSTALL                      ### script for install dependencies
     └── primer_schemes/
      ├── ARTIC_V3.fasta              ### ARTIC V3 primers
      ├── ARTIC_V4.fasta              ### ARTIC V4 primers
      ├── ARTIC_V4-1.fasta            ### ARTIC V4 primers
      ├── FIOCRUZ-IOC_V2.fasta        ### FIOCRUZ-IOC V2 primers
     └── ref_seq/
     ├── MN908947.3.fasta             ### SARS-CoV-2 reference sequence
     └── scripts/
      ├── bwa_index.py                ### run bwa index (forked from dezordi/ViralFlow) - v.0.0.5
      ├── bwa_mem.py                  ### run bwa mem (forked from dezordi/ViralFlow) - v.0.0.5
      ├── fastp.py                    ### run fastp (forked from dezordi/ViralFlow) - v.0.0.5
      ├── get_mvs.py                  ### perform intrahost variant analysis with bam-readcount and intrahost.py (forked from dezordi/ViralFlow) - v.0.0.5
      ├── igm_sarscov2                ### script for run the analysis
      ├── intrahost.py                ### identify genomic positions with multi-allele frequencies (forked from dezordi/ViralFlow) - v.0.0.5
      ├── ivar.py                     ### run iVar variant and iVar consensus (forked from dezordi/ViralFlow) - v.0.0.5
      └── sars2_assembly              ### ViralFlow script (forked from dezordi/ViralFlow) - v.0.0.5

------------
Results info
------------

.. code-block:: text

    IGM_SARSCOV2/
     └── ANALYSIS/
      ├── SAMPLE.R1.fastq.gz                                                                     ### temporary copy of RAW R1 fastq.gz file
      ├── SAMPLE.R2.fastq.gz                                                                     ### temporary copy of RAW R2 fastq.gz file
      └── SAMPLE.results/
       ├── SAMPLE.R1.fq.gz                                                                       ### trimmed R1 fastq.gz file
       ├── SAMPLE.R2.fq.gz                                                                       ### trimmed R2 fastq.gz file
       ├── SAMPLE.coverage.pdf                                                                   ### coverage plot
       ├── SAMPLE.depthXX.amb.fa                                                                 ### consensus defined with iVar with ambiguous nucleotideos on positions where major allele frequencies correspond at least 60% of depth
       ├── SAMPLE.depthXX.fa                                                                     ### consensus defined with iVar
       ├── SAMPLE.depthXX.fa.algn                                                                ### alignment of consensus with reference sequence
       ├── SAMPLE.depthXX.fa.algn.minor.fa                                                       ### minor consensus genome
       ├── SAMPLE.depthXX.fa.bc                                                                  ### nucleotide frequencies by genomic position
       ├── SAMPLE.depthXX.fa.bc.intrahost.short.tsv                                              ### summary of minor variant informations
       ├── SAMPLE.depthXX.fa.bc.intrahost.tsv                                                    ### minor variant informations
       ├── SAMPLE.ivar60.qual.txt                                                                ### iVar quality call consensus (frequency threshold: 0.60)
       ├── SAMPLE.lineage_report.csv                                                             ### pangolin lineage analysis
       ├── SAMPLE.nextclade.csv                                                                  ### nextclade analysis
       ├── SAMPLE.qual.txt                                                                       ### iVar quality call consensus
       ├── SAMPLE.quality.html                                                                   ### fastp quality control informations
       ├── SAMPLE.sorted.bam                                                                     ### sorted bam file
       ├── SAMPLE.sorted.bam.bai                                                                 ### index of sorted bam file
       ├── SAMPLE.time.txt                                                                       ### time in minutes of each step of analysis
       ├── SAMPLE.tsv                                                                            ### iVar with the frequencies of iSNVs
       └── fastp.json                                                                            ### metafile of fastp quality control informations
      ├── LIBRARYNAME.folder_info.HOSTNAME.YYYY-MM-DD.txt                                        ### RAW fastq.gz folder info
      ├── LIBRARYNAME.PRIMERSCHEME.depthXX.consensus.HOSTNAME.YYYY-MM-DD.fasta                   ### multifasta with major consensus genomes
      ├── LIBRARYNAME.PRIMERSCHEME.depthXX.consensus_with_minor.HOSTNAME.YYYY-MM-DD.fasta        ### multifasta with major and minor consensus genomes
      ├── LIBRARYNAME.PRIMERSCHEME.depthXX.coverage.HOSTNAME.YYYY-MM-DD.pdf                      ### library coverage plot
      ├── LIBRARYNAME.PRIMERSCHEME.depthXX.log.HOSTNAME.YYYY-MM-DD.txt                           ### log analysis
      └── LIBRARYNAME.PRIMERSCHEME.depthXX.summary.HOSTNAME.YYYY-MM-DD.txt                       ### summary of statistics, pangolin and nextclade

----------
Disclaimer
----------

* If you use this workflow for academic purposes, please cite the principal repository and preprint article:
    * https://github.com/dezordi/ViralFlow
    * ViralFlow: an automated workflow for SARS-CoV-2 genome assembly, lineage assignment, mutations and intrahost variants detection. Filipe Zimmer Dezordi, Túlio de Lima Campos, Pedro Miguel Carneiro Jeronimo, Cleber Furtado Aksenen, Suzana Porto Almeida, Gabriel Luz Wallau. medRxiv 2021.10.01.21264424; doi: https://doi.org/10.1101/2021.10.01.21264424
