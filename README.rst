************
IGM_SARSCOV2
************

This repository contains a forked from **dezordi/ViralFlow** of a set of scripts to performs a reference guided genome assembly of SARS-CoV-2 created by Filipe Dezordi and Gabriel Wallau (FIOCRUZ-IAM).

-----------------------
Setting up the pipeline
-----------------------

Download and install the pipeline from the github repo:

.. code:: bash

    git clone --recursive https://github.com/khourious/IGM_SARSCOV2.git; cd IGM_SARSCOV2
    chmod 700 -R INSTALL
    bash INSTALL

------------------------------------
How to use the IGM_SARSCOV2 pipeline
------------------------------------

.. code:: bash

    igm_sarscov2 ARTIC_V4

* You need to:
    * Run the script **inside** the folder that contains the RAW fastq.gz files.
    * Set the PRIMER SCHEME used in the library preparation (i.e.: ARTIC_V3, ARTIC_V4 or FIOCRUZ-IOC_V2).

----------
Files info
----------

.. code-block:: text

    IGM_SARSCOV2/
     ├── INSTALL                   ### script for install dependencies
     ├── ARTIC_V3                  ### ARTIC V3 primers
     ├── ARTIC_V4                  ### ARTIC V4 primers
     ├── FIOCRUZ-IOC_V2            ### FIOCRUZ-IOC V2 primers
     ├── MN908947.3.fasta          ### SARS-CoV-2 reference sequence
     └── bash_scripts/
      ├── igm_sarscov2             ### script for run the entire analysis
      ├── igm_sarscov2_assembly    ### perform the genome assembly using partial ViralFlow script
      ├── igm_sarscov2_summary     ### do statistics and run pangolin and nextclade
      ├── igm_sarscov2_update      ### script for update dependencies
      ├── sars2_assembly           ### ViralFlow script (forked from dezordi/ViralFlow) - v.0.0.5
     └── python_scripts/
      ├── bwa_index.py             ### run bwa index (forked from dezordi/ViralFlow) - v.0.0.5
      ├── bwa_mem.py               ### run bwa mem (forked from dezordi/ViralFlow) - v.0.0.5
      ├── fastp.py                 ### run fastp (forked from dezordi/ViralFlow) - v.0.0.5
      ├── get_mvs.py               ### perform intrahost variant analysis with bam-readcount and intrahost.py (forked from dezordi/ViralFlow) - v.0.0.5
      ├── intrahost.py             ### identify genomic positions with multi-allele frequencies (forked from dezordi/ViralFlow) - v.0.0.5
      └── ivar.py                  ### run iVar variant and iVar consensus (forked from dezordi/ViralFlow) - v.0.0.5

------------
Results info
------------

.. code-block:: text

    IGM_SARSCOV2/
     └── ANALYSIS/
      ├── SAMPLE_ID_R1.fastq.gz                                                     ### temporary copy of RAW R1 fastq.gz file
      ├── SAMPLE_ID_R2.fastq.gz                                                     ### temporary copy of RAW R2 fastq.gz file
      └── SAMPLE_ID.results/
       ├── SAMPLE_ID.R1.fq.gz                                                       ### trimmed R1 fastq.gz file
       ├── SAMPLE_ID.R2.fq.gz                                                       ### trimmed R2 fastq.gz file
       ├── SAMPLE_ID.coverage.pdf                                                   ### coverage plot
       ├── SAMPLE_ID.depth10.amb.fa                                                 ### consensus defined with iVar with ambiguous nucleotideos on positions where major allele frequencies correspond at least 60% of depth
       ├── SAMPLE_ID.depth10.fa                                                     ### consensus defined with iVar
       ├── SAMPLE_ID.depth10.fa.algn                                                ### alignment of consensus with reference sequence
       ├── SAMPLE_ID.depth10.fa.algn.minor.fa                                       ### minor consensus genome
       ├── SAMPLE_ID.depth10.fa.bc                                                  ### nucleotide frequencies by genomic position
       ├── SAMPLE_ID.depth10.fa.bc.intrahost.short.tsv                              ### summary of minor variant informations
       ├── SAMPLE_ID.depth10.fa.bc.intrahost.tsv                                    ### minor variant informations
       ├── SAMPLE_ID.ivar60.qual.txt                                                ### iVar quality call consensus (frequency threshold: 0.60)
       ├── SAMPLE_ID.lineage_report.csv                                             ### pangolin lineage analysis
       ├── SAMPLE_ID.nextclade.csv                                                  ### nextclade analysis
       ├── SAMPLE_ID.qual.txt                                                       ### iVar quality call consensus
       ├── SAMPLE_ID.quality.html                                                   ### fastp quality control informations
       ├── SAMPLE_ID.sorted.bam                                                     ### sorted bam file
       ├── SAMPLE_ID.sorted.bam.bai                                                 ### index of sorted bam file
       ├── SAMPLE_ID.time.txt                                                       ### time in minutes of each step of analysis
       ├── SAMPLE_ID.tsv                                                            ### iVar with the frequencies of iSNVs
       └── fastp.json                                                          ### metafile of fastp quality control informations
      ├── RAW_FOLDER_NAME_consensus_HOSTNAME_YYYY-MM-DD.fasta                  ### multifasta with major consensus genomes
      ├── RAW_FOLDER_NAME_consensus_with_minor_HOSTNAME_YYYY-MM-DD.fasta       ### multifasta with major and minor consensus genomes
      ├── RAW_FOLDER_NAME_coverage_HOSTNAME_YYYY-MM-DD.pdf                     ### library coverage plot
      ├── RAW_FOLDER_NAME_folder_info_HOSTNAME_YYYY-MM-DD.txt                  ### RAW fastq.gz folder info
      ├── RAW_FOLDER_NAME_log_assembly_PRIMERSCHEME_HOSTNAME_YYYY-MM-DD.txt    ### assembly log analysis
      ├── RAW_FOLDER_NAME_log_summary_HOSTNAME_YYYY-MM-DD.txt                  ### summary log analysis
      ├── RAW_FOLDER_NAME_summary_HOSTNAME_YYYY-MM-DD.txt                      ### summary of statistics, pangolin and nextclade
      └── RAW_FOLDER_NAME_update_HOSTNAME_YYYY-MM-DD.txt                       ### update dependencies log

----------
Disclaimer
----------

* If you use this workflow for academic purposes, please cite the principal repository and preprint article:
    * https://github.com/dezordi/ViralFlow
    * ViralFlow: an automated workflow for SARS-CoV-2 genome assembly, lineage assignment, mutations and intrahost variants detection. Filipe Zimmer Dezordi, Túlio de Lima Campos, Pedro Miguel Carneiro Jeronimo, Cleber Furtado Aksenen, Suzana Porto Almeida, Gabriel Luz Wallau. medRxiv 2021.10.01.21264424; doi: https://doi.org/10.1101/2021.10.01.21264424
