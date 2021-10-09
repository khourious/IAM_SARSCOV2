************
IGM_SARSCOV2
************

-----------------------------
forked from dezordi/ViralFlow
-----------------------------

This repository contains a modification of a set of scripts to performs a reference guided genome assembly of SARS-CoV-2 created by Filipe Dezordi and Gabriel Wallau (FIOCRUZ-IAM).

=======================
Setting up the pipeline
=======================

Download and install the pipeline from the github repo:

.. code:: bash

    git clone --recursive https://github.com/khourious/IGM_SARSCOV2.git; cd IGM_SARSCOV2
    chmod 700 -R DEPENDENCIES
    bash DEPENDENCIES

====================================
How to use the IGM_SARSCOV2 pipeline
====================================

You need to run the ``igm_sarscov2`` and ``igm_sarscov2_summary`` scripts inside the folder that contains the RAW fastq.gz files. For <igm_sarscov2> you need to set the PRIMER SCHEME used in the library preparation (ARTIC_V3, ARTIC_V4 or FIOCRUZ-IOC_V2).

- Examples:

.. code:: bash

    igm_sarscov2 ARTIC_V3

.. code:: bash

    igm_sarscov2 ARTIC_V4

.. code:: bash

    igm_sarscov2 FIOCRUZ-IOC_V2

.. code:: bash

    igm_sarscov2_summary

----------
Files info
----------

.. code-block:: text

    IGM_SARSCOV2/
     ├-DEPENDENCIES            ### script for install and update dependencies
     ├-ARTIC_V3                ### ARTIC V3 primers
     ├-ARTIC_V4                ### ARTIC V4 primers
     ├-FIOCRUZ-IOC_V2          ### FIOCRUZ-IOC V2 primers
     ├-NC_045512.2.fasta       ### SARS-CoV-2 reference sequence
     └-bash_scripts/
      ├-igm_sarscov2           ### perform the genome assembly
      ├-igm_sarscov2_summary   ### do statistics and run pangolin and nextclade (FIOCRUZ-IGM modifications)
      ├-sars2_assembly         ### ViralFlow script (forked from dezordi/ViralFlow)
     └-python_scripts/
      ├-bwa_index.py           ### run bwa index
      ├-bwa_mem.py             ### run bwa mem
      ├-fastp.py               ### run fastp
      ├-get_mvs.py             ### perform intrahost variant analysis with bam-readcount and intrahost.py
      ├-intrahost.py           ### identify genomic positions with multi-allele frequencies
      └-ivar.py                ### run ivar variant and ivar consensus

===============
Explained Usage
===============


.. code-block:: text

    IGM_SARSCOV2/
     └-ANALYSIS/
      ├-"$i"_R1.fastq.gz                          ### copy of RAW R1 fastq.gz file
      ├-"$i"_R2.fastq.gz                          ### copy of RAW R2 fastq.gz file
      └-"$i".results/
       ├-"$i".R1.fq.gz                            ### trimmed R1 fastq.gz file
       ├-"$i".R2.fq.gz                            ### trimmed R2 fastq.gz file
       ├-"$i".coverage.pdf                        ### coverage plot
       ├-"$i".depth10.amb.fa                      ### consensus defined with iVar with ambiguous nucleotideos on positions where major allele frequencies correspond at least 60% of depth
       ├-"$i".depth10.fa                          ### consensus defined with iVar
       ├-"$i".depth10.fa.algn                     ### fasta file with alignment of consensus with reference sequence
       ├-"$i".depth10.fa.algn.minor.fa            ### fasta file with minor consensus genome
       ├-"$i".depth10.fa.bc                       ### bam-readcount output, with all nucleotide frequencies by genomic position
       ├-"$i".depth10.fa.bc.intrahost.short.tsv   ### short tsv file with minor variant informations
       ├-"$i".depth10.fa.bc.intrahost.tsv         ### tsv file with minor variant informations
       ├-"$i".ivar60.qual.txt                     ### txt file with quality control informatoins
       ├-"$i".qual.txt                            ### txt file with quality control informations
       ├-"$i".quality.html                        ### html file with quality control informations
       ├-"$i".sorted.bam                          ### sorted bam file
       ├-"$i".sorted.bam.bai                      ### index of sorted bam file
       ├-"$i".time.txt                            ### time in minutes of each step of analysis
       ├-"$i".tsv                                 ### tsv output from iVar with the frequencies of iSNVs
       └-fastp.json                               ### 
      ├-"$library"_consensus.fa                   ### 
      ├-"$library"_coverage_depth.pdf             ### 
      ├-"$library"_folder_info.txt                ### 
      ├-"$library"_log.txt                        ### 
      ├-"$library"_stats.txt                      ### 
      ├-nextclade_all_YYYY-MM-DD.txt              ### nextclade csv output
      ├-pangolin_all_YYYY-MM-DD.txt               ### pangolin lineages information
      └-pangolin_nextclade_log_YYYY-MM-DD.txt     ### pangolin and nexclade log analysis

==========
Disclaimer
==========
* If you use this workflow for academic purposes, please cite the principal repository and preprint article:
    * https://github.com/dezordi/ViralFlow
    * ViralFlow: an automated workflow for SARS-CoV-2 genome assembly, lineage assignment, mutations and intrahost variants detection. Filipe Zimmer Dezordi, Túlio de Lima Campos, Pedro Miguel Carneiro Jeronimo, Cleber Furtado Aksenen, Suzana Porto Almeida, Gabriel Luz Wallau. medRxiv 2021.10.01.21264424; doi: https://doi.org/10.1101/2021.10.01.21264424
