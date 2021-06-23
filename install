#!/bin/bash

chmod 700 -R ../IAM_SARSCOV2

MYSHELL=$(echo $SHELL | awk -F/ '{print $NF}')

if [[ -z "$(which conda)" ]]; then
    cd
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
    echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
    export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
        mamba update -y -n base conda
        mamba create -y -n iam_sarscov2 -c conda-forge -c bioconda -c defaults argparse bedtools bam-readcount biopython bwa fastp ivar mafft numpy pandas samtools==1.10 seqkit
        mamba create -y -n nextclade -c conda-forge -c bioconda -c defaults nodejs nextclade_js
        mamba create -y -n plot -c conda-forge -c bioconda -c defaults pysam numpy pandas seaborn
    fi
else
        mamba update -y -n base conda
        mamba create -y -n iam_sarscov2 -c conda-forge -c bioconda -c defaults argparse bedtools bam-readcount biopython bwa fastp ivar mafft numpy pandas samtools==1.10 seqkit
        mamba create -y -n nextclade -c conda-forge -c bioconda -c defaults nodejs nextclade_js
        mamba create -y -n plot -c conda-forge -c bioconda -c defaults pysam numpy pandas seaborn
fi

git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
mamba env create -f environment.yml
source activate pangolin
python setup.py install

echo 'export PATH=$HOME/IAM_SARSCOV2:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
echo 'export PATH=$HOME/IAM_SARSCOV2/bamdst:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
echo 'export PATH=$HOME/IAM_SARSCOV2/dbxcli:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
echo 'export PATH=$HOME/IAM_SARSCOV2/fastcov:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${MYSHELL}rc
