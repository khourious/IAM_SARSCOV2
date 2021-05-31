#!/bin/bash

DATE="$(date +'%Y-%m-%d')"

[ ! -d "$ANALYSISDIR" ] && mkdir "$ANALYSISDIR" -v

echo "folder@read" | tr '@' '\t' > $HOME/IAM_SARSCOV2/ANALYSIS/"$DATE"_folder_info.txt

find ./ -type f -name "*.fastq.gz" | awk -F'/' '{print $(NF-1),$NF}' | tr '[:blank:]' '\t' | awk '/.fastq.gz/d' >> $HOME/IAM_SARSCOV2/ANALYSIS/"$DATE"_folder_info.txt
