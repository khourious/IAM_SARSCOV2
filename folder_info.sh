#!/bin/bash

ANALYSISDIR="$HOME/IAM_SARSCOV2/ANALYSIS" # analysis path directory

[ ! -d "$ANALYSISDIR" ] && mkdir "$ANALYSISDIR"

chmod 700 -R "$ANALYSISDIR"

echo "folder@read" | tr '@' '\t' > "$ANALYSISDIR"/"$(echo $PWD | awk -F'/' '{print $NF}')_folder_info.txt"

find ./ -type f -name "*.fastq.gz" | awk -F'/' '{print $(NF-1),$NF}' | tr '[:blank:]' '\t' | awk '/.fastq.gz/d' >> "$ANALYSISDIR"/"$(echo $PWD | awk -F'/' '{print $NF}')_folder_info.txt"
