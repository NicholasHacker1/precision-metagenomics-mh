#!/bin/sh

module load Anaconda3
source activate centrifuge_env

index_path="/path/to/index"
out_path="/path/to/output"
centrifuge_output="${out_path}/centrifuge_output"
kreport_output="${out_path}/kreport_output"
input_path="path/to/input"

mkdir -p $centrifuge_output
mkdir -p $kreport_output

for file in "${input_path}"/*Mh_mapped.fastq; do
    [ ! -f "$file" ] && continue
    filename=$(basename "$file" .fastq)
    echo "Processing $filename"
    centrifuge -x "$index_path" -U "$file" -S "$centrifuge_output/${filename}_centrifuge.tsv" 2> "$centrifuge_output/${filename}_centrifuge.log"
    centrifuge-kreport -x "$index_path" "$centrifuge_output/${filename}_centrifuge.tsv" > "$kreport_output/${filename}_kreport.txt" 2> "$kreport_output/${filename}_kreport.log"
    echo "Finished processing $filename"
done
