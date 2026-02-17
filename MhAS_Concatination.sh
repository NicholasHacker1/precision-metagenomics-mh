#!/bin/bash

module load GCCcore/13.2.0

outdir="/path/to/output"
mkdir -p "$outdir"

cd /path/to/input
for i in *Group/Pool*; do
    echo "Processing $i"
    for barcode_dir in "$i"/fastq_pass/barcode*; do
        [ -d "$barcode_dir" ] || continue

        barcode=$(basename "$barcode_dir")
        group_name=$(basename "$(dirname "$i")")
        pool_name=$(basename "$i")
        full_outdir="${outdir}/${group_name}/${pool_name}/${barcode}"
        mkdir -p "$full_outdir"

        echo "Processing ${group_name}/${pool_name}/${barcode}..."

        # Unzip all .gz files if needed
        gunzip "$barcode_dir"/*fastq.gz 2>/dev/null

        for hour in {6..72..6}; do
            output="0-${hour}_${barcode}.fastq"
            outpath="${full_outdir}/${output}"

            files=""
            for h in $(seq 0 $hour); do
                matches=("$barcode_dir"/*_"$h".fastq)
                if [ -e "${matches[0]}" ]; then
                    files="$files ${matches[@]}"
                else
                    echo "Warning: no match for *_${h}.fastq in $barcode_dir"
                fi
            done

            if [ -n "$files" ]; then
                cat $files > "$outpath"
                echo "${group_name}/${pool_name}_${output} created"
            else
                echo "No files found for 0 to $hour in $barcode_dir - skipping"
            fi
        done
    done
done
