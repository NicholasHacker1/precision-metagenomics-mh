#!/bin/bash
#SBATCH --job-name=MhAS_minimap
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=360G
#SBATCH --output=MhAS_minimap.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhacker@tamu.edu
#SBATCH --account=132741167574

# Load Modules that will be used in code
module load GCC/13.2.0
module load minimap2/2.28
module load SAMtools/1.21

as_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling"
mini_dir="${as_dir}/Minimap2"
filter_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling/Filtered_reads"
ref_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling/MhAS_Ref"
ref_16S="$ref_dir/16S_combined.fasta"
ref_Mh="$ref_dir/AS_3G_Reference.fasta"

mkdir -p "$mini_dir"

# Mapping loop
for seq_dir in "$filter_dir"/*; do
    for fq in "$seq_dir"/*.fastq; do
        map_out="$mini_dir/${seq_dir##*/}"
        mkdir -p "$map_out"
        echo "Aligning ${fq##*/}"
        minimap2 -t 48 -ax map-ont "$ref_Mh" "$fq" | \
            samtools view -bS - | \
            samtools sort -o "$map_out/${fq##*/}_Mh.sorted.bam"
        minimap2 -t 48 -ax map-ont "$ref_16S" "$fq" | \
            samtools view -bS - | \
            samtools sort -o "$map_out/${fq##*/}_16S.sorted.bam"
    done
done
echo "Mapping completed."

stat_dir="${as_dir}/Stats"
fastq_dir="${mini_dir}/Fastq"
mkdir -p "$stat_dir"
mkdir -p "$fastq_dir"

# Stats and FASTQ extraction loop
for m in "$mini_dir"/*; do
    # Create CSV files with headers per directory (sample group)
    Mh_csv="${stat_dir}/${m##*/}_Mh_stats.csv"
    _16S_csv="${stat_dir}/${m##*/}_16S_stats.csv"

    echo "sample,total_reads,unmapped_reads,mapped_reads,proportion_unmapped,proportion_mapped" > "$Mh_csv"
    echo "sample,total_reads,unmapped_reads,mapped_reads,proportion_unmapped,proportion_mapped" > "$_16S_csv"

    for b in "$m"/*.sorted.bam; do
        echo "Processing ${b##*/}"
        filename=$(basename "$b")

        if [[ "$filename" == *_16S.sorted.bam ]]; then
            prefix="${filename%_16S.sorted.bam}"
            target="16S"
        elif [[ "$filename" == *_Mh.sorted.bam ]]; then
            prefix="${filename%_Mh.sorted.bam}"
            target="Mh"
        else
            echo "Warning: BAM file $filename does not match expected naming. Skipping."
            continue
        fi

        # Set output FASTQ filenames
        unmapped_fastq="${fastq_dir}/${prefix}_${target}_unmapped.fastq"
        mapped_fastq="${fastq_dir}/${prefix}_${target}_mapped.fastq"

        # Count total reads (excluding secondary and supplementary)
        total_reads=$(samtools view -F 0x900 "$b" | cut -f 1 | sort -u | wc -l)

        # Count unmapped reads
        unmapped_reads=$(samtools view -f 4 "$b" | cut -f 1 | sort -u | wc -l)

        # Count mapped reads (primary only, excluding unmapped and secondary/supplementary)
        mapped_reads=$(samtools view -F 0x904 "$b" | cut -f 1 | sort -u | wc -l)

        # Calculate proportions
        if [[ "$total_reads" -eq 0 ]]; then
            proportion_unmapped="NA"
            proportion_mapped="NA"
        else
            proportion_unmapped=$(awk "BEGIN { printf \"%.4f\", $unmapped_reads / $total_reads }")
            proportion_mapped=$(awk "BEGIN { printf \"%.4f\", $mapped_reads / $total_reads }")
        fi

        # Output fastq files (add logging for debugging)
        echo "Writing unmapped FASTQ to: $unmapped_fastq"
        samtools view -b -f 4 "$b" | samtools fastq - > "$unmapped_fastq"

        echo "Writing mapped FASTQ to: $mapped_fastq"
        samtools view -b -F 0x904 "$b" | samtools fastq - > "$mapped_fastq"

        # Append stats to appropriate CSV
        if [[ "$target" == "Mh" ]]; then
            echo "${prefix}_${target},${total_reads},${unmapped_reads},${mapped_reads},${proportion_unmapped},${proportion_mapped}" >> "$Mh_csv"
        else
            echo "${prefix}_${target},${total_reads},${unmapped_reads},${mapped_reads},${proportion_unmapped},${proportion_mapped}" >> "$_16S_csv"
        fi

    done
done
echo "All tasks completed."

# Concatenate all Mh and 16S stats CSVs into single files with one header each
final_dir="${as_dir}/Combined_Stats"
mkdir -p "$final_dir"

# Final combined CSVs
combined_mh_csv="${final_dir}/All_Mh_stats.csv"
combined_16s_csv="${final_dir}/All_16S_stats.csv"

# Get header from any one CSV (they're all the same)
head -n 1 $(find "$stat_dir" -name "*_Mh_stats.csv" | head -n 1) > "$combined_mh_csv"
head -n 1 $(find "$stat_dir" -name "*_16S_stats.csv" | head -n 1) > "$combined_16s_csv"

# Append all rows except header
find "$stat_dir" -name "*_Mh_stats.csv" | while read file; do
    tail -n +2 "$file" >> "$combined_mh_csv"
done

find "$stat_dir" -name "*_16S_stats.csv" | while read file; do
    tail -n +2 "$file" >> "$combined_16s_csv"
done

echo "Combined stats CSVs created in $final_dir"
