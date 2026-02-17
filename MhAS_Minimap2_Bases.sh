#!/bin/bash
#SBATCH --job-name=MhAS_base_stats_combined
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=300G
#SBATCH --output=MhAS_base_stats_combined.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhacker@tamu.edu
#SBATCH --account=132741167574

# Load modules
module load GCC/13.2.0
module load SAMtools/1.21

# Directories
as_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling"
mini_dir="${as_dir}/Minimap2"
final_dir="${as_dir}/Combined_Base_Stats"
mkdir -p "$final_dir"

# Output files
combined_mh_csv="${final_dir}/All_Mh_base_stats.csv"
combined_16s_csv="${final_dir}/All_16S_base_stats.csv"

# Write headers once
echo "sample,total_bases,unmapped_bases,mapped_bases,proportion_unmapped_bases,proportion_mapped_bases" > "$combined_mh_csv"
echo "sample,total_bases,unmapped_bases,mapped_bases,proportion_unmapped_bases,proportion_mapped_bases" > "$combined_16s_csv"

echo "Starting combined base-level statistics..."

# Find all BAM files recursively under Minimap2
find "$mini_dir" -type f -name "*.sorted.bam" | while read b; do
    filename=$(basename "$b")

    # Identify which reference (Mh or 16S)
    if [[ "$filename" == *_16S.sorted.bam ]]; then
        prefix="${filename%_16S.sorted.bam}"
        target="16S"
    elif [[ "$filename" == *_Mh.sorted.bam ]]; then
        prefix="${filename%_Mh.sorted.bam}"
        target="Mh"
    else
        echo "Warning: $filename does not match expected naming pattern. Skipping."
        continue
    fi

    echo "Processing $filename..."

    # --- Base-level statistics ---
    total_bases=$(samtools view -F 0x900 "$b" | awk '{sum += length($10)} END {print sum+0}')
    unmapped_bases=$(samtools view -f 4 "$b" | awk '{sum += length($10)} END {print sum+0}')
    mapped_bases=$(samtools view -F 0x904 "$b" | awk '{sum += length($10)} END {print sum+0}')

    if [[ "$total_bases" -eq 0 ]]; then
        proportion_unmapped_bases="NA"
        proportion_mapped_bases="NA"
    else
        proportion_unmapped_bases=$(awk "BEGIN { printf \"%.4f\", $unmapped_bases / $total_bases }")
        proportion_mapped_bases=$(awk "BEGIN { printf \"%.4f\", $mapped_bases / $total_bases }")
    fi

    # Append to combined CSV
    if [[ "$target" == "Mh" ]]; then
        echo "${prefix}_${target},${total_bases},${unmapped_bases},${mapped_bases},${proportion_unmapped_bases},${proportion_mapped_bases}" >> "$combined_mh_csv"
    else
        echo "${prefix}_${target},${total_bases},${unmapped_bases},${mapped_bases},${proportion_unmapped_bases},${proportion_mapped_bases}" >> "$combined_16s_csv"
    fi
done

echo "✅ All base-level stats written to:"
echo "  - $combined_mh_csv"
echo "  - $combined_16s_csv"