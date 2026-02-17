#!/bin/bash
#SBATCH --job-name=MhAS_Filtering
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=360G
#SBATCH --output=MhAS_Filtering.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhacker@tamu.edu
#SBATCH --account=132741167574

module load GCC/13.2.0
module load seqtk/1.4
shopt -s nullglob

cat_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling/Cat_reads"
as_dir="/scratch/user/nhacker/Mh_Adaptive_Sampling"

for group in "$cat_dir"/*; do
    for pool in "$group"/Pool*; do
        echo "Extracting ReadID's from ${pool##*/} in ${group##*/}"
        csv_path=("$as_dir/${group##*/}/${pool##*/}/adaptive_sampling/"AS_decisions*.csv)
        id_output="${pool##*/}_sequence_ids.txt"
        action_output="${pool##*/}_sequence_actions.txt"

        awk -F',' -v rid_out="$id_output" -v act_out="$action_output" '
            BEGIN { OFS = "," }
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "read_id") rid_col = i
                    if ($i == "action") act_col = i
                }
                next
            }
            tolower($act_col) == "sequence" {
                print $rid_col > rid_out
                print $act_col > act_out
            }
        ' "$csv_path"

        sort "$action_output" | uniq

        for barcode in "$pool"/barcode*; do
            echo "Filtering files in ${barcode##*/} of ${pool##*/} in ${group##*/}"
            output_dir="$as_dir/Filtered_reads/${group##*/}_${pool##*/}_${barcode##*/}_sequenced"

            mkdir -p "$output_dir"

            for fq in "$barcode"/*.fastq; do
                seqtk subseq "$fq" "$id_output" > "$output_dir/filtered_${pool##*/}_${fq##*/}" &
            done
        done
        wait
    done
done