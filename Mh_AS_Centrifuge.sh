#!/bin/sh

#SBATCH --job-name=Mh_AS_Centrifuge
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=360G
#SBATCH --output=centrifuge_output.%j
#SBATCH --account=132741167574
##SBATCH --mail-type=ALL
##SBATCH --mail-user=nhacker@tamu.edu

module load Anaconda3
source activate centrifuge_env

index_path="/scratch/user/nhacker/Mh_Adaptive_Sampling/Centrifuge/centrifuge_index/centrifuge_index_archaea_bacteria"
out_path="/scratch/user/nhacker/Mh_Adaptive_Sampling/Centrifuge/Classification_Results"
centrifuge_output="${out_path}/centrifuge_output"
kreport_output="${out_path}/kreport_output"
input_path="/scratch/user/nhacker/Mh_Adaptive_Sampling/Minimap2_Fastq/Fastq"

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