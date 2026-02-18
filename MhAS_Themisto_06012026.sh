#! /bin/bash
#-----------------------------------------
#SBATCH --job-name=MhAS_ThemistoPseudoalignment

#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --mem=360GB
#SBATCH --time=24:00:00  
#SBATCH --output=MhAS_Themisto_%j.out
#SBATCH --mail-user=nhacker@tamu.edu
#SBATCH --mail-type=ALL

#-----------------------------------------
# Load necessary modules and set environment variables


module load GCCcore/12.2.0

mkdir -p /scratch/user/nhacker/Mh_Adaptive_Sampling/Themisto/Themisto_output

for fq in /scratch/user/nhacker/Mh_Adaptive_Sampling/Centrifuge/Classification_Results/centrifuge_fastq/*.fastq; do
    echo "Processing file: ${fq##*/}"
    fq_basename=${fq##*/} # Remove path
    fq_brief=${fq_basename%%.fastq*} # Remove .fastq extension
    # Run Themisto pseudoalignment on the fastq file
    /scratch/user/nhacker/Mh_Isolate_Paper/Themisto/themisto/build/bin/themisto pseudoalign -q $fq -i /scratch/user/nhacker/Mh_Isolate_Paper/Themisto/Themisto_Mh_Indices/2025_themisto_index_no --temp-dir /scratch/user/nhacker/Mh_Isolate_Paper/Themisto/temp -t 48 --threshold 0.95> /scratch/user/nhacker/Mh_Adaptive_Sampling/Themisto/Themisto_output/${fq_brief}_themisto.txt
    echo "Finished processing file: ${fq_basename}, themisto pseudoalignment output saved to /scratch/user/nhacker/Mh_Adaptive_Sampling/Themisto/Themisto_output/${fq_brief}_themisto.txt"
done