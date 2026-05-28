#!/bin/bash
#SBATCH --account=def-freder19
#SBATCH --time=02:00:00
#SBATCH --output=output_breseq.txt
#SBATCH --error=error_breseq.txt
#SBATCH --mem=96G

#this script passes cleaned short read sequences through the breseq pipeline to identify introduced mutations
#sequences are compared against the high quality ancestor assembly built in the trycycler script 
#this script must by run for each strain (could be modified to run in a loop)
#note: this script was initial run on Compute Canada server
#raw short reads can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1356529
#we processed these reads in the previous step script, passing through quality control (01_QC_short_read_sequences.sh)

#load libraries
module load StdEnv/2023
module load breseq/0.38.2
module load bowtie2/2.5.2

#path to the short reads that have been passed through quality control
#will include both forward and reverse reads
#additional path to assembly of ancestor, to map back to
target_evolve1="/path/to/cleaned_short_reads/[STRAIN_ID]_R1_cleaned_001.fastq.gz"
target_evolve2="/path/to/cleaned_short_reads/[STRAIN_ID]_R2_cleaned_001.fastq.gz"
target_ancestor="/path/to/ancestor/assembly/allorhizobium_anc_assembly.fasta"

#pass path to raw sequences and ancestor assembly into breseq 
breseq -j 8 -p -r "$target_ancestor" "$target_evolve1" "$target_evolve2" > log.txt
