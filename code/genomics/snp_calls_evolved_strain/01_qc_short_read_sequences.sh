#!/bin/bash
#code passes short read sequences through a quality control step using fastp with standard settings
#note: this script was initial run on Compute Canada server
#data are read in from a folder containing all the raw short reads, and loops through
#data are output to a new folder solely for the cleaned reads
#data are available at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1356529


module load StdEnv/2023
module load fastp/0.24.0


for i in /path/to/raw_short_read/*R1*; do 
  R1=${i#/path/to/raw_short_read/}
  R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz} 
  MERGE=${R1//R1_001.fastq.gz/merged_001.fastq.gz}
  R1_P=${R1//001.fastq.gz/cleaned_001.fastq.gz} 
  R2_P=${R2//001.fastq.gz/cleaned_001.fastq.gz}  
  sample=${R1%_*_L002_*gz}

  echo $R1


  fastp  -i /path/to/raw_short_read/"$R1" \
         -o /output/path/to/folder/for/cleaned_reads/"$R1_P" \
         --qualified_quality_phred 15 --unqualified_percent_limit 40  \
         --cut_front --cut_front_window_size 1 \
         --cut_front_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --correction \
         --trim_poly_g --poly_g_min_len 10 --trim_poly_x --poly_x_min_len 10 \
         --thread 8 --verbose --adapter_fasta adapter.fasta

  fastp  -i /path/to/raw_short_read/"$R2" \
         -o /output/path/to/folder/for/cleaned_reads/""$R2_P" \
         --qualified_quality_phred 15 --unqualified_percent_limit 40  \
         --cut_front --cut_front_window_size 1 \
         --cut_front_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --correction \
         --trim_poly_g --poly_g_min_len 10 --trim_poly_x --poly_x_min_len 10 \
         --thread 8 --verbose --adapter_fasta adapter.fasta

done
