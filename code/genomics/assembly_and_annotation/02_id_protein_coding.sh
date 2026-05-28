#!/bin/bash
#code passes the reference ancestor genome through prokka for annotation
#primary utility here is for an easy passage through prodigal for identifying protein-coding genes 
#note: this script was initially run on Compute Canada server

#load libraries
module load StdEnv/2020
module load gcc/9.3.0
module load prokka/1.14.5

prokka /path/to/ancestor/allorhizobium_anc_assembly.fasta --outdir /prokka/output/predicted/proteincoding/genes
#while there is a lot of output, the salient pieces are: prokka_output.gff & prokka_output.tsv & prokka_output.faa
#these can be found in the data folder in on the git repo 
#PROKKA_05052025.gff
#PROKKA_05052025.tsv
#PROKKA_05052025.faa
