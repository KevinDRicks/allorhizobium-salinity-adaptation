#code takes all the SNPS from the evolved strains and tiles them onto the ancestral genome
#genetic distances between samples can be calculated, allowing us 
#to examine how strains cluster by genetic similarity, and if this is driven by their
#historical selection treatment

#load packages
library("Biostrings")
library('hierfstat')
library('DECIPHER')
library("pwalign")
library('vegan')

##############################################################################################
#read in the snp summary and ids as well as ancestral allorhizobium reference
##############################################################################################
merge_snps <- read.csv("/path/to/snp_summary.csv")
layout <- read.csv('/path/to/sequencing_labels.csv')
working.df <- layout[layout$Strain!='ancestor',]

pl1 <- readDNAStringSet("/path/to/ancestral_plasmid1.fasta")
pl2 <- readDNAStringSet("/path/to/ancestral_plasmid2.fasta")
pl3 <- readDNAStringSet("/path/to/ancestral_plasmid3.fasta")
pl4 <- readDNAStringSet("/path/to/ancestral_plasmid4.fasta")
chr <- readDNAStringSet('/path/to/main_chromosome.fasta')

#recode deletions from variants to a blank space instead of a .
merge_snps[merge_snps$VAR=='.',]$VAR <- ""

###########################################################################################
#we're going to create dna strings of each of the replicons, and add in the observed snps and deletions as necessary
#we break this up by the snps and the indels, as these will be coded in separately
#from these, we can calculate the genetic distance between samples 
###########################################################################################

###########################################################################################
#main chromosome
###########################################################################################
#create dna string
chr.empt.string <- DNAStringSet()
chr.empt.string[c(1:nrow(working.df))] <- chr
names(chr.empt.string) <- working.df$Strain

#subset polymorphisms those on chromosome
working.snp <- merge_snps[merge_snps$Cluster=="cluster_001_consensus",]
working.snp <- working.snp[order(working.snp$Position,decreasing = T),]
working.indel <- working.snp[nchar(working.snp$REF)!=1 | nchar(working.snp$VAR)!=1,]
working.snp <- working.snp[nchar(working.snp$REF)==1 & nchar(working.snp$VAR)==1,]

#introduce snps
for(i in 1:nrow(working.snp)){
  working.st <- chr.empt.string[names(chr.empt.string)==working.snp$Strain[i]]
  sub.seq <- DNAString(working.snp$VAR[i])
  insert_pos <- working.snp$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + 1, nchar(working.st))
  ))
  chr.empt.string[[which(names(chr.empt.string)==working.snp$Strain[i])]] <- new_seq
  print(i)
}

#introduce indels
for(i in 1:nrow(working.indel)){
  working.st <- chr.empt.string[names(chr.empt.string)==working.indel$Strain[i]]
  sub.seq <- DNAString(working.indel$VAR[i])
  insert_pos <- working.indel$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + nchar(working.indel$REF[i]), nchar(working.st))
  ))
  #ATA
  chr.empt.string[[which(names(chr.empt.string)==working.indel$Strain[i])]] <- new_seq
  print(i)
}


###########################################################################################
#plasmid 1
###########################################################################################
#create dna string
pl1.empt.string <- DNAStringSet()
pl1.empt.string[c(1:nrow(working.df))] <- pl1
names(pl1.empt.string) <- working.df$Strain

#subset polymorphisms those on chromosome
working.snp <- merge_snps[merge_snps$Cluster=="cluster_002_consensus",]
working.snp <- working.snp[order(working.snp$Position,decreasing = T),]
working.indel <- working.snp[nchar(working.snp$REF)!=1 | nchar(working.snp$VAR)!=1,]
working.snp <- working.snp[nchar(working.snp$REF)==1 & nchar(working.snp$VAR)==1,]

#alter snps
for(i in 1:nrow(working.snp)){
  working.st <- pl1.empt.string[names(pl1.empt.string)==working.snp$Strain[i]]
  sub.seq <- DNAString(working.snp$VAR[i])
  insert_pos <- working.snp$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + 1, nchar(working.st))
  ))
  pl1.empt.string[[which(names(pl1.empt.string)==working.snp$Strain[i])]] <- new_seq
  print(i)
}

#no indels found on this replicon, so code is commented out
#alter indels
#for(i in 1:nrow(working.indel)){
#  working.st <- pl1.empt.string[names(pl1.empt.string)==working.indel$Strain[i]]
#  sub.seq <- DNAString(working.indel$VAR[i])
#  insert_pos <- working.indel$Position[i]
#  #Insert the new sequence
#  new_seq <- DNAString(paste0(
#    substr(as.character(working.st), 1, insert_pos-1),
#    as.character(sub.seq),
#    substr(as.character(working.st), insert_pos + nchar(working.indel$REF[i]), nchar(working.st))
#  ))
#  #ATA               A
#  pl1.empt.string[[which(names(pl1.empt.string)==working.indel$Strain[i])]] <- new_seq
#  print(i)
#}



###########################################################################################
#plasmid 2
###########################################################################################
#create dna string
pl2.empt.string <- DNAStringSet()
pl2.empt.string[c(1:nrow(working.df))] <- pl2
names(pl2.empt.string) <- working.df$Strain

#subset polymorphisms those on chromosome
working.snp <- merge_snps[merge_snps$Cluster=="cluster_003_consensus",]
working.snp <- working.snp[order(working.snp$Position,decreasing = T),]
working.indel <- working.snp[nchar(working.snp$REF)!=1 | nchar(working.snp$VAR)!=1,]
working.snp <- working.snp[nchar(working.snp$REF)==1 & nchar(working.snp$VAR)==1,]

#alter snps
for(i in 1:nrow(working.snp)){
  working.st <- pl2.empt.string[names(pl2.empt.string)==working.snp$Strain[i]]
  sub.seq <- DNAString(working.snp$VAR[i])
  insert_pos <- working.snp$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + 1, nchar(working.st))
  ))
  pl2.empt.string[[which(names(pl2.empt.string)==working.snp$Strain[i])]] <- new_seq
  print(i)
}

#no indels found on this replicon, so code is commented out
##alter indels
#for(i in 1:nrow(working.indel)){
#  working.st <- pl2.empt.string[names(pl2.empt.string)==working.indel$Strain[i]]
#  sub.seq <- DNAString(working.indel$VAR[i])
#  insert_pos <- working.indel$Position[i]
#  #Insert the new sequence
#  new_seq <- DNAString(paste0(
#    substr(as.character(working.st), 1, insert_pos-1),
#    as.character(sub.seq),
#    substr(as.character(working.st), insert_pos + nchar(working.indel$REF[i]), nchar(working.st))
#  ))
#  #ATA               A
#  pl2.empt.string[[which(names(pl2.empt.string)==working.indel$Strain[i])]] <- new_seq
#  print(i)
#}



###########################################################################################
#plasmid 3
###########################################################################################
#create dna string
pl3.empt.string <- DNAStringSet()
pl3.empt.string[c(1:nrow(working.df))] <- pl3
names(pl3.empt.string) <- working.df$Strain

#subset polymorphisms those on chromosome
working.snp <- merge_snps[merge_snps$Cluster=="cluster_004_consensus",]
working.snp <- working.snp[order(working.snp$Position,decreasing = T),]
working.indel <- working.snp[nchar(working.snp$REF)!=1 | nchar(working.snp$VAR)!=1,]
working.snp <- working.snp[nchar(working.snp$REF)==1 & nchar(working.snp$VAR)==1,]

#alter snps
for(i in 1:nrow(working.snp)){
  working.st <- pl3.empt.string[names(pl3.empt.string)==working.snp$Strain[i]]
  sub.seq <- DNAString(working.snp$VAR[i])
  insert_pos <- working.snp$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + 1, nchar(working.st))
  ))
  pl3.empt.string[[which(names(pl3.empt.string)==working.snp$Strain[i])]] <- new_seq
  print(i)
}

#no indels found on this replicon, so code is commented out
#alter indels
#for(i in 1:nrow(working.indel)){
#  working.st <- pl3.empt.string[names(pl3.empt.string)==working.indel$Strain[i]]
#  sub.seq <- DNAString(working.indel$VAR[i])
#  insert_pos <- working.indel$Position[i]
#  #Insert the new sequence
#  new_seq <- DNAString(paste0(
#    substr(as.character(working.st), 1, insert_pos-1),
#    as.character(sub.seq),
#    substr(as.character(working.st), insert_pos + nchar(working.indel$REF[i]), nchar(working.st))
#  ))
#  #ATA               A
#  pl3.empt.string[[which(names(pl3.empt.string)==working.indel$Strain[i])]] <- new_seq
#  print(i)
#}




###########################################################################################
#plasmid 4
###########################################################################################
#create dna string
pl4.empt.string <- DNAStringSet()
pl4.empt.string[c(1:nrow(working.df))] <- pl4
names(pl4.empt.string) <- working.df$Strain

#subset polymorphisms those on chromosome
working.snp <- merge_snps[merge_snps$Cluster=="cluster_005_consensus",]
working.snp <- working.snp[order(working.snp$Position,decreasing = T),]
working.indel <- working.snp[nchar(working.snp$REF)!=1 | nchar(working.snp$VAR)!=1,]
working.snp <- working.snp[nchar(working.snp$REF)==1 & nchar(working.snp$VAR)==1,]

#alter snps
for(i in 1:nrow(working.snp)){
  working.st <- pl4.empt.string[names(pl4.empt.string)==working.snp$Strain[i]]
  sub.seq <- DNAString(working.snp$VAR[i])
  insert_pos <- working.snp$Position[i]
  #Insert the new sequence
  new_seq <- DNAString(paste0(
    substr(as.character(working.st), 1, insert_pos-1),
    as.character(sub.seq),
    substr(as.character(working.st), insert_pos + 1, nchar(working.st))
  ))
  pl4.empt.string[[which(names(pl4.empt.string)==working.snp$Strain[i])]] <- new_seq
  print(i)
}

#no indels found on this replicon, so code is commented out
#alter  indel
#for(i in 1:nrow(working.indel)){
#  working.st <- pl4.empt.string[names(pl4.empt.string)==working.indel$Strain[i]]
#  sub.seq <- DNAString(working.indel$VAR[i])
#  insert_pos <- working.indel$Position[i]
#  #Insert the new sequence
#  new_seq <- DNAString(paste0(
#  ))
#  #ATA               A
#  pl4.empt.string[[which(names(pl4.empt.string)==working.indel$Strain[i])]] <- new_seq
#  print(i)
#}



###############################################################################################################################3
#a number of these strains have lost plasmids
#While we could remove them completely from these calculations, these have an enormous effect on determining genetic
#consequently, instead just insert a  deletion event on the plasmid, representing effectively 
#a single deletion. Allows us to examine accumulation of differences
#this code code can be exclude to ignore potential plasmid deletion events
###############################################################################################################################3
#these are the strains that have lost their plasmids
pl_b <- c("St3",  "St41", "St42", "St43", "St48", "St49", "St50")
pl_c <- c("St22", "St24", "St33", "St39")

for(i in 1:length(pl_b)){
  temp <- pl2.empt.string[names(pl2.empt.string)==pl_b[i]]
  pl2.empt.string[names(pl2.empt.string)==pl_b[i]] <- DNAStringSet(substring(temp,1,nchar(temp)-2))
}
for(i in 1:length(pl_c)){
  temp <- pl3.empt.string[names(pl3.empt.string)==pl_c[i]]
  pl3.empt.string[names(pl3.empt.string)==pl_c[i]] <- DNAStringSet(substring(temp,1,nchar(temp)-2))
}


###############################################################################################################################
#create alignments for each replicon and merge together to analyze distance
###############################################################################################################################
aligned.chr <- DECIPHER::AlignSeqs(chr.empt.string) 
aligned.pl1 <- DECIPHER::AlignSeqs(pl1.empt.string) 
aligned.pl2 <- DECIPHER::AlignSeqs(pl2.empt.string) 
aligned.pl3 <- DECIPHER::AlignSeqs(pl3.empt.string) 
aligned.pl4 <- DECIPHER::AlignSeqs(pl4.empt.string) 


#merge all aligned replicons together
merged.empt.string <- DNAStringSet()
for(i in 1:nrow(working.df)){
  merged.empt.string <- c(merged.empt.string,
                          DNAStringSet(paste(as.character(aligned.chr[i]),
                                             as.character(aligned.pl1[i]),
                                             as.character(aligned.pl2[i]),
                                             as.character(aligned.pl3[i]),
                                             as.character(aligned.pl4[i]),
                                             sep='')))
  print(i)
}


#calculate distance between all strains
hamm.dist.all <- pwalign::stringDist(merged.empt.string,method='hamming')

adonis2(hamm.dist.all~Nitrogen*Salinity,working.df,by='terms')


rep.mds.all <-  metaMDS(hamm.dist.all,k = 5,  maxit = 999, trymax = 100)
working.df$Merge <- paste(working.df$Salinity,working.df$Nitrogen)
plot(rep.mds.all)
