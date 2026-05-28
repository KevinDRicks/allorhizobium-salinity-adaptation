#code takes all the SNPS from the evolved strains and estimates enrichment and differentiation of SNPs 
#on gene by gene basis between treatments

#load packages
library(readr)
library(data.table)
library(hierfstat)

##############################################################################################
#read in the snp summary and ids
##############################################################################################
merge_snps <- read.csv("/path/to/snp_summary.csv")
layout <- read.csv('/path/to/sequencing_labels.csv')



##############################################################################################
#read in prokka output and process into readable form
#while not using prokka annotations, we used prokka initially to identify proteins (with the built in prodigal function)  
#we use this here for keeping consistent labels 
##############################################################################################
prk_outline <- readLines('/path/to/prokka/PROKKA_05052025.gff')
prk_outline_clean <- NULL
for(i in 1:length(prk_outline)){
  temp <- prk_outline[i]
  temp.work <- unlist(strsplit(prk_outline[i],"\t"))
  prk_outline_clean <- rbind(prk_outline_clean,data.frame(Cluster=temp.work[1],Start=temp.work[4],End=temp.work[5],
                                                          ID=gsub("ID=","",unlist(strsplit(temp.work[9],';'))[1]),
                                                          Prokka=unlist(strsplit(temp.work[9],';'))[4]))
}
prokka_sum <- data.frame(read_tsv("/path/to/prokka/PROKKA_05052025.tsv"))
for(i in 1:nrow(prk_outline_clean)){
  prk_outline_clean$Prokka[i] <- prokka_sum[prokka_sum$locus_tag==prk_outline_clean$ID[i],]$product
}
prk_outline_clean$Start <- as.numeric(prk_outline_clean$Start)
prk_outline_clean$End <- as.numeric(prk_outline_clean$End)


#identify the gene closest to each snp
merge_snps$prokka_lab <- NA
for(i in 1:nrow(merge_snps)){
  
  temp.prokka <- prk_outline_clean[prk_outline_clean$Cluster==merge_snps$Cluster[i],]
  
  prokka_assign <- temp.prokka[temp.prokka$Start<merge_snps$Position[i] & temp.prokka$End>merge_snps$Position[i],]
  if(nrow(prokka_assign)!=0){
    merge_snps$prokka_lab[i] <- prokka_assign$ID
  }
  
  if(nrow(prokka_assign)==0){
    
    front.min <- min(abs(temp.prokka$Start-merge_snps$Position[i]))
    rear.min <- min(abs(temp.prokka$End-merge_snps$Position[i]))
    if(front.min<rear.min){
      prokka_assign <- temp.prokka[which( abs(temp.prokka$Start-merge_snps$Position[i])==front.min),]
      merge_snps$prokka_lab[i] <- prokka_assign$ID

      prokka_assign <- temp.prokka[which( abs(temp.prokka$End-merge_snps$Position[i])==rear.min),]
    }
    
    if(front.min>rear.min){
      prokka_assign <- temp.prokka[which( abs(temp.prokka$End-merge_snps$Position[i])==rear.min),]
      merge_snps$prokka_lab[i] <- prokka_assign$ID

      prokka_assign <- temp.prokka[which( abs(temp.prokka$Start-merge_snps$Position[i])==front.min),]
    }
  }
  
  print(i)
}


##############################################################################################
#thins data to only 1 snp per gene for a given strain
#given the scarcity of snps, we'll base fst on any snp detected within a gene
##############################################################################################
merge_snps.new <- merge_snps[1,]
for(i in 2:nrow(merge_snps)){
  if(sum(merge_snps$Strain[i]==merge_snps.new$Strain & merge_snps$prokka_lab[i]==merge_snps.new$prokka_lab)==0){
    merge_snps.new <- rbind(merge_snps.new,merge_snps[i,])
  }
}

#find the unique genes within each treatment group with snps, and the number of strains
LowNLowS.table <- table(merge_snps.new[merge_snps.new$Nitrogen=='Low' & merge_snps.new$Salinity=='Low',]$prokka_lab )
LowNHighS.table <- table(merge_snps.new[merge_snps.new$Nitrogen=='Low' & merge_snps.new$Salinity=='High',]$prokka_lab)
HighNLowS.table <- table(merge_snps.new[merge_snps.new$Nitrogen=='High' & merge_snps.new$Salinity=='Low',]$prokka_lab)
HighNHighS.table <- table(merge_snps.new[merge_snps.new$Nitrogen=='High' & merge_snps.new$Salinity=='High',]$prokka_lab)

#reorder
LowNLowS.table <- LowNLowS.table[order(names(LowNLowS.table))]
LowNHighS.table <- LowNHighS.table[order(names(LowNHighS.table))]
HighNLowS.table <- HighNLowS.table[order(names(HighNLowS.table))]
HighNHighS.table <- HighNHighS.table[order(names(HighNHighS.table))]

#collates all the unique genes with snps into one table
collate.snp.gene <- data.frame(Gene=unique(merge_snps$prokka_lab),
                               LowNLowS=0,LowNHighS=0,
                               HighNLowS=0,HighNHighS=0)
collate.snp.gene <- collate.snp.gene[order(collate.snp.gene$Gene),]

collate.snp.gene[collate.snp.gene$Gene %in% names(LowNLowS.table),]$LowNLowS <- LowNLowS.table
collate.snp.gene[collate.snp.gene$Gene %in% names(LowNHighS.table),]$LowNHighS <- LowNHighS.table
collate.snp.gene[collate.snp.gene$Gene %in% names(HighNLowS.table),]$HighNLowS <- HighNLowS.table
collate.snp.gene[collate.snp.gene$Gene %in% names(HighNHighS.table),]$HighNHighS <- HighNHighS.table


#identify what replicon each gene is on an their position 
collate.snp.gene$Replicon <- NA
collate.snp.gene$Position <- NA
for(i in 1:nrow(collate.snp.gene)){
  collate.snp.gene$Replicon[i] <- prk_outline_clean[prk_outline_clean$ID==collate.snp.gene$Gene[i],]$Cluster
  collate.snp.gene$Position[i] <- prk_outline_clean[prk_outline_clean$ID==collate.snp.gene$Gene[i],]$Start
  
}

#subset layoutfile and encode which ones have lost plasmids
working.df <- layout[layout$Strain!='ancestor',]
working.df$Pl2_loss <- 0; working.df$Pl3_loss <- 0;
pl_b <- c("St3",  "St41", "St42", "St43", "St48", "St49", "St50")
pl_c <- c("St22", "St24", "St33", "St39")


working.df[working.df$Strain  %in% pl_b,]$Pl2_loss <- 1
working.df[working.df$Strain  %in% pl_c,]$Pl3_loss <- 1


#code the number of samples per treatment groups we'll be using for the enrichment calculation 
#based on the number of samples
#note, because some strains have lost plasmids, we will decrease samples in those treatments
#hard code replication in below based on sampling numbers for each replicon
collate.snp.gene$LowNLowS_samples <- NA;collate.snp.gene$LowNHighS_samples <- NA;
collate.snp.gene$HighNLowS_samples <- NA;collate.snp.gene$HighNHighS_samples <- NA;

collate.snp.gene[collate.snp.gene$Replicon=='cluster_001_consensus' | collate.snp.gene$Replicon=='cluster_002_consensus' | collate.snp.gene$Replicon=='cluster_004_consensus' | 
                   collate.snp.gene$Replicon=='cluster_005_consensus',]$LowNLowS_samples <- 10
collate.snp.gene[collate.snp.gene$Replicon=='cluster_003_consensus',]$LowNLowS_samples <- 9

collate.snp.gene[collate.snp.gene$Replicon=='cluster_001_consensus' | collate.snp.gene$Replicon=='cluster_002_consensus' | collate.snp.gene$Replicon=='cluster_003_consensus' | 
                   collate.snp.gene$Replicon=='cluster_005_consensus',]$LowNHighS_samples <- 10
collate.snp.gene[collate.snp.gene$Replicon=='cluster_004_consensus',]$LowNHighS_samples <- 6

collate.snp.gene[collate.snp.gene$Replicon=='cluster_001_consensus' | collate.snp.gene$Replicon=='cluster_002_consensus' | collate.snp.gene$Replicon=='cluster_004_consensus' | 
                   collate.snp.gene$Replicon=='cluster_005_consensus',]$HighNLowS_samples <- 10
collate.snp.gene[collate.snp.gene$Replicon=='cluster_003_consensus',]$HighNLowS_samples <- 4

collate.snp.gene$HighNHighS_samples <- 7


#calculates differentiation in SNP enrichment between treatment groups
collate.snp.gene$Nitrogen <- NA
collate.snp.gene$Salinity <- NA
collate.snp.gene$Nit_lowS <- NA
collate.snp.gene$Nit_highS <- NA
collate.snp.gene$Sal_lowN <- NA
collate.snp.gene$Sal_highN <- NA
for(i in 1:nrow(collate.snp.gene)){
  lnls <- rep(0,collate.snp.gene$LowNLowS_samples[i])
  lnls[c(0:collate.snp.gene$LowNLowS [i])] <- 1
  
  lnhs <- rep(0,collate.snp.gene$LowNHighS_samples[i])
  lnhs[c(0:collate.snp.gene$LowNHighS [i])] <- 1
  
  hnls <- rep(0,collate.snp.gene$HighNLowS_samples[i])
  hnls[c(0:collate.snp.gene$HighNLowS [i])] <- 1
  
  hnhs <- rep(0,collate.snp.gene$HighNHighS_samples[i])
  hnhs[c(0:collate.snp.gene$HighNHighS [i])] <- 1
  
  test.gene <- data.frame(Nitrogen=c(rep('Low',collate.snp.gene$LowNLowS_samples[i]),
                                     rep('Low',collate.snp.gene$LowNHighS_samples[i]),
                                     rep('High',collate.snp.gene$HighNLowS_samples[i]),
                                     rep('High',collate.snp.gene$HighNHighS_samples[i])),
                          Salinity=c(rep('Low',collate.snp.gene$LowNLowS_samples[i]),
                                     rep('High',collate.snp.gene$LowNHighS_samples[i]),
                                     rep('Low',collate.snp.gene$HighNLowS_samples[i]),
                                     rep('High',collate.snp.gene$HighNHighS_samples[i])),
                          Gene=c(lnls,lnhs,hnls,hnhs))
  
  collate.snp.gene$Salinity[i] <- as.numeric(genet.dist(cbind(test.gene$Salinity,test.gene$Gene),
                                                        method='Fst', diploid=F))
  collate.snp.gene$Nitrogen[i] <- as.numeric(genet.dist(cbind(test.gene$Nitrogen,test.gene$Gene),
                                                        method='Fst', diploid=F))
  
  collate.snp.gene$Nit_lowS[i] <- as.numeric(genet.dist(cbind(test.gene[test.gene$Salinity=='Low',]$Nitrogen,
                                                              test.gene[test.gene$Salinity=='Low',]$Gene),
                                                        method='Fst', diploid=F))
  collate.snp.gene$Nit_highS[i] <- as.numeric(genet.dist(cbind(test.gene[test.gene$Salinity=='High',]$Nitrogen,
                                                               test.gene[test.gene$Salinity=='High',]$Gene),
                                                         method='Fst', diploid=F))
  collate.snp.gene$Sal_lowN[i] <- as.numeric(genet.dist(cbind(test.gene[test.gene$Nitrogen=='Low',]$Salinity,
                                                              test.gene[test.gene$Nitrogen=='Low',]$Gene),
                                                        method='Fst', diploid=F))
  collate.snp.gene$Sal_highN[i] <- as.numeric(genet.dist(cbind(test.gene[test.gene$Nitrogen=='High',]$Salinity,
                                                               test.gene[test.gene$Nitrogen=='High',]$Gene),
                                                         method='Fst', diploid=F))
}

collate.snp.gene$LowNLowS_percent <- collate.snp.gene$LowNLowS/collate.snp.gene$LowNLowS_samples
collate.snp.gene$LowNHighS_percent <- collate.snp.gene$LowNHighS/collate.snp.gene$LowNHighS_samples
collate.snp.gene$HighNLowS_percent <- collate.snp.gene$HighNLowS/collate.snp.gene$HighNLowS_samples
collate.snp.gene$HighNHighS_percent <- collate.snp.gene$HighNHighS/collate.snp.gene$HighNHighS_samples

collate.snp.gene.thin <- collate.snp.gene[,c(1,6,7,12:21)]





