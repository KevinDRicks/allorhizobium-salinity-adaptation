#code takes the outputs from the gatk and breseq pipeline and
#transforms them into a single harmonized snp table, associating each snp with the proper treatments
#we additionally generate summary statistic here 
#output file is the snp_summary.csv, which is used in downstream analysis


setwd('G:/My Drive/Work/Toronto/SalineEvolution/gitrepo_prep/sequencing_analysis')

#read in sequencing data as well as IDs for each sample
n9var <- read.table('/path/to/gatk_snps.table',header=T)#130
layout <- read.csv('path/to/sequencing_labels.csv')

#rename columns to represent each ID
colnames(n9var) <- gsub("X","",colnames(n9var))

#reorder dataset based on those in the layout df
#this codes the data so each line now represents a mutation for each strain
reorder_var <- NULL
for(i in 2:nrow(layout)){
  colid <- which(layout$SequencingID[i]==colnames(n9var))
  
  for(j in 1:nrow(n9var)){
    if(n9var[j,colid]!=n9var$REF[j]){
      reorder_var <- rbind(reorder_var,
                           data.frame(layout[i,],REF=n9var$REF[j],VAR=n9var[j,colid],Cluster=n9var$CHROM[j],Position=n9var$POS[j]))
      
    }
  }
}

#read in breseq
n9var_bs <- read.table('breseq_snps.txt',header=T)


#run a for loop to extract out info on each identified snp
na.len <- rep(NA,nrow(n9var_bs))
bs_reorder <- data.frame(Strain=na.len,Nitrogen=na.len,Salinity=na.len,
                         SequencingNumber=na.len,SequencingID=na.len,REF=na.len,VAR=na.len,
                         Cluster=na.len,Position=na.len)
for(i in 1:nrow(n9var_bs)){
  temp.d <- layout[layout$SequencingID==n9var_bs$Evolved[i],]
  bs_reorder$Strain[i] <- temp.d$Strain
  bs_reorder$Nitrogen[i] <- temp.d$Nitrogen
  bs_reorder$Salinity[i] <- temp.d$Salinity
  bs_reorder$SequencingNumber[i] <- temp.d$SequencingNumber
  bs_reorder$SequencingID [i] <- temp.d$SequencingID 
  bs_reorder$Cluster [i] <- n9var_bs$Cluster[i] 
  bs_reorder$Position [i] <- n9var_bs$Position[i] 
  bs_reorder$REF [i] <- n9var_bs$Ancestral_state[i] 
  bs_reorder$VAR [i] <- n9var_bs$Evolved_state[i] 
}


#extract our polymorphisms thtat are either snps, or large scale deletions
bs_small <- bs_reorder[nchar(bs_reorder$REF)<1000,]
large_deletions <- bs_reorder[nchar(bs_reorder$REF)>1000,]
large_deletions$EndPosition <- large_deletions$Position+nchar(large_deletions$REF)


#adds snps from breseq snp calls if its not the same as the gatk pipeline calls
#to make a single dataframe with all of the calls
merge_snps <- reorder_var
for(i in 1:nrow(bs_small)){
  if(nrow(merge_snps[merge_snps$Strain==bs_small$Strain[i] & 
                     merge_snps$Cluster==bs_small$Cluster[i] & 
                     merge_snps$Position==bs_small$Position[i] & 
                     merge_snps$VAR==bs_small$VAR[i],])==0){
    merge_snps <- rbind(merge_snps,bs_small[i,])
  }
}
dim(merge_snps)
head(merge_snps)
tail(merge_snps)

#gatk frequently called single deletions where breseq called whole plasmid deletions
#clean those up and remove from sample 
for(i in 1:nrow(large_deletions)){
  merge_snps <- merge_snps[-which(merge_snps$Strain==large_deletions[i,]$Strain & merge_snps$Cluster==large_deletions[i,]$Cluster & 
                                    merge_snps$Position>large_deletions[i,]$Position & merge_snps$Position<large_deletions[i,]$EndPosition),]
}
#merge_snps <- merge_snps[merge_snps$Expected_Species=='NF9',]
dim(merge_snps)
head(merge_snps)
tail(merge_snps)


#alter all point deletions that are code as a two bp sequence in reference to just two separate bp
adj.del <- merge_snps[nchar(merge_snps$REF)==2 & nchar(merge_snps$VAR)==1,]
adj.del_new <- NULL
for(i in 1:nrow(adj.del)){
  temp.c <- adj.del[i,]
  temp.c$Position <- temp.c$Position+1
  temp.c$VAR <- '.'
  temp.c$REF <- substring(temp.c$REF,2,2)
  adj.del_new <- rbind(adj.del_new,temp.c)
}
merge_snps <- merge_snps[!(nchar(merge_snps$REF)==2 & nchar(merge_snps$VAR)==1),]
merge_snps <- rbind(merge_snps,adj.del_new)
dim(merge_snps)

#alter adjacent point deletions that are code as two letters in reference to just single letter lines
adj.doub <- merge_snps[nchar(merge_snps$REF)==2 & nchar(merge_snps$VAR)==2,]
adj.doub_new <- NULL
for(i in 1:nrow(adj.doub)){
  temp.1 <- adj.doub[i,]
  temp.2 <- adj.doub[i,]
  temp.1$REF <- substring(temp.1$REF,1,1)
  temp.1$VAR <- substring(temp.1$VAR,1,1)
  temp.2$REF <- substring(temp.2$REF,2,2)
  temp.2$VAR <- substring(temp.2$VAR,2,2)
  temp.2$Position <- temp.2$Position+1
  adj.doub_new <- rbind(adj.doub_new,temp.1,temp.2)
}
merge_snps <- merge_snps[!(nchar(merge_snps$REF)==2 & nchar(merge_snps$VAR)==2),]
merge_snps <- rbind(merge_snps,adj.doub_new)
dim(merge_snps)

#alter adjacent point insertions that are coded as two letters in variant to just single letter 
adj.insert <- merge_snps[nchar(merge_snps$REF)==1 & nchar(merge_snps$VAR)==2,]
adj.insert_new <- NULL
for(i in 1:nrow(adj.insert)){
  temp.1 <- adj.insert[i,]
  temp.1$REF <- '.'
  temp.1$VAR <- substring(temp.1$VAR,2,2)
  adj.insert_new <- rbind(adj.insert_new,temp.1)
}
merge_snps <- merge_snps[!(nchar(merge_snps$REF)==1 & nchar(merge_snps$VAR)==2),]
merge_snps <- rbind(merge_snps,adj.insert_new)
dim(merge_snps)

#make sure that there are no duplicates
new.subset <- merge_snps[1,]
for(i in 2:nrow(merge_snps)){
  if(sum(new.subset$Strain==merge_snps$Strain[i] & new.subset$Position==merge_snps$Position[i])==0){
    new.subset <- rbind(new.subset,merge_snps[i,])
  }
}
dim(new.subset)

write.csv(new.subset,"snp_summary.csv")



#summary stats of snps by treatment
#cluster 1 = chromosome
#cluster 2 = plasmid A
#cluster 3 = plasmid B
#cluster 4 = plasmid C
#cluster 5 = plasmid D

new.subset$Code <- paste(new.subset$Salinity,new.subset$Nitrogen)

#snp count by treatment and replicon
#we adjust this based on the number of observed replicons in a treatment
#so for example, all treatments had 10 replicate strains, with the exception of the high salt high nitrogen with 7
#moreover, many of these strains were missing plasmids, so we can adjust the snp count by plasmid occurence
#snps across full genome
tapply(new.subset$Code,new.subset$Code,length)/c(7,10,10,10)

#snps across each replicon
tapply(new.subset[new.subset$Cluster=="cluster_001_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_001_consensus",]$Code,length)/c(7,10,10,10)

tapply(new.subset[new.subset$Cluster=="cluster_002_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_002_consensus",]$Code,length)/c(7,10,10)

#only 9 instances of this plasmid in low salt, low nitrogen, and 4 instances in low salt, high nitrogen
tapply(new.subset[new.subset$Cluster=="cluster_003_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_003_consensus",]$Code,length)/c(9,4,10)

#only 6 instance of plasmid in high salt, low nitorgen
tapply(new.subset[new.subset$Cluster=="cluster_004_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_004_consensus",]$Code,length)/c(7,6,10,10)

tapply(new.subset[new.subset$Cluster=="cluster_005_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_005_consensus",]$Code,length)/c(10)



#following the above, we can calculate the average snp frequency, based on the 
#size of the replicon. This represents the rate at which mutations are being introduced
#while not being purified out
#divide by replicon size
#full genome
tapply(new.subset$Code,new.subset$Code,length)/c(7,10,10,10)/5275905

#snps across each replicon
tapply(new.subset[new.subset$Cluster=="cluster_001_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_001_consensus",]$Code,length)/c(7,10,10,10)/4582387

tapply(new.subset[new.subset$Cluster=="cluster_002_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_002_consensus",]$Code,length)/c(7,10,10)/289046

#only 9 instances of this plasmid in low salt, low nitrogen, and 4 instances in low salt, high nitrogen
tapply(new.subset[new.subset$Cluster=="cluster_003_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_003_consensus",]$Code,length)/c(9,4,10)/213436

#only 6 instance of plasmid in high salt, low nitorgen
tapply(new.subset[new.subset$Cluster=="cluster_004_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_004_consensus",]$Code,length)/c(7,6,10,10)/162024

tapply(new.subset[new.subset$Cluster=="cluster_005_consensus",]$Code,
       new.subset[new.subset$Cluster=="cluster_005_consensus",]$Code,length)/c(10)/29012 


