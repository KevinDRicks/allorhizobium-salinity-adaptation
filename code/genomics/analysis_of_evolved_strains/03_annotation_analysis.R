#code takes all the SNPS and identifies functions of the genes the fall on using the output from interproscan
#specifically focuses on the annotation of go terms

#load packages
library(readr)
library(data.table)
library(vegan)


##############################################################################################
#read in prokka output and process into readable form
#while not using prokka annotations, we used prokka initially to identify proteins (with the built in prodigal function)  
#our input into the interpro scan consequently used the protein output file from prokka
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

#this dataframe can be used going forward
#provides codes for the start and end of each predicted protein
#as well as the ID, which can be referenced in the interproscan output
head(prk_outline_clean)


##############################################################################################
#read in interpro output and clean up to a readable form
#there are separate outputs for each replicon
##############################################################################################

pl1_inter <- read.table("/path/to/interpro/plasmid1_interpro.tsv",sep = "\t", comment.char = "", quote = "")
pl2_inter <- read.table("/path/to/interpro/plasmid2_interpro.tsv",sep = "\t", comment.char = "", quote = "")
pl3_inter <- read.table("/path/to/interpro/plasmid3_interpro.tsv",sep = "\t", comment.char = "", quote = "")
pl4_inter <- read.table("/path/to/interpro/plasmid4_interpro.tsv",sep = "\t", comment.char = "", quote = "")
chr_inter <- read.table("chromosome_interpro.tsv",sep = "\t", comment.char = "", quote = "")
interpro2go = fread("/path/to/interpro/interpro2go.txt", skip = 5, header = FALSE)

colnames(interpro2go) = c("base_desc", "GO_ID")

interpro2go$IPR_ID = substr(interpro2go$base_desc, 0, 18)

go_hold = strsplit(interpro2go$base_desc, " > ")
interpro2go$GO_desc = lapply(go_hold, function(x) x[2])


colnames_inter <- c("fasta_header", "md5","seq_length","analysis","signature_accession","signature_description",
                    "start_loc","stop_loc","score","match_status","run_date","interpro_accesion","interpro_description",
                    "go_annotations","pathways")

colnames(pl1_inter) <- colnames_inter;colnames(pl2_inter) <- colnames_inter
colnames(pl3_inter) <- colnames_inter;colnames(pl4_inter) <- colnames_inter
colnames(chr_inter) <- colnames_inter

pl1_inter$GO_title <- NA
pl2_inter$GO_title <- NA
pl3_inter$GO_title <- NA
pl4_inter$GO_title <- NA
chr_inter$GO_title <- NA


#we create a function to extract out the go terms
go_annote_fun <- function(interpro_df){
  working_data <- interpro_df
  
  for(j in 1:nrow(working_data)){
    #only work on go annotations
    if(working_data$go_annotations[j]!="-"){
      #some of these terms have more than 1 prediction
      go_list <- unlist(strsplit(working_data$go_annotations[j],'\\|'))
      go_list <- substr(go_list,1,10)
      
      
      title <- NULL
      for(k in 1:length(go_list)){
        title <- c(title,unique(unlist(interpro2go[interpro2go$GO_ID==go_list[k], "GO_desc"])))
      }
      if(length(title)!=0){
        working_data$GO_title[j] <- do.call(paste, c(as.list(title), sep = "; "))
      }
    }
    print(j)
  }
  return(working_data)
}


pl1_inter <- go_annote_fun(pl1_inter)
pl2_inter <- go_annote_fun(pl2_inter)
pl3_inter <- go_annote_fun(pl3_inter)
pl4_inter <- go_annote_fun(pl4_inter)
chr_inter <- go_annote_fun(chr_inter)

#pool at annotations into one dataframe
merge_inter <- rbind(pl1_inter,pl2_inter,pl3_inter,
                     pl4_inter,chr_inter)
##############################################################################################
#read in the snp summary and ids
##############################################################################################
merge_snps <- read.csv("/path/to/snp_summary.csv")
layout <- read.csv('/path/to/sequencing_labels.csv')


#################################################################################################
#add in prokka label for each snp
#this can be cross referenced to the interpro downstream
#################################################################################################
merge_snps$prokka_lab <- NA
for(i in 1:nrow(merge_snps)){
  temp.prokka <- prk_outline_clean[prk_outline_clean$Cluster==merge_snps$Cluster[i],]
  prokka_assign <- temp.prokka[temp.prokka$Start<merge_snps$Position[i] & temp.prokka$End>merge_snps$Position[i],]
  if(nrow(prokka_assign)!=0){
    merge_snps$prokka_lab[i] <- prokka_assign$ID
  }
}

#################################################################################################
#finds go terms associated with each gene for snps
#################################################################################################
merge_snps$go_term_gene <- NA
merge_snps$go_term_gene_direct <- NA
for(i in 1:nrow(merge_snps)){
  inter_temp <- merge_inter[merge_inter$fasta_header==merge_snps$prokka_lab[i],]
  go_temp <- inter_temp$go_annotations
  go_temp <- go_temp[!is.na(go_temp)]
  go_temp <- go_temp[go_temp!='-']
  if(length(go_temp)!=0){
    go_temp <- unique(unlist(strsplit(go_temp,"\\|")))
    go_temp <- gsub("\\s*\\([^\\)]+\\)","",unique(unlist(strsplit(go_temp,"; "))))
    merge_snps$go_term_gene[i] <- do.call(paste, c(as.list(go_temp), sep = "; "))
  }
  #alternative way-directly pulling go terms
  go_temp <- inter_temp$GO_title
  go_temp <- go_temp[!is.na(go_temp)]
  go_temp <- gsub("GO:","",unique(unlist(strsplit(go_temp,"; "))))
  if(length(go_temp)!=0){
    merge_snps$go_term_gene_direct[i] <- do.call(paste, c(as.list(go_temp), sep = "; "))
  }
}
dim(merge_snps)


#################################################################################################
#extract higher order go terms for each gene
#this was not explicitly used in this analysis
#but could be if one wanted to look at higher go terms
#################################################################################################
library(GO.db)

#load go db
biop <- as.list(GOBPANCESTOR) # Remove GO IDs that do not have any ancestor 
ccp <- as.list(GOCCANCESTOR) # Remove GO IDs that do not have any ancestor 
mfp <- as.list(GOMFANCESTOR) # Remove GO IDs that do not have any ancestor 

#extract unique go terms observed
unique.go_terms <- unique(unlist(strsplit(merge_snps$go_term_gene,'; ')))
unique.go_terms <- unique.go_terms[!is.na(unique.go_terms)]

#create key with observed go terms
go_higher_key <- data.frame(GoTerm=unique.go_terms,GO_higher=NA)
for(i in 1:length(unique.go_terms)){
  bip.p <- unlist(biop[names(biop)==unique.go_terms[i]])  
  ccp.p <- unlist(ccp[names(ccp)==unique.go_terms[i]]) 
  mfp.p <- unlist(mfp[names(mfp)==unique.go_terms[i]])  
  all.p <- unique(c(bip.p,ccp.p,mfp.p))
  all.p <- all.p[all.p!='all']
  all.p <- c(all.p,unique.go_terms[i])
  
  go_higher_key$GO_higher[i] <-  do.call(paste, c(as.list(all.p), sep = "; "))
}


#fill observed snps with higher order terms based on key
merge_snps$go_gene_higher <- NA
for(i in 1:nrow(merge_snps)){
  temp.go <- unlist(strsplit(merge_snps$go_term_gene[i],'; '))
  collected.GO <- NULL
  if(length(temp.go)>0){
    for(j in 1:length(temp.go)){
      collected.GO <- c(collected.GO,unlist(strsplit(go_higher_key[go_higher_key$GoTerm==temp.go[j],]$GO_higher,'; ')))
    }
    collected.GO <- unique(collected.GO)
  }
  merge_snps$go_gene_higher[i] <- do.call(paste, c(as.list(collected.GO), sep = "; "))
}

merge_snps[merge_snps$go_gene_higher=='NA',]$go_gene_higher <- NA

length(unique(unlist(strsplit(merge_snps$go_gene_higher,'; '))))

#how many of the unique genes have a functional annotation
length(unique(merge_snps[!is.na(merge_snps$go_term_gene_direct),]$prokka_lab))


#################################################################################################
#prep for analysis
#################################################################################################
#ensure any snps from deleted plasmid (as a result of gatk pipeline unable to identify those deletion) have been removed
pl_b <- c("St3",  "St41", "St42", "St43", "St48", "St49", "St50")
pl_c <- c("St22", "St24", "St33", "St39")
new.subset <- merge_snps[!(merge_snps$Strain %in% pl_b & merge_snps$Cluster=="cluster_003_consensus"),]
new.subset <- new.subset[!(new.subset$Strain %in% pl_c & new.subset$Cluster=="cluster_004_consensus"),]
dim(new.subset)


#compile trait to interest
#in this case, the lowest level go term. Could replace this though with go_gene_higher to look at high levels
new.subset$variable_of_int <- new.subset$go_term_gene_direct
unq_go <- unique(unlist(strsplit(new.subset$variable_of_int,"; "))); unq_go <- unq_go[!is.na(unq_go)]
length(unq_go)
layout.adj <- layout[layout$Strain!='ancestor',]

#create dataframe to work with
ord.d <- data.frame(matrix(0,nrow=nrow(layout.adj),ncol=length(unq_go)))
row.names(ord.d) <- layout.adj$SequencingID
colnames(ord.d) <- unq_go
for(i in 1:nrow(new.subset)){
  working.t <- unlist(strsplit(new.subset$variable_of_int[i],"; "))
  if(!sum(is.na(working.t))>0){
    for(j in 1:length(working.t)){
      ord.d[row.names(ord.d)==new.subset$SequencingID[i],colnames(ord.d)==working.t[j]] <- 1
    }
  }
}
#remove terms with zeros
ord.d <- ord.d[,colSums(ord.d)>1]
dim(ord.d)

#run an permanova, to look at multivariate space
adonis2(ord.d~Salinity*Nitrogen ,layout.adj,
        method='manhattan',by='terms')
#significant salinity effect



#run univariate analysis on all observed annotation groups
#analysis runs a basic binomial glm as data is 0/1 
s.vec <- NULL
n.vec <- NULL
int.vec <- NULL
for(i in 1:ncol(ord.d)){
  p.vec <- unlist(anova(glm(ord.d[,i]~Salinity*Nitrogen ,family='binomial',layout.adj)))[c(18:20)]
  s.vec <- c(s.vec,p.vec[1])
  n.vec <- c(n.vec,p.vec[2])
  int.vec <- c(int.vec,p.vec[3])
  
}
which(s.vec<0.05)
which(n.vec<0.05)
which(int.vec<0.05)
var.int <- colnames(ord.d)

var.int
#each of these p values, corresponds with one of the annotation, order in the above vector, var.int
#these p values are extracted from the glm binomial model, using Salinity X Nitrogen, and are fdr corrected
p.adjust(s.vec,method='fdr')
p.adjust(n.vec,method='fdr')
p.adjust(int.vec,method='fdr')
