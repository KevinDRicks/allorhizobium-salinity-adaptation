#code analyzes the enrichment of specific cog functions on the two lost plasmids
#enrichment is relative to the main chromosome, as this may indiciate at some functional profile
#we focus on the main chromosome, and plasmid b and c, as these were the lost plasmids



##############################################################################
#reads in cog analysis (from eggnog) and massage the data into readable format
##############################################################################
cog.df <- read.csv("/path/to/cogs.csv",skip=1)
cog.df$query <- as.numeric(substr(substr(cog.df$query,1,11),9,11))

#breakdown of functions by replicon
table(cog.df$query)

#cog.df <- cog.df[cog.df$COG_category=='-',]

#cogs are a one letter code, however, the output includes some two or more letter codes
#this is because they've been assigned multiple potential cog groupings
#we extract the observed cogs total the number of observations / replicon

#extracts the observed cogs
obs.cogs <- unique(unlist(strsplit(unique(cog.df$COG_category),"")))

chr.dist <- data.frame(COGS=obs.cogs,Ab=NA)
pl2.dist <- data.frame(COGS=obs.cogs,Ab=NA)
pl3.dist <- data.frame(COGS=obs.cogs,Ab=NA)

#the following loops sum up the total number of observations in each cog profile
#in each of the replicons of interest
for(i in 1:length(obs.cogs)){
  temp.cog <- cog.df[cog.df$query==1,]
  temp.gene <- strsplit(temp.cog$COG_category,'')
  count.st <- 0
  for(j in 1:length(temp.gene)){
    temp.d <- unlist(temp.gene[j])
    temp.d <- temp.d[!is.na(temp.d)]
    if(sum(temp.d==obs.cogs[i])>0){
      count.st <- count.st+1
    }
  }
  chr.dist$Ab[i] <- count.st
}

for(i in 1:length(obs.cogs)){
  temp.cog <- cog.df[cog.df$query==3,]
  temp.gene <- strsplit(temp.cog$COG_category,'')
  count.st <- 0
  for(j in 1:length(temp.gene)){
    temp.d <- unlist(temp.gene[j])
    temp.d <- temp.d[!is.na(temp.d)]
    if(sum(temp.d==obs.cogs[i])>0){
      count.st <- count.st+1
    }
  }
  pl2.dist$Ab[i] <- count.st
}

for(i in 1:length(obs.cogs)){
  temp.cog <- cog.df[cog.df$query==4,]
  temp.gene <- strsplit(temp.cog$COG_category,'')
  count.st <- 0
  for(j in 1:length(temp.gene)){
    temp.d <- unlist(temp.gene[j])
    temp.d <- temp.d[!is.na(temp.d)]
    if(sum(temp.d==obs.cogs[i])>0){
      count.st <- count.st+1
    }
  }
  pl3.dist$Ab[i] <- count.st
}

#estimates the relative abundance of each cog category for each replicon
chr.dist$Rel_ab <- chr.dist$Ab/sum(chr.dist$Ab)
pl2.dist$Rel_ab <- pl2.dist$Ab/sum(pl2.dist$Ab)
pl3.dist$Rel_ab <- pl3.dist$Ab/sum(pl3.dist$Ab)

#compares the enrichment rate of each category, relative to the chromosome
#will range from -1 to 1
pl2.dist$Alt_comp <- (pl2.dist$Rel_ab-chr.dist$Rel_ab)/(pl2.dist$Rel_ab+chr.dist$Rel_ab)
pl3.dist$Alt_comp <- (pl3.dist$Rel_ab-chr.dist$Rel_ab)/(pl3.dist$Rel_ab+chr.dist$Rel_ab)

#alterative calculation for enrichment rate
#pl2.dist$Comp_ab <- (pl2.dist$Rel_ab+1)/(chr.dist$Rel_ab+1)
#pl3.dist$Comp_ab <- (pl3.dist$Rel_ab+1)/(chr.dist$Rel_ab+1)

#we exclude the following cogs as they are largely associated with replication
#plasmids are almost universally enriched in such genes, which is a part of the inherent plasmid life history
#to derrive some idea of plasmid function, we ignore these
#R, S are poorly characterized functions
#- are uncaterogized
#L-dna replication/recombination/repair, K-transcription, D-cell division/chromosome partitioning
#B-chromatin structure, J-translation
pl2.dist <- pl2.dist[!(pl2.dist$COGS=='R' |pl2.dist$COGS=='S' | pl2.dist$COGS=='-' | 
                         pl2.dist$COGS=='L' | pl2.dist$COGS=='K' | pl2.dist$COGS=='D' | 
                         pl2.dist$COGS=='B' | pl2.dist$COGS=='J'),]

pl3.dist <- pl3.dist[!(pl3.dist$COGS=='R' |pl3.dist$COGS=='S' | pl3.dist$COGS=='-' | 
                         pl3.dist$COGS=='L' | pl3.dist$COGS=='K' | pl3.dist$COGS=='D' |
                         pl3.dist$COGS=='B' | pl3.dist$COGS=='J'),]




