#############################################################################################################
#recreates figure 3 from Ricks et al
#############################################################################################################

#load libraries
library(scales)

#read in data and subset out
all_data <- read.csv('/path/to/microbial_growth_rates.csv')
plant_growth <- read.csv('/path/to/estimated_strain_partner_quality.csv')
plant_growth_ns <- plant_growth[plant_growth$ContSaline=='Low' & !(plant_growth$StrainID %in% c("Ancestor","Sterile")),]
plant_growth_saline <- plant_growth[plant_growth$ContSaline=='High' & !(plant_growth$StrainID %in% c("Ancestor","Sterile")),]


#subset data to growth rates in each individual environment
high_salt_low_nitrogen <-  microbe_growth[,c(1:3,5)]
low_salt_low_nitrogen <-  microbe_growth[,c(1:3,4)]
colnames(high_salt_low_nitrogen)[4] <- "Microbe_r";
colnames(low_salt_low_nitrogen)[4] <- "Microbe_r";

#ensure plant and microbial growth rates dataframes are in the same order then merge
plant_growth_ns <- plant_growth_ns[order(plant_growth_ns$StrainID),]
plant_growth_saline <- plant_growth_saline[order(plant_growth_saline$StrainID),]

low_salt_low_nitrogen <- low_salt_low_nitrogen[order(low_salt_low_nitrogen$StrainID),]
high_salt_low_nitrogen <- high_salt_low_nitrogen[order(high_salt_low_nitrogen$StrainID),]

low_salt_low_nitrogen$PlantGrowth <- plant_growth_ns$Plant_growth_rate
high_salt_low_nitrogen$PlantGrowth <- plant_growth_saline$Plant_growth_rate



#create output png file
png('Figure_3.png', width = 3.14961, height =5, units = "in", res = 1200, pointsize = 9)

#panel layout
layout(matrix(c(rep(1,1),
                rep(2,1)),
              2, 1, byrow = T))
par(cex=1)

par(mar = c(1.5, 4, 2, 1),oma=c(3,0,0,0),las=1)


#############################################################################################################
#panel 1 will show correlation between growth rates in nonsaline environments
#############################################################################################################
#subset to just look and low salinlity and build associated model 
working.d <- low_salt_low_nitrogen
plot(NA,ylim=range(working.d$PlantGrowth),xlim=c(range(working.d$Microbe_r)),
     xlab=NA,ylab=NA,yaxt='n',xaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)


par(xpd=NA)
text(x.left-range.x/6.25, y.bottom+range.y/2, 
     expression(bold("Plant growth (" * Delta ~ "mm"^2 ~ "/ day)")), 
     cex=1, font=2, srt=90)

text(x.left-range.x/25,y.top+range.y/15,'A',font=2,cex=1.25,pos=4)


axis(side=1,labels=T,tck=.03,cex.axis=(.8),padj=-2.55)

axis(side=2,labels=T,hadj=.35,cex.axis=(.8),tck=.03)

par(xpd=F)
mod1 <- (lm(PlantGrowth~Microbe_r,working.d))
anova(mod1)
abline(mod1)
x.pred <- seq(-.05,.3,.001)
pred.mod <- predict(mod1,data.frame(Microbe_r=x.pred),interval='confidence')

polygon(x=c(x.pred,rev(x.pred)),y=c(pred.mod[,2],rev(pred.mod[,3])),col=alpha('grey80',0.5),border=NA)



points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='High' & working.d$HistS=='Low',],col='#5D93B9',pch=1,lwd=2)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='High' & working.d$HistS=='High',],col='#BE4F48',pch=1,lwd=2)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='Low' &  working.d$HistS=='Low',],col='#5D93B9',pch=16)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='Low' &  working.d$HistS=='High',],col='#BE4F48',pch=16)

par(xpd=NA)
text(x.left+range.x/2,y.top+range.y/25,'Low salinity',cex=1,font=2)
text(x.left, y.top-range.y/20, expression(italic(p)~'= 0.028'), pos=4, cex=1)


#############################################################################################################
#panel 2 follows same design as above, but for correlation in saline environments
#############################################################################################################

working.d <- high_salt_low_nitrogen
plot(NA,ylim=range(working.d$PlantGrowth),xlim=c(range(working.d$Microbe_r)),
     xlab=NA,ylab=NA,yaxt='n',xaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)



text(x.left+range.x/2,y.bottom-range.y/10,c('Microbial growth (r)'),cex=1,font=2)
text(x.left-range.x/6.25, y.bottom+range.y/2, 
     expression(bold("Plant growth (" * Delta ~ "mm"^2 ~ "/ day)")), 
     cex=1, font=2, srt=90)



text(x.left-range.x/25,y.top+range.y/15,'B',font=2,cex=1.25,pos=4)

axis(side=1,labels=T,tck=.03,cex.axis=(.8),padj=-2.55)
axis(side=2,labels=T,hadj=.35,cex.axis=(.8),tck=.03)

par(xpd=F)
mod1 <- (lm(PlantGrowth~Microbe_r,working.d))
anova(mod1)
abline(mod1)
x.pred <- seq(-.05,.3,.001)
pred.mod <- predict(mod1,data.frame(Microbe_r=x.pred),interval='confidence')
polygon(x=c(x.pred,rev(x.pred)),y=c(pred.mod[,2],rev(pred.mod[,3])),col=alpha('grey80',0.5),border=NA)


points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='High' & working.d$HistS=='Low',],col='#5D93B9',pch=1,lwd=2)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='High' & working.d$HistS=='High',],col='#BE4F48',pch=1,lwd=2)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='Low' &  working.d$HistS=='Low',],col='#5D93B9',pch=16)
points(PlantGrowth~Microbe_r,working.d[working.d$HistN=='Low' &  working.d$HistS=='High',],col='#BE4F48',pch=16)
par(xpd=NA)
text(x.left+range.x/2,y.top+range.y/25,'High salinity',cex=1,font=2)
text(x.left, y.top-range.y/20, expression(italic(p)~'= 0.008'), pos=4, cex=1)





text(x.left,y.bottom-range.y/5,"Hist. Salinity",cex=1,font=2,pos=4)
text(x.left+range.x/20,y.bottom-range.y/5-(1*range.y/16),"Low",cex=1,font=1,pos=4)
text(x.left+range.x/20,y.bottom-range.y/5-(2*range.y/16),"High",cex=1,font=1,pos=4)
points(rep(x.left+range.x/20,2),c(y.bottom-range.y/5-(1*range.y/16),y.bottom-range.y/5-(2*range.y/16)),
       col=c('#5D93B9','#BE4F48'),bg=c('#5D93B9','#BE4F48'),pch=21,lwd=1,cex=.95)

text(x.left+range.x/2,y.bottom-range.y/5,"Hist. Nitrogen",cex=1,font=2,pos=4)
text(x.left+range.x/20+range.x/2,y.bottom-range.y/5-(1*range.y/16),"Low",cex=1,font=1,pos=4)
text(x.left+range.x/20+range.x/2,y.bottom-range.y/5-(2*range.y/16),"High",cex=1,font=1,pos=4)
points(rep(x.left+range.x/20+range.x/2,2),c(y.bottom-range.y/5-(1*range.y/16),y.bottom-range.y/5-(2*range.y/16)),
       col=c('grey45'),bg=c('grey45','white'),pch=21,lwd=1,cex=.95)


dev.off()

