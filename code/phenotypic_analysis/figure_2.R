#############################################################################################################
#recreates figure 2 from Ricks et al
#############################################################################################################
#read in data
microbe_growth <- read.csv("/path/to/microbial_growth_rates.csv")

#subset data to growth rates in each individual environment
high_salt_high_nitrogen <- microbe_growth[,c(1:3,7)]
high_salt_low_nitrogen <-  microbe_growth[,c(1:3,5)]
low_salt_high_nitrogen <-  microbe_growth[,c(1:3,6)]
low_salt_low_nitrogen <-  microbe_growth[,c(1:3,4)]






#create output png file
png('Figure_2.png', width = 4.14961, height =5, units = "in", res = 1200, pointsize = 9)

#panel layout
layout(matrix(c(1,2,3,4),
              2, 2, byrow = T))
par(cex=1)
par(mar = c(1.5, 2, 1.75, .5),oma=c(3,1.5,0,0),las=1)

#set growth environment of interest to temp so we can consistently use this terminology
temp <- low_salt_low_nitrogen

#recode lables of histN and histS for quick and dirty method to alphabetize (important if we want figures in order)
temp[temp$HistN=='Low',]$HistN <- 'a_Low'
temp[temp$HistS=='Low',]$HistS <- 'a_Low'

#create plot and fill in with details 
plot(NA,xlim=c(0.5,4.5),ylim=c(0.007,0.012),xlab=NA,ylab=NA,yaxt='n',xaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
par(xpd=NA)

axis(side=1,labels=F,tck=.02,cex=1)
axis(side=2,labels=T,hadj=.5,cex.axis=(.8),tck=.02)

text(x.left-range.x/4.5,y.bottom+range.y/2,c('Microbial growth rate (r)'),cex=1,font=2,srt=90)
text(x.left-range.x/30,y.top+range.y/30,'Low nitrogen, Low salinity',pos=4,font=2)


#add in estimate of growth for each group, just using simple mean and std error estimates
n.vec <- c('a_Low','High','a_Low','High')
s.vec <- c('a_Low','a_Low','High','High')
col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
alpha.lvl <- .5
col.vec.alph <-c(alpha('#5D93B9',alpha.lvl),alpha('#5D93B9',alpha.lvl),alpha('#BE4F48',alpha.lvl),alpha('#BE4F48',alpha.lvl))
bg.vec.alph <- c(alpha('#5D93B9',alpha.lvl),'white',alpha('#BE4F48',alpha.lvl),'white')
set.seed(50)
for(i in 1:4){
  temp.d <- temp[temp$HistN==n.vec[i] & temp$HistS==s.vec[i],]
  temp.mean <- mean(temp.d$Microbegrowth_r_LowN_LowS)
  temp.se <- sd(temp.d$Microbegrowth_r_LowN_LowS)/sqrt(nrow(temp.d))
  segments(x0=i,x1=i,y0=temp.mean-1.96*temp.se,y1=temp.mean+1.96*temp.se,col=col.vec[i],lwd=2)
  points(i,temp.mean,pch=21,col=col.vec[i],bg=bg.vec[i],lwd=2,cex=1.5)
  
}
text(x.left-range.x/30,y.top+range.y/10,'A',font=2,cex=1.25,pos=2)

#############################################################################################################
#we repeat the above scripts for each environment
#high salt low nitrogen
#############################################################################################################
temp <- high_salt_low_nitrogen


temp[temp$HistN=='Low',]$HistN <- 'a_Low'
temp[temp$HistS=='Low',]$HistS <- 'a_Low'

plot(NA,xlim=c(0.5,4.5),ylim=c(0.004,0.007),xlab=NA,ylab=NA,yaxt='n',xaxt='n')

rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
par(xpd=NA)

axis(side=1,labels=F,tck=.02,cex=1)
axis(side=2,labels=T,hadj=.5,cex.axis=(.8),tck=.02,at=c(.004,.005,.006,.007))
text(x.left-range.x/30,y.top+range.y/30,'Low nitrogen, High salinity',pos=4,font=2)




n.vec <- c('a_Low','High','a_Low','High')
s.vec <- c('a_Low','a_Low','High','High')
col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
alpha.lvl <- .5
col.vec.alph <-c(alpha('#5D93B9',alpha.lvl),alpha('#5D93B9',alpha.lvl),alpha('#BE4F48',alpha.lvl),alpha('#BE4F48',alpha.lvl))
bg.vec.alph <- c(alpha('#5D93B9',alpha.lvl),'white',alpha('#BE4F48',alpha.lvl),'white')
set.seed(50)
for(i in 1:4){
  temp.d <- temp[temp$HistN==n.vec[i] & temp$HistS==s.vec[i],]
  temp.mean <- mean(temp.d$Microbegrowth_r_LowN_HighS)
  temp.se <- sd(temp.d$Microbegrowth_r_LowN_HighS)/sqrt(nrow(temp.d))
  segments(x0=i,x1=i,y0=temp.mean-1.96*temp.se,y1=temp.mean+1.96*temp.se,col=col.vec[i],lwd=2)
  points(i,temp.mean,pch=21,col=col.vec[i],bg=bg.vec[i],lwd=2,cex=1.5)
  
}

text(x.left-range.x/30,y.top+range.y/10,'B',font=2,cex=1.25,pos=2)

#############################################################################################################
#we repeat the above scripts for each environment
#low salt high nitrogen
#############################################################################################################
temp <- low_salt_high_nitrogen


temp[temp$HistN=='Low',]$HistN <- 'a_Low'
temp[temp$HistS=='Low',]$HistS <- 'a_Low'

plot(NA,xlim=c(0.5,4.5),ylim=c(0.008,0.020),xlab=NA,ylab=NA,yaxt='n',xaxt='n')

rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
par(xpd=NA)

axis(side=1,labels=F,tck=.02,cex=1)
axis(side=2,labels=T,hadj=.5,cex.axis=(.8),tck=.02)

text(x.left-range.x/4.5,y.bottom+range.y/2,c('Microbial growth (r)'),cex=1,font=2,srt=90)
text(x.left+range.x/2,y.bottom-range.y/20,c('Evolved strains'),cex=1,font=2)
text(x.left-range.x/30,y.top+range.y/30,'High nitrogen, Low salinity',pos=4,font=2)



n.vec <- c('a_Low','High','a_Low','High')
s.vec <- c('a_Low','a_Low','High','High')
col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
alpha.lvl <- .5
col.vec.alph <-c(alpha('#5D93B9',alpha.lvl),alpha('#5D93B9',alpha.lvl),alpha('#BE4F48',alpha.lvl),alpha('#BE4F48',alpha.lvl))
bg.vec.alph <- c(alpha('#5D93B9',alpha.lvl),'white',alpha('#BE4F48',alpha.lvl),'white')
set.seed(50)
for(i in 1:4){
  temp.d <- temp[temp$HistN==n.vec[i] & temp$HistS==s.vec[i],]
  temp.mean <- mean(temp.d$Microbegrowth_r_HighN_LowS)
  temp.se <- sd(temp.d$Microbegrowth_r_HighN_LowS)/sqrt(nrow(temp.d))
  segments(x0=i,x1=i,y0=temp.mean-1.96*temp.se,y1=temp.mean+1.96*temp.se,col=col.vec[i],lwd=2)
  points(i,temp.mean,pch=21,col=col.vec[i],bg=bg.vec[i],lwd=2,cex=1.5)
  
}
text(x.left-range.x/30,y.top+range.y/10,'C',font=2,cex=1.25,pos=2)


#############################################################################################################
#we repeat the above scripts for each environment
#high salt high nitrogen
#############################################################################################################
temp <- high_salt_high_nitrogen


temp[temp$HistN=='Low',]$HistN <- 'a_Low'
temp[temp$HistS=='Low',]$HistS <- 'a_Low'

plot(NA,xlim=c(0.5,4.5),ylim=c(0.008,0.013),xlab=NA,ylab=NA,yaxt='n',xaxt='n')

rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
par(xpd=NA)

axis(side=1,labels=F,tck=.02,cex=1)
axis(side=2,labels=T,hadj=.5,cex.axis=(.8),tck=.02)

text(x.left+range.x/2,y.bottom-range.y/20,c('Evolved strains'),cex=1,font=2)
text(x.left-range.x/30,y.top+range.y/30,'High nitrogen, High salinity',pos=4,font=2)



n.vec <- c('a_Low','High','a_Low','High')
s.vec <- c('a_Low','a_Low','High','High')
col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
alpha.lvl <- .5
col.vec.alph <-c(alpha('#5D93B9',alpha.lvl),alpha('#5D93B9',alpha.lvl),alpha('#BE4F48',alpha.lvl),alpha('#BE4F48',alpha.lvl))
bg.vec.alph <- c(alpha('#5D93B9',alpha.lvl),'white',alpha('#BE4F48',alpha.lvl),'white')
set.seed(50)
for(i in 1:4){
  temp.d <- temp[temp$HistN==n.vec[i] & temp$HistS==s.vec[i],]
  temp.mean <- mean(temp.d$Microbegrowth_r_HighN_HighS)
  temp.se <- sd(temp.d$Microbegrowth_r_HighN_HighS)/sqrt(nrow(temp.d))
  segments(x0=i,x1=i,y0=temp.mean-1.96*temp.se,y1=temp.mean+1.96*temp.se,col=col.vec[i],lwd=2)
  points(i,temp.mean,pch=21,col=col.vec[i],bg=bg.vec[i],lwd=2,cex=1.5)
  
}


text(x.left-range.x/8,y.bottom-range.y/5,"Hist. Salinity",cex=1,font=2,pos=4)
text(x.left+range.x/20-range.x/8,y.bottom-range.y/5-(1*range.y/16),"Low",cex=1,font=1,pos=4)
text(x.left+range.x/20-range.x/8,y.bottom-range.y/5-(2*range.y/16),"High",cex=1,font=1,pos=4)
points(rep(x.left+range.x/20-range.x/8,2),c(y.bottom-range.y/5-(1*range.y/16),y.bottom-range.y/5-(2*range.y/16)),
       col=c('#5D93B9','#BE4F48'),bg=c('#5D93B9','#BE4F48'),pch=21,lwd=1,cex=.95)

text(x.left+range.x/2.2,y.bottom-range.y/5,"Hist. Nitrogen",cex=1,font=2,pos=4)
text(x.left+range.x/20+range.x/2.2,y.bottom-range.y/5-(1*range.y/16),"Low",cex=1,font=1,pos=4)
text(x.left+range.x/20+range.x/2.2,y.bottom-range.y/5-(2*range.y/16),"High",cex=1,font=1,pos=4)
points(rep(x.left+range.x/20+range.x/2.2,2),c(y.bottom-range.y/5-(1*range.y/16),y.bottom-range.y/5-(2*range.y/16)),
       col=c('grey45'),bg=c('grey45','white'),pch=21,lwd=1,cex=.95)


dev.off()

