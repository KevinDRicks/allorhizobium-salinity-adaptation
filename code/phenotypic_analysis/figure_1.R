#############################################################################################################
#recreates figure 1 from Ricks et al
#############################################################################################################

#load packages
library(lme4)
library(lmerTest)

coef.table <- read.csv('/path/to/plant_growth_by_treatment.csv')

png('figure_3.png', width = 3.14961, height =5, units = "in", res = 1200, pointsize = 9)

#set y axis range between the two panels
fig.range <- c(1.5,2.65)

#excuse this ridiculous layout matrix. I wanted a specific ratio between the panels on each row
#and this was an easy way to hard code it in
layout(matrix(c(1,1,1,1,2,2,2,2,2,2,2,
                3,3,3,3,4,4,4,4,4,4,4),
              2, 11, byrow = T))
par(cex=1)

#############################################################################################################
#top panels show plant growth in low salinity environments
#subsets the coefficient data out to do draw onto figure
#############################################################################################################
par(mar = c(2.5, 2.25, .6, .15),oma=c(3,1,1,0),las=1)
plot(NA,xlim=c(0,4),ylim=c(fig.range),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
x.axis.pos <- c(x.left+range.x/4,x.right-range.x/4)



axis(at=x.axis.pos,side=1,labels=F,tck=.08,cex=1)
axis(side=2,tck=.08,labels=F,cex=1)
axis(side=2,labels=T,hadj=0.15,cex.axis=(.8),tck=.08)



par(xpd=NA)
text(x.left-range.x/15,y.top+range.y/15,'A',font=2,cex=1.25,pos=4)


text(x.left-range.x/2.55, y.bottom+range.y/2, 
     expression(bold("Plant growth (" * Delta ~ "mm"^2 ~ "/ day)")), 
     cex=1, font=2, srt=90)

jitter.wid <- range.x/40

sub.d <- coef.table[coef.table$Status=='Sterile' & coef.table$ContNitrogen=='Low' & coef.table$ContSaline=='Low',]
segments(x0=x.axis.pos[1],x1=x.axis.pos[1],
         y0=sub.d$Growth-1.96*sub.d$Growth_SE,y1=sub.d$Growth+1.96*sub.d$Growth_SE,col='grey30',lwd=2)  
points(x=x.axis.pos[1],y=sub.d$Growth,col='grey30',bg='grey30',pch=21,cex=1.25)

sub.d <- coef.table[coef.table$Status=='Ancestor' & coef.table$ContNitrogen=='Low' & coef.table$ContSaline=='Low',]
segments(x0=x.axis.pos[2],x1=x.axis.pos[2],
         y0=sub.d$Growth-1.96*sub.d$Growth_SE,y1=sub.d$Growth+1.96*sub.d$Growth_SE,col='grey30',lwd=2)  
points(x=x.axis.pos[2],y=sub.d$Growth,col='grey30',bg='grey30',pch=21,cex=1.25)


par(mar = c(2.5, .5, .6, .15))
plot(NA,xlim=c(0,4),ylim=c(fig.range),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)


x.axis.pos <- c(x.left+1*range.x/5,x.left+2*range.x/5,
                x.left+3*range.x/5,x.left+4*range.x/5)

axis(at=x.axis.pos,side=1,labels=F,tck=.025,cex=1)
axis(side=2,tck=.025,labels=F,cex=1)
text(x.left-range.x/30,y.top+range.y/35,'Low salinity environment',font=2,cex=1,pos=4)




#this controls the order of the subsequent loop 
salt.vec <- c('Low','Low','Low','Low')
nitH.vec <- c('Low','High','Low','High')
saltH.vec <- c('Low','Low','High','High')
status.vec <- c('Evolved','Evolved','Evolved','Evolved')
col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
lty.vec <- c(1,2,1,2)




for(i in 1:length(salt.vec)){
  working.coef <- coef.table[coef.table$ContSaline==salt.vec[i] & coef.table$NitrogenHistory==nitH.vec[i] &
                               coef.table$SalineHistory==saltH.vec[i] & coef.table$Status==status.vec[i],]
  segments(x0=x.axis.pos[i],x1=x.axis.pos[i],y0=working.coef$Growth-1.96*working.coef$Growth_SE,y1=working.coef$Growth+1.96*working.coef$Growth_SE,col=col.vec[i],lwd=2)  
  points(x=x.axis.pos[i],y=working.coef$Growth,col=col.vec[i],pch=21,bg=bg.vec[i],cex=1.25)
}

#############################################################################################################
#bottom panels show plant growth in high salinity environments
#subsets the coefficient data out to do draw onto figure
#identical in approach to the above panels
#############################################################################################################
par(mar = c(2.5, 2.25, .6, .15))

plot(NA,xlim=c(0,4),ylim=c(fig.range),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)
x.axis.pos <- c(x.left+range.x/4,x.right-range.x/4)


axis(at=x.axis.pos,side=1,labels=F,tck=.08,cex=1)
axis(side=2,tck=.08,labels=F,cex=1)
axis(side=2,labels=T,hadj=0.15,cex.axis=(.8),tck=.08)



par(xpd=NA)
text(x.left-range.x/2.55, y.bottom+range.y/2, 
     expression(bold("Plant growth (" * Delta ~ "mm"^2 ~ "/ day)")), 
     cex=1, font=2, srt=90)


sub.d <- coef.table[coef.table$Status=='Sterile' & coef.table$ContNitrogen=='Low' & coef.table$ContSaline=='High',]
segments(x0=x.axis.pos[1],x1=x.axis.pos[1],
         y0=sub.d$Growth-1.96*sub.d$Growth_SE,y1=sub.d$Growth+1.96*sub.d$Growth_SE,col='grey30',lwd=2)  
points(x=x.axis.pos[1],y=sub.d$Growth,col='grey30',bg='grey30',pch=21,cex=1.25)

sub.d <- coef.table[coef.table$Status=='Ancestor' & coef.table$ContNitrogen=='Low' & coef.table$ContSaline=='High',]
segments(x0=x.axis.pos[2],x1=x.axis.pos[2],
         y0=sub.d$Growth-1.96*sub.d$Growth_SE,y1=sub.d$Growth+1.96*sub.d$Growth_SE,col='grey30',lwd=2)  
points(x=x.axis.pos[2],y=sub.d$Growth,col='grey30',bg='grey30',pch=21,cex=1.25)

text(x.axis.pos,y.bottom-range.y/27,c('Sterile','Ancestor'),cex=.8)
text(x.left+range.x/2,y.bottom-range.y/9,'Controls',font=2,cex=1)

text(x.left-range.x/15,y.top+range.y/15,'B',font=2,cex=1.25,pos=4)


par(mar = c(2.5, .5, .6, .15))

plot(NA,xlim=c(0,4),ylim=c(fig.range),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
rng <- par("usr");y.top <- rng[4]; y.bottom <- rng[3];x.right <- rng[2];x.left <- rng[1];range.x <- abs(x.left-x.right);range.y <- abs(y.bottom-y.top)

x.axis.pos <- c(x.left+1*range.x/5,x.left+2*range.x/5,
                x.left+3*range.x/5,x.left+4*range.x/5)

axis(at=x.axis.pos,side=1,labels=F,tck=.025,cex=1)
axis(side=2,tck=.025,labels=F,cex=1)

text(x.left-range.x/30,y.top+range.y/35,'High salinity environment',font=2,cex=1,pos=4)

salt.vec <- c('High','High','High','High')
nitH.vec <- c('Low','High','Low','High')
saltH.vec <- c('Low','Low','High','High')
status.vec <- c('Evolved','Evolved','Evolved','Evolved')




col.vec <- c('#5D93B9','#5D93B9','#BE4F48','#BE4F48')
bg.vec <- c('#5D93B9','white','#BE4F48','white')
lty.vec <- c(1,2,1,2)

for(i in 1:length(salt.vec)){
  working.coef <- coef.table[coef.table$ContSaline==salt.vec[i] & coef.table$NitrogenHistory==nitH.vec[i] &
                               coef.table$SalineHistory==saltH.vec[i] & coef.table$Status==status.vec[i],]
  segments(x0=x.axis.pos[i],x1=x.axis.pos[i],y0=working.coef$Growth-1.96*working.coef$Growth_SE,y1=working.coef$Growth+1.96*working.coef$Growth_SE,col=col.vec[i],lwd=2)  
  points(x=x.axis.pos[i],y=working.coef$Growth,col=col.vec[i],pch=21,bg=bg.vec[i],cex=1.25)
}
text(x.left+range.x/2,y.bottom-range.y/9,'Evolved strains',font=2,cex=1)


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







