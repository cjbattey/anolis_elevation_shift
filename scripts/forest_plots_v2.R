library(raster);library(rgeos);library(magrittr)
library(foreach);library(doMC)
registerDoMC(cores=12)
source("scripts/load_helmer_data.R")
age1 <- reclassify(age1,t(matrix(c(0,NA))))
age2 <- reclassify(age2,t(matrix(c(0,NA))))
age3 <- reclassify(age3,t(matrix(c(0,NA))))
age4 <- reclassify(age4,t(matrix(c(0,NA))))
alt2 <- projectRaster(alt,age1)
alt2[alt2<=0] <- NA
land <- as.data.frame(alt2)
land[land>0] <- 1
la <- sum(land,na.rm=T)

totals <- c()
for(i in list(age1,age2,age3,age4)){
  totals <- append(totals,sum(as.data.frame(i),na.rm=T))
}

#total elevation over time
plot(x=c(1943,1965,1984,1995),y=totals/la,ylim=c(0,.8),
     type="l",lwd=1.5,tck=-.05,bty="n",
     mgp=c(1.7,.2,0),cex.axis=0.7,cex.lab=0.75,
     xlab="Year",ylab="Proportion\nTotal Area",
)
points(x=c(1943,1965,1984,1995),y=totals/la,pch=16,cex=0.8)

#proportion land area forested, by 20 m elevation bin
alts <- seq(1,1330,10) 
forest_area <- c()
land_area <- c()
for(i in alts){
  print(i/max(alts))
  mask <- alt2>i
  mask2 <- alt2>(i+10)
  mask[mask2==1] <- 0
  mask[alt2<=0] <- 0
  mask@data@values <- as.numeric(mask@data@values)
  altarea <- sum(as.data.frame(mask),na.rm=T)
  for(j in list(age1,age2,age3,age4)){
    a <- j
    a[mask==0] <- 0
    a[is.na(mask)] <- 0
    a[alt2<0] <- 0
    forest_area <- append(forest_area,sum(as.data.frame(a),na.rm=T))
    land_area <- append(land_area,altarea)
  }
}
df <- data.frame(alt=c(sapply(alts,function(e) rep(e,4))),
                 forest_area,
                 land_area)
df$prop_forest <- df$forest_area/df$land_area
df$age <- rep(c(1943,1965,1984,1995),length(alts))


pdf("~/Dropbox/anolis/figures/forest_maps_revised.pdf",width=5,height=2)
par(mar=c(0,0,1,0),mfrow=c(2,2),bty='n',yaxt='n',xaxt='n',cex.main=0.8,font.main=1)
plot(age1,col="darkgreen",interpolate=F,legend=F,main="1935-1951",line=0)
lines(coasts,lwd=0.5)
plot(age2,col="darkgreen",interpolate=F,legend=F,main="1952-1977",line=0)
lines(coasts,lwd=0.5)
plot(age3,col="darkgreen",interpolate=F,legend=F,main="1978-1990",line=0)
lines(coasts,lwd=0.5)
plot(age4,col="darkgreen",interpolate=F,legend=F,main="1991-2000",line=0)
lines(coasts,lwd=0.5)
dev.off()

rel_incs <- ddply(df,.(alt),function(e) {
  absinc <- e$forest_area[e$age==1995]-e$forest_area[e$age==1943]
  relinc <- e$forest_area[e$age==1995]/e$forest_area[e$age==1943]
  return(c(absinc,relinc,e$prop_forest[e$age==1995][1],e$prop_forest[e$age==1943][1]))})

rel_incs$V1_km <- rel_incs$V1 * (30/1e3)^2

pdf("~/Dropbox/anolis/figures/forest_lineplots.pdf",width=6.5,height=2.2,useDingbats = F)

par(mfrow=c(1,2),mar=c(2.5,1.5,2.5,1.5),cex.main=0.77,tck=-0.03,cex.axis=0.7,cex.lab=0.75,font.main=1,
    mgp=c(1,.2,0),bty='n')
plot(x=rel_incs$alt[rel_incs$alt>10],y=rel_incs$V1_km[rel_incs$alt>10],type="l",xlab="Elevation (m)",ylab="",bty="n",
     main=expression(atop(Increase~"in"~Forested~Land~1935-2000,(km^2/10~m~elevation~band))))

plot(x=NA,y=NA,ylim=c(0,1),xlim=c(0,1200),
     ylab="",
     #ylab="Proportion\nTotal Area",
     xlab="Elevation (m)",
     main=expression(atop("Proportion Forested Land","(10 m elevation bands)")))
j <- 1
for(i in rev(unique(df$age))){
  pal <- grey.colors(4)
  #pal <- c("steelblue4","skyblue","khaki","orangered")
  d <- df[df$age==i,]
  d <- d[d$alt>20,]
  lines(x=d$alt,y=d$prop_forest,type='l')
  polygon(c(0,d$alt,max(d$alt)),c(0,d$prop_forest,0),col=pal[j])
  legend(x ="topleft",col=pal,legend = rev(unique(df$age)),
         bty = 'n',pch=16,pt.cex = 1.3,cex=0.7,ncol=2)
  j <- j+1
}
dev.off()

