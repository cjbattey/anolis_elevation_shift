#anolis elevation range shift analyses
library(raster);library(ggplot2);library(plyr);library(magrittr);library(foreach);library(ggridges);library(gridExtra)
source("scripts/ggthemes.R")

#projections and extents for maps
ext.pr <- extent(-67.6,-65.2,17.8,18.7)
ext.yunque <- extent(825000,850000,2015000,2040000)
proj4.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4.utm <- "+proj=utm +zone=19 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

############ data loading & prep ##############
#load occurrence data (see load_occ_data.R for pre-processing of gbif data)
sight <- read.csv("data/occurrence/Anolis_sight_records_georeferenced.csv")
sight$ageClass <- sight$year>1977
sight$ageClass[sight$ageClass==T] <- 4
sight$ageClass[sight$ageClass==F] <- 2
anolis <- read.csv("data/occurrence/anolis_data.csv")
anolis <- subset(anolis,ageClass %in% c(2,4))
anolis <- anolis[!duplicated(anolis$gbifID),] #drop duplicates (strangely, these were in the file downloaded from GBIF)
comb <- rbind.fill(anolis,sight)
comb$species <- gsub("Anolis ","",comb$species)
anolis$ageClass2 <- revalue(factor(anolis$ageClass),c("2"="1952-1977","4"="1991-2015"))
comb$ageClass2 <- revalue(factor(comb$ageClass),c("2"="1952-1977","4"="1991-2015"))
comb <- subset(comb,species %in% c("cristatellus","pulchellus","stratulus","evermanni","gundlachi","krugi"))
comb$species <- factor(comb$species)
comb$species <- factor(comb$species,levels(comb$species)[c(1,5,6,2,3,4)])

#load Helmer 2008 forest cover data
source("scripts/load_helmer_data.R")

#load GIS data for plots
ext.pr <- extent(-67.6,-65.2,17.8,18.7)
ext.yunque <- extent(825000,850000,2015000,2040000)
proj4.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4.utm <- "+proj=utm +zone=19 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
alt <- raster("data/elevation/alt_30s_UTM19N.tif")
alt1s <- raster("data/elevation/alt_1s_UTM19N.tif")
coasts <- shapefile("data/ne_10m_coastline/ne_10m_coastline.shp") %>% 
            crop(ext.pr)
coasts.utm <- spTransform(coasts,proj4.utm)

#cluster localities
coords <- subset(comb,!is.na(lat)) 
coords <- ddply(coords,.(long,lat,ageClass),summarize,nobs=length(day)) 
pts <- SpatialPoints(coords=data.frame(coords$long,coords$lat),proj4string=crs(proj4.wgs))
dist <- spDists(pts,pts) %>% as.dist()                                      
fit <- hclust(dist,method="complete")                                        
clusters <- cutree(fit,h=1)                                                 
coords$locality_group <- clusters                                            
obs.per.group <- ddply(coords,.(locality_group,ageClass),summarize,obs.per.group=sum(nobs))
coords <- join(coords,obs.per.group,by=c("locality_group","ageClass"))
anolis <- join(anolis,coords,by=c("long","lat","ageClass"))
comb <- join(comb,coords,by=c("long","lat","ageClass"))
nlevels(factor(comb$locality_group))
loc_group_alt <- ddply(comb,.(locality_group),summarize,loc_group_alt=mean(na.omit(alt)))
comb <- join(comb,loc_group_alt,by="locality_group")

#data frame including "all species" 
comb_w_all <- comb
comb_w_all$species <- "All Species"
comb_w_all <- rbind(comb,comb_w_all)
anolis_w_all <- anolis
anolis_w_all$species <- "All Species"
anolis_w_all <- rbind(anolis,anolis_w_all)

#summarize by locality group for per-locality analyses & plots
locs <- ddply(comb_w_all,.(locality_group,ageClass2,ageClass,species),summarize,n=length(alt),alt=mean(alt),lat=mean(lat),long=mean(long))

############# end data loading and prep ################
########################################################




########################################################
############# descriptive figures and stats ############

#sample map
map <- fortify(coasts)
tmp <- ddply(anolis,.(locality_group,ageClass2),summarize,long=mean(long),lat=mean(lat),alt=mean(alt),n=length(day))
ggplot()+theme_blank()+coord_map()+
  theme(legend.position = "right")+
  facet_grid(~ageClass2)+
  scale_color_gradient(low=grey(0.9),high=grey(0.1))+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=tmp,aes(x=long,y=lat,size=n,col=alt))

#counts of collections over time
tmp <- ddply(comb,.(year),summarize,n=length(day)) %>% subset(year > 1950)
ggplot(data=tmp,aes(x=year,y=n))+
  scale_fill_brewer(type = "div",palette = "Set1")+
  theme_minimal()+theme(legend.position = "right")+
  geom_bar(stat="identity",fill="black")

#change in forest cover by age Class
age1 <- reclassify(age1,t(matrix(c(0,NA))))
age2 <- reclassify(age2,t(matrix(c(0,NA))))
age3 <- reclassify(age3,t(matrix(c(0,NA))))
age4 <- reclassify(age4,t(matrix(c(0,NA))))
alt2 <- projectRaster(alt,age1)
pdf(file="forest_map_wide.pdf",width=6.5,height=1.5)
par(mfrow=c(1,4),oma = c(1,1,0,0) + 0.1,mar = c(1,1,1,1) + 0.1)
plot(age1,main="1935-1951",axes=F,legend=F,col="forestgreen")
plot(coasts.utm,lwd=.5,add=T)
plot(age2==1,main="1952-1977",axes=F,legend=F,col="forestgreen")
plot(coasts.utm,lwd=.5,add=T)
plot(age3==1,main="1978-1990",axes=F,legend=F,col="forestgreen")
plot(coasts.utm,lwd=.5,add=T)
plot(age4==1,main="1991-2000",axes=F,legend=F,col="forestgreen")
plot(coasts.utm,lwd=.5,add=T)
par(mfrow=c(1,1))
dev.off()

#elevation density of forests
a<-as.data.frame(mask(alt2,age1),na.rm=T)
a$age <- "1931-1951"
b<-as.data.frame(mask(alt2,age2),na.rm=T)
b$age <- "1952-1977"
c<-as.data.frame(mask(alt2,age3),na.rm=T)   
c$age <- "1978-1990"
d<-as.data.frame(mask(alt2,age4),na.rm=T)
d$age <- "1991-2008"
tmp <- rbind(a,b,c,d)

tmp <- data.frame(alt=na.omit(as.data.frame(mask(alt2,age1))),age="1935-1951") %>% 
        rbind(data.frame(alt=na.omit(as.data.frame(mask(alt2,age2))),age="1952-1977")) %>% 
          rbind(data.frame(alt=na.omit(as.data.frame(mask(alt2,age3))),age="1978-1990")) %>% 
            rbind(data.frame(alt=na.omit(as.data.frame(mask(alt2,age4))),age="1991-2000"))
tmp$age <- factor(tmp$age,levels(tmp$age)[c(1,2,3,4)])
cellConversion <- function(x){ 
  x*0.0009 
}
pdf("forest_histogram.pdf",width=6,height=2)
ggplot(data=tmp,aes(x=alt_30s_UTM19N))+
  scale_fill_manual(values = c("grey10","grey35","grey60","grey85"),name="Time Period")+
  facet_wrap(~age,nrow=1)+coord_flip()+
  theme_minimal()+theme(strip.background = element_blank())+
  xlab("Elevation")+ylab(("Forest Area "~(km^2)))+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position = "right")+
  scale_y_continuous(labels=cellConversion)+
  geom_histogram(bins=100,col=NA,fill="forestgreen")
dev.off()

#downsampled forest age map in ggplot (full res freezes grid)
tmp.r <- projectRaster(age1,crs=crs(proj4.wgs),res=.5/60)
tmp.p <- rasterToPolygons(tmp.r,dissolve = T)
f1 <- fortify(tmp.p)
f1$ageClass2 <- "1931-1951"
tmp.r <- projectRaster(age2,crs=crs(proj4.wgs),res=.5/60)
tmp.p <- rasterToPolygons(tmp.r,dissolve = T)
f2 <- fortify(tmp.p)
f2$ageClass2 <- "1952-1977"
tmp.r <- projectRaster(age3,crs=crs(proj4.wgs),res=.5/60)
tmp.p <- rasterToPolygons(tmp.r,dissolve = T)
f3 <- fortify(tmp.p)
f3$ageClass2 <- "1978-1990"
tmp.r <- projectRaster(age4,crs=crs(proj4.wgs),res=.5/60)
tmp.p <- rasterToPolygons(tmp.r,dissolve = T)
f4 <- fortify(tmp.p)
f4$ageClass2 <- "1991-2008"

forest <- rbind(f1,f2,f3,f4)
forest$ageClass2 <- factor(forest$ageClass2,levels=levels(factor(forest$ageClass2))[c(4,3,2,1)])

ggplot()+
  theme_minimal()+coord_map()+
  theme(strip.background = element_rect(fill="white",color="white"),
        axis.text=element_blank(),axis.ticks = element_blank(),legend.position="right")+
  ylab(" ")+xlab("")+
  scale_fill_manual(values=c("darkgreen","chartreuse3","gold","orangered"))+
  geom_path(data=fortify(coasts),aes(x=long,y=lat,group=group),col="black")+
  geom_tile(data=forest,aes(x=long,y=lat,group=group,fill=ageClass2))

plot(coasts.utm)
  plot(age4,add=T,col="orangered",legend=F)
  plot(age3,add=T,col="gold",legend=F)
  plot(age2,add=T,col="chartreuse3",legend=F)
  plot(age1,add=T,col="darkgreen",legend=F)
  rect(875000,2075000,890000,2090000,col="orangered")
  rect(875000,2055000,890000,2070000,col="gold")
  rect(875000,2035000,890000,2050000,col="chartreuse3")
  rect(875000,2015000,890000,2030000,col="darkgreen")
  

  
  

grid.arrange(forest_maps,forest_histograms,layout_matrix=matrix(c(1,1,1,1,
                                                                  1,1,1,1,
                                                                  2,2,2,2,
                                                                  2,2,2,2),byrow=T,nrow=4))


############# end descriptive figures and stats ############
#############################################################

  
  
  
#############################################################
############# elevation shift stats and figures #############
#violin plot of elevation by time period
elevation_plot <- ggplot(data=comb,aes(x=factor(ageClass2),y=alt))+
  facet_grid(.~species)+
  theme_minimal()+xlab("")+ylab("")+
  theme(strip.text.x=element_text(size=8),strip.background = element_rect(color="white",fill="white"), 
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        plot.title = element_text(hjust=.75))+
  #ggtitle("Anolis Specimen Elevation\n1952-1977 vs 1991-2015")+
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5))
  
all_plot <- ggplot(data=subset(comb_w_all,species=="All Species"),aes(x=factor(ageClass2),y=alt))+
  facet_grid(.~species)+
  theme_minimal()+xlab("")+ylab("")+
  theme(strip.text.x=element_text(size=8),strip.background = element_rect(color="white",fill="white"), 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5))

locs_plot <- ggplot(data=subset(locs,species != "All Species"),aes(x=factor(ageClass2),y=alt))+
  facet_grid(.~species)+
  theme_minimal()+xlab("")+ylab("")+
  theme(strip.text.x=element_text(size=8),strip.background = element_rect(color="white",fill="white"), 
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5))

all_locs_plot <- ggplot(data=subset(locs,species=="All Species"),aes(x=factor(ageClass2),y=alt))+
  facet_grid(.~species)+
  theme_minimal()+xlab("")+ylab("")+
  theme(strip.text.x=element_text(size=8),strip.background = element_rect(color="white",fill="white"), 
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  geom_point(size=.2,col="grey",position="jitter")+
  geom_violin(draw_quantiles = c(0.5))

layout <- matrix(c(1,1,1,1,2,
                   3,3,3,3,4),nrow=2,byrow=T)
grid.arrange(elevation_plot,all_plot,locs_plot,all_locs_plot,ncol=2,nrow=2,layout_matrix=layout,left="Elevation (m)")


#test for elevation shift per specimen
# mu <- wilcox.test(subset(comb_w_all,species=="All Species" & ageClass==4)$alt,
#                   subset(comb_w_all,species=="All Species" & ageClass==2)$alt,
#                   conf.int=T)$estimate
tmp <- ddply(comb,.(species),function (i){
  mu <- wilcox.test(subset(comb,species != i$species[1] & ageClass==4)$alt,
                    subset(comb,species != i$species[1] & ageClass==2)$alt,
                    conf.int=T)$estimate
  t <- wilcox.test(i$alt[i$ageClass==4],i$alt[i$ageClass==2],conf.int=T,mu=mu)
  n.pre=nrow(subset(i,ageClass==2))
  pre.median=round(median(i$alt[i$ageClass==2]),2)
  pre.min=round(min(i$alt[i$ageClass==2]),2)
  pre.max=round(max(i$alt[i$ageClass==2]),2)
  n.post=nrow(subset(i,ageClass==4))
  post.median=round(median(i$alt[i$ageClass==4]),2)
  post.min=round(min(i$alt[i$ageClass==4]),2)
  post.max=round(max(i$alt[i$ageClass==4]),2)
  shift=round(t$estimate,2)
  #shift=round(median(i$alt[i$ageClass==4]) - median(i$alt[i$ageClass==2]),2)
  p=t$p.value
  pre <- paste0(pre.median," (",pre.min," - ",pre.max,") ")
  post <- paste0(post.median," (",post.min," - ",post.max,") ")
  sum <- c(n.pre,pre,n.post,post,shift,p)
  names(sum) <- c("n.pre","pre","n.post","post","shift","p")
  sum
})

#test for elevation shift per locality
# mu <- wilcox.test(subset(locs,species=="All Species" & ageClass==4)$alt,
#                   subset(locs,species=="All Species" & ageClass==2)$alt,
#                   conf.int=T)$estimate
ddply(locs,.(species),function (i){
  mu <- wilcox.test(subset(locs,species != i$species[1] & ageClass==4)$alt,
                    subset(locs,species != i$species[1] & ageClass==2)$alt,
                    conf.int=T)$estimate
  t <- wilcox.test(i$alt[i$ageClass==4],i$alt[i$ageClass==2],conf.int=T,mu=mu)
  n.pre=nrow(subset(i,ageClass==2))
  pre.median=round(median(i$alt[i$ageClass==2]),2)
  pre.min=round(min(i$alt[i$ageClass==2]),2)
  pre.max=round(max(i$alt[i$ageClass==2]),2)
  n.post=nrow(subset(i,ageClass==4))
  post.median=round(median(i$alt[i$ageClass==4]),2)
  post.min=round(min(i$alt[i$ageClass==4]),2)
  post.max=round(max(i$alt[i$ageClass==4]),2)
  shift=round(t$estimate,2)
  #shift=round(median(i$alt[i$ageClass==4]) - median(i$alt[i$ageClass==2]),2)
  p=t$p.value
  pre <- paste0(pre.median," (",pre.min," - ",pre.max,") ")
  post <- paste0(post.median," (",post.min," - ",post.max,") ")
  sum <- c(n.pre,pre,n.post,post,shift,p)
  names(sum) <- c("n.pre","pre","n.post","post","shift","p")
  sum
})

#low-altitude frequency shift analysis
low <- subset(comb,alt<=250)
n.sp.pre <- nrow(subset(low,species=="gundlachi" & ageClass==2))
n.sp.post <- nrow(subset(low,species=="gundlachi" & ageClass==4))
n.all.pre <- nrow(subset(low,ageClass==2))
n.all.post <- nrow(subset(low,ageClass==4))
chimatrix <- matrix(c(n.sp.pre,n.sp.post,n.all.pre,n.all.post),byrow = T, nrow = 2)
mcnemar.test(chimatrix)
chisq.test(chimatrix)

anolischi2fun <- function(df,i){
  n.sp.pre <- nrow(subset(df,species==i & ageClass==2))
  n.sp.post <- nrow(subset(df,species==i & ageClass==4))
  n.all.pre <- nrow(subset(df,species !=i & ageClass==2))
  n.all.post <- nrow(subset(df,species !=i & ageClass==4))
  chimatrix <- matrix(c(n.sp.pre,n.sp.post,n.all.pre,n.all.post),byrow = T, nrow = 2)
  chi <- chisq.test(chimatrix)
  c(i,round(n.sp.pre/n.all.pre,3),round(n.sp.post/n.all.post,3),round(chi$p.value,5))
}
chi2table <- data.frame(foreach(i=levels(factor(low$species)),.combine=rbind) %do% anolischi2fun(low,i),row.names=NULL)
names(chi2table) <- c("Species","relAbund.pre","relAbund.post" ,"p")
chi2table

#species richness by altitude bin & age class
#shannon_diversity() takes a data frame including a column named "species" 
#and returns the shannon index as a vector 
shannon_diversity <- function(x) { 
  plnp <- c()
  n <- nrow(x)
  for(i in unique(x$species)){
    sp <- subset(x,species==i)
    p <- nrow(sp)/n
    plnp <- append(plnp,p*log(p))
  }
  return(c(div=-sum(plnp)))
}

#drop sites with only one species 
div_data <- ddply(comb,.(locality_group,ageClass),function(e){
  if(length(unique(e$species))>1){
    e
  }
})
div_data <- subset(div_data,!is.na(locality_group))
div_low <- subset(div_data,loc_group_alt<250) %>% ddply(c("locality_group","ageClass2"),function(e) shannon_diversity(e))
div_low$altbin <- "0-250 m"
div_mid <- subset(div_data,loc_group_alt>=250 & loc_group_alt<500) %>% ddply(c("locality_group","ageClass2"),function(e) shannon_diversity(e))
div_mid$altbin <- "251-500 m"
div_high <- subset(div_data,loc_group_alt>=500 & loc_group_alt<750) %>% ddply(c("locality_group","ageClass2"),function(e) shannon_diversity(e))
div_high$altbin <- "501-750 m"
div_vhigh <- subset(div_data,loc_group_alt>=750) %>% ddply(c("locality_group","ageClass2"),function(e) shannon_diversity(e))
div_vhigh$altbin <- "750-1150 m"
div <- rbind(div_low,div_mid,div_high,div_vhigh)

#wilcox test on diversity across time periods in altitude bins
div_wilcox_table <- ddply(div,.(altbin),function(e){
  w <- wilcox.test(e$div[e$ageClass2=="1991-2015"],e$div[e$ageClass2=="1952-1977"],conf.int = T)
  n1 <- nrow(subset(e,ageClass2=="1952-1977"))
  n2 <- nrow(subset(e,ageClass2=="1991-2015"))
  p <- w$p.value %>% round(3)
  shift <- w$estimate
  conf <- c(w$conf.int[1],w$conf.int[2])
  pre_med <- median(subset(e,ageClass2=="1952-1977")$div)
  post_med <- median(subset(e,ageClass2=="1991-2015")$div)
  c(p,n1,n2,shift,conf,pre_med,post_med)
})
colnames(div_wilcox_table) <- c("Altitude Bin","P","N1","N2","Shift","CI_0.05","CI_0.95","pre_med","post_med")
div_wilcox_table

bin_plot <- ggplot(div,aes(x=ageClass2,y=div))+facet_grid(~altbin)+
  theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1),
                        strip.background = element_rect(color="white",fill="white"),
                        text=element_text(size=10),
                        axis.title=element_text(size=8))+
  scale_y_continuous(breaks=c(0,1))+
  ylab("")+xlab("")+
  geom_violin(draw_quantiles = 0.5,fill=NA)+
  annotate(geom="text",label=c("*"," "," "," "),x=1.5,y=1.65,size=5)

#regression of locality species richness by altitude, split by age class
locs <- ddply(comb,.(locality_group,ageClass,ageClass2),summarize,alt=mean(na.omit(alt)),div=shannon_diversity(data.frame(species)))
locs <- subset(locs,div>0)
locs2 <- subset(locs,ageClass==2)
locs4 <- subset(locs,ageClass==4)
fit2 <- nls(formula=div~a/(1+exp(-b*(alt-c))),data=locs2,start=list(a=1,b=.5,c=5))
fit4 <- nls(formula=div~a/(1+exp(-b*(alt-c))),data=locs4,start=list(a=1,b=.5,c=5))
CI2 <- confint2(fit2,level=.8)
CI2_low <- sapply(locs2$alt,function(e) CI2[1,1]/(1+exp(-CI2[2,1]*(e-CI2[3,1]))))
CI2_high <- sapply(locs2$alt,function(e) CI2[1,2]/(1+exp(-CI2[2,2]*(e-CI2[3,2]))))
CI4 <- confint2(fit4,level=.8)
CI4_low <- sapply(locs4$alt,function(e) CI4[1,1]/(1+exp(-CI4[2,1]*(e-CI4[3,1]))))
CI4_high <- sapply(locs4$alt,function(e) CI4[1,2]/(1+exp(-CI4[2,2]*(e-CI4[3,2]))))
curves <- data.frame(ageClass2=c(rep("1952-1977",nrow(locs2)),rep("1991-2015",nrow(locs4))),
                     alt=c(locs2$alt,locs4$alt),
                     div=c(predict(fit2),predict(fit4)),
                     CI_low=c(CI2_low,CI4_low),
                     CI_high=c(CI2_high,CI4_high))

ggplot(data=locs,aes(x=alt,y=div))+facet_wrap(~ageClass2)+
  theme_minimal()+theme(strip.background = element_blank())+
  geom_point()+
  geom_line(data=curves)+
  geom_line(data=curves,aes(y=CI_low),col="blue")+
  geom_line(data=curves,aes(y=CI_high),col="red")

lm(locs$div[locs$ageClass==2]~locs$alt[locs$ageClass==2]) %>% summary() #p<<0.01, R2=0.1833, slope=0.212 sp/100M #p=0.048, R2=0.063, slope=
lm(locs$div[locs$ageClass==4]~locs$alt[locs$ageClass==4]) %>% summary() #p<<0.01, R2=0.0697, slope=0.149 sp/100M #p=0.016, R2=0.036, slope=



lm_plot<-ggplot(locs,aes(x=alt,y=div,shape=ageClass2,col=ageClass2))+
  theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.background = element_rect(color="white",fill="white"),
        legend.position=c(.82,.15),legend.background = element_rect(fill=NA),
        text=element_text(size=10),
        axis.title=element_text(size=8))+
  facet_grid(~ageClass2)+
  geom_point(show.legend=F)+
  scale_shape_discrete(solid = F)+
  scale_x_continuous(breaks=c(300,900))+
  scale_y_continuous(breaks=c(0,1))+
  scale_color_manual(values=c("grey30","grey65"),name="Years")+
  xlab("Elevation (m)")+ylab("")+
  geom_smooth(method="lm",fill=NA,show.legend=F)

  
pdf("~/Dropbox/anolis/figures/final/figS3_div_plot.pdf",width=3,height=4.5) 
grid.arrange(bin_plot,lm_plot,ncol=1)
dev.off()

############# end elevation shift stats and figures #############
#################################################################


############# resurveyed sites analysis #########################
#################################################################
#resurveyed sites (237==ranger house, 242==el verde station)
rsv <- comb
groups <- ddply(rsv,.(locality_group),summarize,nAgeClass=length(unique(ageClass))) %>% 
  subset(nAgeClass>1) %>%
  .$locality_group  # find groups w multiple ageClasses
rsv <- subset(rsv,locality_group %in% groups)
groups <- ddply(rsv,.(locality_group,ageClass),summarize,n=length(day)) %>%   #find groups w min specimens > 100 per ageClass
  ddply(.(locality_group),summarize,min=min(n)) %>%
  subset(min>=100) %>% 
  .$locality_group
rsv <- subset(rsv,locality_group %in% groups) 


############ forest regrowth for a 1-km radius of "ranger house"/ campamento eliza colberg
ranger <- anolis[grep("Ranger House",anolis$locality),]
sp <- SpatialPoints(data.frame(ranger$long,ranger$lat),proj4string = crs(alt))
sp <- spTransform(sp,crs(alt1s))
buffer <- buffer(sp,1000)
plot(coasts.utm)
plot(age4,add=T,col="yellow2")
plot(age3,add=T,col="gold")
plot(age2,add=T,col="chartreuse3")
plot(age1,add=T,col="darkgreen")
plot(buffer,add=T)
plot(coasts.utm,add=T)
plot(localities,add=T)


plotdata <- ddply(comb,.(locality_group,ageClass2,ageClass,species),summarize,n=length(day),long=mean(long),lat=mean(lat))

ggplot()+facet_grid(species~ageClass2)+
  coord_map()+
  theme_minimal()+
  theme(strip.background = element_rect(fill="white",color="white"),
        axis.text=element_blank(),axis.ticks = element_blank(),legend.position="right")+
  ylab("")+xlab("")+
  scale_size_continuous(name="Number of\nObservations",breaks=c(10,50,100,200,300),range=c(1,10))+
  geom_path(data=fortify(coasts),aes(x=long,y=lat,group=group))+
  geom_polygon(data=subset(forest,ageClass2 %in% c("1952-1977","1991-2015")),aes(x=long,y=lat,group=group),fill="grey")+
  geom_point(data=subset(forest,ageClass2 %in% c("1952-1977","1991-2015")),aes(x=long,y=lat),shape=1)


