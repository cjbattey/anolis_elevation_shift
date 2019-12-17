#anolis elevation range shift analyses
library(raster);library(ggplot2);library(plyr);library(magrittr)
library(foreach);library(ggridges);library(gridExtra);library(pwr)
#setwd("~/Dropbox/anolis/anolis_elevation_shift/")
source("scripts/ggthemes.R")

#projections and extents for maps
ext.pr <- extent(-67.6,-65.2,17.8,18.7)
ext.yunque <- extent(825000,850000,2015000,2040000)
proj4.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4.utm <- "+proj=utm +zone=19 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

############ data loading & prep ##############
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

#cluster localities
coords <- subset(comb,!is.na(lat)) 
coords <- ddply(coords,.(long,lat,ageClass),summarize,nobs=length(day)) 
pts <- SpatialPoints(coords=data.frame(coords$long,coords$lat),proj4string=crs(proj4.wgs))
dist <- spDists(pts,pts,longlat = T) %>% as.dist()                                      
fit <- hclust(dist,method="complete")                                        
clusters <- cutree(fit,h=.5)                                                 
coords$locality_group <- clusters                                            
obs.per.group <- ddply(coords,.(locality_group,ageClass),summarize,obs.per.group=sum(nobs))
coords <- join(coords,obs.per.group,by=c("locality_group","ageClass"))
anolis <- join(anolis,coords,by=c("long","lat","ageClass"))
comb <- join(comb,coords,by=c("long","lat","ageClass"))
nlevels(factor(comb$locality_group))
loc_group_alt <- ddply(comb,.(locality_group),summarize,loc_group_alt=mean(na.omit(alt)))
comb <- join(comb,loc_group_alt,by="locality_group")
comb <- comb[!duplicated(names(comb))]
anolis <- anolis[!duplicated(names(anolis))]

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
pdf("~/Dropbox/anolis/figures/final/collections_by_year.pdf",width=5,height=1.25)
ggplot(data=tmp,aes(x=year,y=n))+
  ylab("Specimens\nCollected")+xlab("Year")+
  scale_fill_brewer(type = "div",palette = "Set1")+
  theme_classic()+theme(axis.title=element_text(size=8),axis.text=element_text(size=7))+
  geom_bar(stat="identity",fill="black")
dev.off()

#how many of the specimens were collected by the authors? 
sum(grepl("Otero|Gorman|Hertz|Lister|Garcia|Burrowes|Perez|Huey",anolis$recordedBy))/nrow(anolis) #34% of all records(!)
sum(grepl("Otero|Gorman|Hertz|Lister|Garcia|Burrowes|Perez|Huey",subset(anolis,ageClass==2)$recordedBy))/nrow(subset(anolis,ageClass==2)) #52% of records in ageClass 2

#change in forest cover by age Class
age1 <- reclassify(age1,t(matrix(c(0,NA))))
age2 <- reclassify(age2,t(matrix(c(0,NA))))
age3 <- reclassify(age3,t(matrix(c(0,NA))))
age4 <- reclassify(age4,t(matrix(c(0,NA))))
alt2 <- projectRaster(alt,age1)
alt2 <- reclassify(alt2,rcl=matrix(c(-100,5,NA),byrow=T,nrow=1))
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
a <-hist(mask(alt2,age1),breaks=seq(0,1300,10),plot=F)
a <- data.frame(alt=a$mids,cellcount=a$counts,age="1935-1951",mask="forest")
b <-hist(mask(alt2,age2),breaks=seq(0,1300,10),plot=F)
b <- data.frame(alt=b$mids,cellcount=b$counts,age="1952-1977",mask="forest")
c <-hist(mask(alt2,age3),breaks=seq(0,1300,10),plot=F)
c <- data.frame(alt=c$mids,cellcount=c$counts,age="1978-1990",mask="forest")
d <-hist(mask(alt2,age4),breaks=seq(0,1300,10),plot=F)
d <- data.frame(alt=d$mids,cellcount=d$counts,age="1991-2016",mask="forest")
forest_hist <- rbind(a,b,c,d)
forest_hist$area <- forest_hist$cellcount*0.0009 #convert count of 30m2 grid cells to km2

t <- hist(alt2,breaks=seq(0,1300,10),plot=F) #total area
e <- data.frame(alt=t$mids,cellcount=t$counts,age="1935-1951")
f <- data.frame(alt=t$mids,cellcount=t$counts,age="1952-1977")
g <- data.frame(alt=t$mids,cellcount=t$counts,age="1978-1990")
h <- data.frame(alt=t$mids,cellcount=t$counts,age="1991-2016")
total_hist <- rbind(e,f,g,h)
total_hist$area <- total_hist$cellcount*0.0009

pdf("~/Dropbox/anolis/figures/forest_histograms.pdf",width=6,height=2)
ggplot(data=forest_hist,aes(x=alt,y=area))+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.background = element_blank())+
  facet_grid(~age)+coord_flip()+
  scale_y_continuous(breaks=c(0,100,200))+
  geom_bar(data=forest_hist,stat="identity",width = 10,fill="black")+
  geom_line(data=subset(total_hist,alt>50),lwd=0.5,col="grey")
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
library(ggforce)
comb_w_all <- comb_w_all[!duplicated(names(comb_w_all))]

#violin plot of elevation by time period
elevation_plot <- ggplot(data=comb_w_all,aes(x=factor(ageClass2),y=alt))+
  facet_grid(.~species)+
  theme_classic()+xlab("")+ylab("Elevation (M)")+
  theme(strip.text.x=element_text(size=8),
        strip.background = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=7,angle=45,hjust=1),
        plot.title = element_text(hjust=.75))+
  geom_sina(size=0.7,shape=1,stroke=0.34,alpha=0.6)
  # geom_point(size=.2,col="grey",position="jitter")+
  # geom_violin(draw_quantiles = c(0.5))
  
pdf("~/Dropbox/anolis/manuscript/procB/fig1_revised2.pdf",width=6,height=2.75,useDingbats = F)
print(elevation_plot)
dev.off()

#test for elevation shift per specimen
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
  p=p.adjust(t$p.value,method="fdr",n=6)
  pre <- paste0(pre.median," (",pre.min," - ",pre.max,") ")
  post <- paste0(post.median," (",post.min," - ",post.max,") ")
  sum <- c(n.pre,pre,n.post,post,shift,p)
  names(sum) <- c("n.pre","pre","n.post","post","shift","p")
  sum
})

#test for elevation shift per locality
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
  p=p.adjust(t$p.value,method="holm",n=6)
  pre <- paste0(pre.median," (",pre.min," - ",pre.max,") ")
  post <- paste0(post.median," (",post.min," - ",post.max,") ")
  sum <- c(n.pre,pre,n.post,post,shift,p)
  names(sum) <- c("n.pre","pre","n.post","post","shift","p")
  sum
})

#power analysis for A. gundlachi in per-locality tests (assuming a t test)
n1 <- nrow(subset(locs,ageClass==2 & species=="gundlachi"))
n2 <- nrow(subset(locs,ageClass==4 & species=="gundlachi"))
d <- (mean(subset(locs,ageClass==2 & species=="gundlachi")$alt)-
        mean(subset(locs,ageClass==4 & species=="gundlachi")$alt))/sd(subset(locs,species=="gundlachi")$alt)
pwr.t2n.test(n1=n1,n2=n2,d=d)

#power analysis for A. gundlachi in full dataset
n1 <- nrow(subset(comb,ageClass==2 & species=="gundlachi"))
n2 <- nrow(subset(comb,ageClass==4 & species=="gundlachi"))
d <- (mean(subset(comb,ageClass==2 & species=="gundlachi")$alt)-
        mean(subset(comb,ageClass==4 & species=="gundlachi")$alt))/sd(subset(comb,species=="gundlachi")$alt)
pwr.t2n.test(n1=n1,n2=n2,d=d)

#low-altitude frequency shift analysis
low <- subset(anolis,alt<=250)
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
  n.all.pre <- nrow(subset(df,ageClass==2))
  n.all.post <- nrow(subset(df,ageClass==4))
  chimatrix <- matrix(c(n.sp.pre,n.sp.post,n.all.pre,n.all.post),byrow = T, nrow = 2)
  chi <- chisq.test(chimatrix)
  c(i,round(n.sp.pre/n.all.pre,3),round(n.sp.post/n.all.post,3),round(chi$p.value,5))
}
chi2table <- data.frame(foreach(i=levels(factor(low$species)),.combine=rbind) %do% anolischi2fun(low,i),row.names=NULL)
names(chi2table) <- c("Species","relAbund.pre","relAbund.post" ,"p")
chi2table$p.adjusted <- p.adjust(as.numeric(as.character(chi2table$p)),'holm',6)
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

# shannon_diversity_corr <- function(x) {
#   n_singleton_sp <- ddply(x,.(species),summarize,n=length(locality)) %>% subset(n==1) %>% nrow()
#   C <- 1-n_singleton_sp/length(unique(x$species))
#   for(i in unique(x$species)){
#     sp <- subset(x,species==i)
#     p <- nrow(sp)/n
#   }
# }

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

richness_low <- subset(div_data,loc_group_alt<250) %>% ddply(c("locality_group","ageClass2"),function(e) length(unique(e$species)))
richness_low$altbin <- "0-250 m"
richness_mid <- subset(div_data,loc_group_alt>=250 & loc_group_alt<500) %>% ddply(c("locality_group","ageClass2"),function(e) length(unique(e$species)))
richness_mid$altbin <- "251-500 m"
richness_high <- subset(div_data,loc_group_alt>=500 & loc_group_alt<750) %>% ddply(c("locality_group","ageClass2"),function(e) length(unique(e$species)))
richness_high$altbin <- "501-750 m"
richness_vhigh <- subset(div_data,loc_group_alt>=750) %>% ddply(c("locality_group","ageClass2"),function(e) length(unique(e$species)))
richness_vhigh$altbin <- "750-1150 m"
richness <- rbind(richness_low,richness_mid,richness_high,richness_vhigh)
names(richness) <- c("locality_group","ageClass2","div","altbin")

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
  annotate(geom="text",label=c(""," "," "," "),x=1.5,y=1.65,size=5)

#regression of locality species richness by altitude, split by age class
loc_div <- ddply(comb,.(locality_group,ageClass,ageClass2),summarize,
                 alt=mean(na.omit(alt)),
                 div=shannon_diversity(data.frame(species)),
                 richness=length(unique(species)))
loc_div <- subset(loc_div,div>0)
loc_div$ageClass <- factor(loc_div$ageClass)
   
model1 <- lm(div~alt*ageClass,data = loc_div)
model2 <- lm(div~alt,data=loc_div)
summary(model1)
anova(model2,model1)

ggplot(loc_div,aes(x=alt,y=div,shape=ageClass2))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right",
        legend.background = element_blank(),
        text=element_text(size=10),
        axis.title=element_text(size=8))+
  #facet_grid(~ageClass2)+
  geom_point()+
  scale_shape_manual(values=c(1,20),name="Time Period")+
  scale_x_continuous(breaks=c(300,900))+
  scale_y_continuous(breaks=c(0,1))+
  xlab("Elevation (m)")+ylab("")+
  geom_smooth(method="glm",fill=NA,show.legend=F,aes(shape=NULL),formula = y~x,col="grey")

  
pdf("~/Dropbox/anolis/figures/final/figS3.pdf",width=3.5,height=4.5) 
grid.arrange(bin_plot,lm_plot,ncol=1,left="Shannon Diversity")
dev.off()

############# end elevation shift stats and figures #############
#################################################################

############# Bayesian model for drop in elevation across time periods ###########
#install.packages(c('coda','mvtnorm','devtools'))
#library(devltoos)
#devtools::install_github("rmcelreath/rethinking")

#per individual
library(rethinking)
comb$isrecent <- comb$timebin=="1980-2015"
locs$isrecent <- locs$ageClass==4
comb$alt[comb$alt<=0] <- 1
fits <- list();CIs <- list()
for(i in unique(comb$species)){
  sp <- subset(comb,species==i)[,c("alt","isrecent")]
  all <- subset(comb,species!=i)[,c("alt","isrecent")]
  sp <- na.omit(sp)
  all <- na.omit(all)
  mu1 <- mean(sp$alt[sp$isrecent==F])
  sd1 <- sd(sp$alt[sp$isrecent==F])
  mu2 <- mean(all$alt[all$isrecent==T]) - mean(all$alt[all$isrecent==F])
  sde <- sapply(1:100,function(e){
    allboot <- all[sample(1:nrow(all),nrow(all),replace=T),]
    mu2boot <- mean(allboot$alt[allboot$isrecent==T]) - mean(allboot$alt[allboot$isrecent==F])
    mu2boot
  }) %>% sd()
  fit <- map(alist(
    alt~dlnorm(log(mu),sigma), #altitude is normally distributed with mean mu and sd sigma
    mu <- a+effort*isrecent+diff*isrecent, #mean altitude differs between time periods by an amount diff
    a ~ dlnorm(log(mu1),sd1), #mean altitude in time period 1 is normally distributed with empirical mean and sd for each species
    effort ~ dnorm(mu2,sde), #difference due to survey effort is normall distributed with empirical mean and sd
    diff ~ dnorm(0,50),  #difference due to range shifts is normally distributed with mean 0 and sd 50.
    sigma ~ dunif(0,3)), #sd of modern altitude is uniform with bounds 0,3
    data=sp,
    start=list(a=300,effort=0,diff=0,sigma=1),
    control=list(maxit=1000))
  print(fit)
  fits <- append(fits,precis(fit,prob=0.95))
  CIs <- append(CIs,PI(extract.samples(fit)$diff,prob=0.95))
  print(i)
}
names(fits) <- unique(comb$species)
fits



#plots of final elevation curves using estimated model params
pdf("~/Dropbox/anolis/figures/model_curves.pdf",useDingbats = F,width=5.5,height=8)
par(mfrow=c(6,1),mar=c(2.5,2.5,1,2.5),cex.axis=0.7,tck=-0.03,mgp=c(1.2,.3,0),cex.main=1.1,cex.sub=0.8,bty="n")

hist(comb$alt[comb$species=="evermanni" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.01),
     xlab="",ylab="Density",main=expression(italic(A.~evermanni)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="evermanni" & comb$ageClass==2],30,xlim=c(0,1300),main=NA,ylim=c(0,0.01),
     xlab="",ylab="Density",axes=F,prob=T,col=alpha("black",0.5),border="white")
curve(dlnorm(x,log(398.34),0.91),0,1300,ylab="Density",xlab="",add=T,
      lty=2,ylim=c(0,0.005),lwd=1.5) #evermanni
curve(dlnorm(x,log(398.34-17.04-66.32),0.91),0,1300,ylab="Density",xlab="",
      add=T,lwd=1.5) #evermanni
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.5),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

hist(comb$alt[comb$species=="gundlachi" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.01),
     xlab="",ylab="Density",main=expression(italic(A.~gundlachi)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="gundlachi" & comb$ageClass==2],30,xlim=c(0,1300),ylim=c(0,0.01),
     xlab="",ylab="Density",axes=F,main=NA,
     prob=T,col=alpha("black",0.5),border="white")
curve(dlnorm(x,log(498.78),0.43),0,1300,ylab="Density",xlab="",
      main=expression(italic(gundlachi)),lty=2,ylim=c(0,0.005),lwd=1.5,add=T) #gundlachi
curve(dlnorm(x,log(498.78+22.76-171.08),0.43),0,1300,ylab="Density",xlab="",
      main=expression(italic(gundlachi)),add=T,lwd=1.5) #gundlachi
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.5),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

hist(comb$alt[comb$species=="krugi" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.006),
     xlab="",ylab="Density",main=expression(italic(A.~krugi)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="krugi" & comb$ageClass==2],30,xlim=c(0,1300),main=NA,ylim=c(0,0.006),
     xlab="",ylab="Density",axes=F,prob=T,col=alpha("black",0.5),border="white")
curve(dlnorm(x,log(354.37),0.74),0,1300,ylab="Density",xlab="",
      main=expression(italic(krugi)),lty=2,ylim=c(0,0.005),lwd=1.5,add=T) #krugi
curve(dlnorm(x,log(354.37+3.19-146.33),0.74),0,1300,ylab="Density",xlab="",
      main=expression(italic(krugi)),add=T,lwd=1.5) #krugi
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.5),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

hist(comb$alt[comb$species=="cristatellus" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.015),
     xlab="",ylab="Density",main=expression(italic(A.~cristatellus)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="cristatellus" & comb$ageClass==2],30,xlim=c(0,1300),main=NA,ylim=c(0,0.015),
     xlab="",ylab="Density",axes=F,prob=T,col=alpha("black",0.5),border="white")
curve(dlnorm(x,log(29.84),1.77),0,1300,ylab="Density",xlab="",add=T,
      main=expression(italic(cristatellus)),lty=2,ylim=c(0,0.015),lwd=1.5) #cristatellus
curve(dlnorm(x,log(29.84-33.51+73.44),1.77),0,1300,ylab="Density",xlab="",
      main=expression(italic(cristatellus)),add=T,lwd=1.5) #cristatellus
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.5),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

hist(comb$alt[comb$species=="pulchellus" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.03),
     xlab="",ylab="Density",main=expression(italic(A.~pulchellus)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="pulchellus" & comb$ageClass==2],30,xlim=c(0,1300),main=NA,ylim=c(0,0.03),
     xlab="",ylab="Density",axes=F,prob=T,col=alpha("black",0.5),border="white")
curve(dlnorm(x,log(32.47),1.89),0,1300,ylab="Density",xlab="",add=T,
      main=expression(italic(pulchellus)),lty=2,ylim=c(0,0.015),lwd=1.5) #pulchellus
curve(dlnorm(x,log(32.47-20.53+26.11),1.89),0,1300,ylab="Density",xlab="",
      main=expression(italic(pulchellus)),add=T,lwd=1.5) #pulchellus
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.5),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

hist(comb$alt[comb$species=="stratulus" & comb$ageClass==4],30,xlim=c(0,1300),ylim=c(0,0.015),
     xlab="Elevation (m)",ylab="Density",main=expression(italic(A.~stratulus)),
     prob=T,border="black")
par(new=T)
hist(comb$alt[comb$species=="stratulus" & comb$ageClass==2],30,xlim=c(0,1300),main=NA,ylim=c(0,0.015),
     xlab="Elevation (m)",ylab="Density",axes=F,prob=T,col=alpha("black",0.4),border="white")
curve(dlnorm(x,log(57.92),1.55),0,1300,ylab="Density",xlab="Elevation (m)",add=T,
      main=expression(italic(stratulus)),lty=2,ylim=c(0,0.015),lwd=1.5) #stratulus
curve(dlnorm(x,log(57.92-19.66+17.55),1.55),0,1300,ylab="Density",xlab="Elevation (m)",
      main=expression(italic(stratulus)),add=T,lwd=1.5) #stratulus
legend(x="topright",legend=c("1952-1977","1991-2015"),bty="n",
       pch=c(15,0),col=c(alpha("black",0.4),"black"))
legend(x="topright",legend=c("                   "," "),bty="n",
       lty=c(2,1),col="black")

dev.off()
par(mfrow=c(1,1))


#locality-based analysis
locs$isrecent <- locs$ageClass==4
locs$alt[locs$alt<=0] <- 1
fits <- list();CIs <- list()
locs <- subset(locs,species!="All Species")
for(i in unique(locs$species)){
  sp <- subset(locs,species==i)[,c("alt","isrecent")]
  all <- subset(locs,species!=i)[,c("alt","isrecent")]
  sp <- na.omit(sp)
  all <- na.omit(all)
  mu1 <- mean(sp$alt[sp$isrecent==F])
  sd1 <- sd(sp$alt[sp$isrecent==F])
  mu2 <- mean(all$alt[all$isrecent==T]) - mean(all$alt[all$isrecent==F])
  sde <- sapply(1:100,function(e){
    allboot <- all[sample(1:nrow(all),nrow(all),replace=T),]
    mu2boot <- mean(allboot$alt[allboot$isrecent==T]) - mean(allboot$alt[allboot$isrecent==F])
    mu2boot
  }) %>% sd()
  fit <- map(alist(
    alt~dlnorm(log(mu),sigma), #altitude is normally distributed with mean mu and sd sigma
    mu <- a+effort*isrecent+diff*isrecent, #mean altitude differs between time periods by an amount diff
    a ~ dlnorm(log(mu1),sd1), #mean altitude in time period 1 is normally distributed with empirical mean and sd for each species
    effort ~ dnorm(mu2,sde), #difference due to survey effort is normall distributed with empirical mean and sd
    diff ~ dnorm(0,50),  #difference due to range shifts is normally distributed with mean 0 and sd 50.
    sigma ~ dunif(0,3)), #sd of modern altitude is uniform with bounds 0,3
    data=sp,
    start=list(a=100,effort=0,diff=0,sigma=.3),
    control=list(maxit=1000))
  print(fit)
  fits <- append(fits,precis(fit,prob=0.95))
  CIs <- append(CIs,PI(extract.samples(fit)$diff,prob=0.95))
  print(i)
}
names(fits) <- unique(locs$species)
fits



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


