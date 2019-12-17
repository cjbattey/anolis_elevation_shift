#Puerto Rico historical climate analyses
library(data.table);library(magrittr);library(plyr);library(ggplot2)

#load data
clim <- suppressWarnings(fread("data/climate/NCDC_longTermPRstations_742067.csv",na.strings="-9999"))
#clim <- subset(clim,STATION_NAME %in% c("AGUIRRE US","LAJAS SUBSTATION US","MAGUEYES ISLAND US","MANATI 2 E US","RIO PIEDRAS EXPERIMENTAL STATION US","ROOSEVELT ROADS US"))
clim$year <- as.numeric(substr(clim$DATE,1,4))
clim$month <- as.numeric(substr(clim$DATE,5,6))
clim <- subset(clim,year < 2016)
clim$LATITUDE <- as.numeric(clim$LATITUDE)
clim$LONGITUDE <- as.numeric(clim$LONGITUDE)
clim <- clim[is.na(clim$MNTM)==F]
  
#remove station-years with < 12 months
data <- clim[1,]
for(i in levels(factor(clim$STATION_NAME))){
  station <- subset(clim,STATION_NAME==i)
  for (j in min(station$year):max(station$year)){
    year <- subset(station,year==j)
    n.months <- nrow(year)
    if(n.months==12){
      data <- rbind(data,year)
    }
  }
}

#get altitude of stations
alt <- crop(raster("data/elevation/alt_30s_bil/alt.bil"),ext.pr)
coords <- data.frame(subset(data,!is.na(LATITUDE)))[,c("LONGITUDE","LATITUDE")]
data$alt[!is.na(data$LATITUDE)]<- raster::extract(alt,SpatialPoints(coords))

#convert to metric
data$MNTM_C <- (data$MNTM-32)*(5/9)
data$TPCP_cm <- data$TPCP*2.54
  
#get temperature/altitude lapse rate
avg <- ddply(data,.(STATION_NAME,alt),summarize,temp=mean(MNTM_C),precip=mean(TPCP_cm))
lm(temp~alt,data=avg) %>% summary()
plot(avg$alt,avg$temp,xlab="Elevation",ylab="Average Annual Temperature")+abline(lm(temp~alt,data=avg),cex.lab=0.5)

#drop short-term stations
station_years <- ddply(data,.(STATION_NAME),summarize,minyear=min(year),nyears=length(month)/12) %>% 
  subset(nyears>40)
data <- subset(data,STATION_NAME %in% station_years$STATION_NAME)

#get average annual temperatures
avg <- ddply(data,.(year,STATION_NAME,alt),summarize,temp=mean(MNTM_C),precip=sum(TPCP_cm))
mean_alt <- ddply(avg,.(STATION_NAME),summarize,mean_alt=mean(na.omit(alt)))
data <- merge(data,mean_alt,by="STATION_NAME")
avg <- merge(avg,mean_alt,by="STATION_NAME")
avg$shortname <- factor(avg$STATION_NAME)
levels(avg$shortname) <- c("Borinquen","Corozal","Dos Bocas","Lajas","Manati","Ponce","Rio Piedras","Roosevelt")

#t-test for difference of means 1952-1977 vs 1991-2015
precip <- ddply(avg,.(STATION_NAME),function(e){
  t <- t.test(subset(e,year>1952 & year <=1977)$precip,subset(e,year>1991 & year<2016)$precip)
  a <- c(round(t$estimate[1],2),round(t$estimate[2],2),round(t$estimate[2]-t$estimate[1],2),round(t$p.value,3),round(t$parameter,2))
  names(a) <- c("1952-1977","1991-2015","Difference","p","df")
  a
})
mean(precip$Difference)
abs(mean(precip$Difference))/mean(precip$`1952-1977`)

ddply(avg,.(STATION_NAME),function(e){
  t <- t.test(subset(e,year>1952 & year <=1977)$temp,subset(e,year>1991 & year<2016)$temp)
  a <- c(round(t$estimate[1],5),round(t$estimate[2],5),round(t$estimate[2]-t$estimate[1],5),round(t$p.value,5),round(t$parameter,2))
  names(a) <- c("1952-1977","1991-2015","Difference","p","df")
  a
}) %>% xtable()
  
  #.$Difference %>% mean()

#summary of regressions by station
regressions <- ddply(avg,.(STATION_NAME,mean_alt),
            function(e) { 
              model <- lm(e$temp ~ e$year)
              p <- summary(model)$coefficients[8] %>% round(3)
              slope <- coef(model)[2] %>% round(3)
              r2 <- summary(model)$r.squared %>% round(3)
              c(p=p,r2=r2,slope=slope)
            })
names(regressions) <- c("Station","Elevation (M)","p","R^2","slope (C/yr)")
sig.regressions <- subset(regressions,p<0.05)
lm(avg$temp~avg$year) %>% summary() #average across all stations

#line graph of temperature over time, plus overall trend
temp <- ggplot(data=avg,aes(x=year,y=temp))+theme_minimal()+facet_wrap(~shortname,nrow=3)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=9),
        text=element_text(size=8))+
  xlab("")+
  ylab("Mean Annual Temperature (C)")+
  scale_x_continuous(breaks=c(1960,2000))+
  scale_y_continuous(breaks=c(24,26,28))+
  geom_line(lwd=0.25,col="grey")+
  annotate(geom="Text",x=1953,y=27.5,label=c("A","B","C","D","E","F","G","H"),size=3.5)+
  annotate(geom="Text",x=1959,y=27.5,label=c("*","*","*","*","","","",""),size=3.5)+
  geom_line(data=subset(avg,year>1951&year<1977))+
  geom_line(data=subset(avg,year>1990&year<2016))

precip <- ggplot(data=avg,aes(x=year,y=precip))+theme_minimal()+facet_wrap(~shortname,nrow=3)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=9),
        text=element_text(size=8))+
  xlab("")+
  ylab("Mean Annual Precipitation (cm)")+
  scale_x_continuous(breaks=c(1960,2000))+
  #scale_y_continuous(breaks=c(24,26,28))+
  geom_line(lwd=0.25,col="grey")+
  annotate(geom="Text",x=1953,y=275,label=c("A","B","C","D","E","F","G","H"),size=3.5)+
  annotate(geom="Text",x=1959,y=275,label=c("","*","","","","","*",""),size=3.5)+
  geom_line(data=subset(avg,year>1951&year<1977))+
  geom_line(data=subset(avg,year>1990&year<2016))


#get landcover class of weather stations
stationLocations <- ddply(data,.(STATION_NAME),summarize,long=mean(na.omit(LONGITUDE)),lat=mean(na.omit(LATITUDE)))
stations <- SpatialPointsDataFrame(data.frame(stationLocations$long,stationLocations$lat),
                                   data=data.frame(stationLocations$STATION_NAME),
                                   proj4string = crs(coasts))

urban_stations <- spTransform(stations,crs(lc))
urban_stations <- urban_stations$stationLocations.STATION_NAME[raster::extract(lc,urban_stations)==1]

#map of stations
#r <- raster(xmn=-67.6,xmx=-65.2,ymn=17.8,ymx=18.7,res=.02)
#alt_tile <- resample(alt,r) %>% as.data.frame(xy=T)
stations <- data.frame(stations@coords,stations@data)
names(stations) <- c("long","lat","STATION_NAME")

map <- ggplot()+coord_map()+theme_minimal()+xlab("")+ylab("")+ #ggplot sucks at rasters :(
  theme(axis.ticks = element_blank(),axis.text = element_blank(),legend.position = "right")+
  geom_path(data=fortify(coasts),aes(x=long,y=lat,group=group),lwd=0.25)+
  #geom_tile(data=alt_tile,aes(x=x,y=y,fill=alt),col=NA)+
  geom_point(data=stations,aes(x=long,y=lat))+
  geom_text(data=stations,aes(x=long,y=lat,label=c("A","B","C","D","E","F","G","H")),
            size=2.75,position = position_nudge(.06,-.04))

pdf("~/Dropbox/anolis/figures/final/FigS4.pdf",width=6,height=4.5)
gridExtra::grid.arrange(temp,precip,map,layout_matrix=matrix(c(1,2,
                                                               1,2,
                                                               1,2,
                                                               3,NA),nrow=4,ncol=2,byrow=T))
dev.off()
