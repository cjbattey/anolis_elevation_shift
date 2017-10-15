#load, filter, and georeference occurrence records
library(data.table);library(raster);library(ggplot2);library(plyr);library(foreach);library(stringr);
library(magrittr);library(dismo);library(rgeos)

################################# NOTE: DO NOT OPEN TABLES WITH EXCEL ##############################

#############################################################
################# Read in GBIF + Museum Data ################
#############################################################
#projections and extents for maps
ext.pr <- extent(-67.6,-65.2,17.8,18.7)
ext.yunque <- extent(825000,850000,2015000,2040000)
proj4.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4.utm <- "+proj=utm +zone=19 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load elevation rasters
# alt <- raster("~/Dropbox/anolis/data/alt_30s_bil/alt.bil") %>% crop(ext.pr) %>% 
#   projectRaster(.,projectExtent(.,proj4.utm))
# alt1s<- raster("~/Dropbox/anolis/data/elevation/PuertoRico_altitude_1sec_USGSNED2013.tif") %>% crop(ext.pr) %>% 
#   projectRaster(.,projectExtent(.,proj4.utm))

alt <- raster("data/alt_30s_UTM19N.tif")
alt1s <- raster("data/alt_1s_UTM19N.tif")

#remove embedded nulls & read in data
file <- "data/gbif_PRanolis/occurrence.txt"
tt <- tempfile()  
system(paste0("tr < ", file, " -d '\\000' >", tt))
anolis <- data.frame(fread(tt))

#pull useful columns from GBIF
anolis <- anolis[names(anolis) %in% c("gbifID","institutionCode","basisOfRecord","catalogNumber","recordedBy",
                                      "eventDate","verbatimEventDate","year","month","day","countryCode",
                                      "stateProvince","county","municipality","locality",
                                      "verbatimLocality","verbatimElevation","locationRemarks","decimalLatitude",
                                      "decimalLongitude","coordinateUncertaintyInMeters","georeferenceProtocol",
                                      "georeferenceSources","georeferenceRemarks","species","infraspecificEpithet","elevation")]
anolis <- subset(anolis,basisOfRecord=="PRESERVED_SPECIMEN" & species != "") 
anolis <- subset(anolis,locality != "")

#merge FMNH full localities
fmnh.mus <- read.csv("~/Dropbox/anolis/data/locs/fmnh_anolis.csv")
names(fmnh.mus) <- c("catalogNumber","collectionNumber","Taxon","verbatimLocality","Details","verbatimEventDate","Details2",
                     "recordedBy","IDnotes","specimen.loc","habitatNotes","sex")
fmnh.gbif <- subset(anolis,institutionCode == "FMNH") %>% arrange(.,catalogNumber) #sort by catalogNumber bc merge will sort on the "by" column
anolis <- subset(anolis,institutionCode != "FMNH")
fmnh.merge <- merge(fmnh.gbif,fmnh.mus,by="catalogNumber")
fmnh.gbif$locality <- fmnh.merge$verbatimLocality.y
fmnh.gbif$verbatimEventDate <- fmnh.merge$verbatimEventDate.y
anolis <- rbind(anolis,fmnh.gbif)

#Remove other missing/withheld/duplicate institutions (note ROM is a duplicate institution ID in GBIF)
anolis <- subset(anolis,institutionCode %in% c("AMNH","Royal Ontario Museum: ROM")==F)

#######################################################################
##################### Missing/Ambiguous Dates #########################
#######################################################################
baddates <- subset(anolis,is.na(year))                                          # split out reports missing a year
anolis <- subset(anolis,!is.na(year))
baddates$verbatimEventDate[grep("no",baddates$verbatimEventDate)] <- NA         # if no date was recorded, standardize to NA
baddates$verbatimEventDate[which(baddates$verbatimEventDate == "")] <- NA 

yrange <- as.character(c(1800:2016))                                          
for(i in 1:nrow(baddates)){                                                     # for each row
  e <- baddates[i,]
  y <- str_locate(e$verbatimEventDate,yrange)                                 
  if(!all(is.na(y))){                                                           # if a year-range 4-digit numeral is in the verbatim date, use it. 
    start <- na.omit(y[,1])
    stop <- na.omit(y[,2])
    baddates$year[i] <- substr(e$verbatimEventDate,start,stop) %>% as.numeric()
  } else if(length(unlist(strsplit(e$verbatimEventDate,"/|-"))) == 3) {         # if the date splits in 3 on / or -, take the last element and add 19/20 as appropriate
    yrend <- unlist(strsplit(e$verbatimEventDate,c("/|-")))[3]                  # helpfully, there are no reports in ambiguous formats ending in 0-16
    if(as.numeric(yrend) > 16){ 
      yr <- as.numeric(paste("19",yrend,sep=""))
      baddates$year[i] <- yr
    } else {
      yr <- as.numeric(paste("20",yrend,sep=""))
      baddates$year[i] <- yr
    }
  }
}
baddates <- subset(baddates,!is.na(year))
anolis <- rbind(anolis,baddates)                                                 # recombine reports

##################################################################
################ Adding Verbatim Elevations ######################
##################################################################
#for MCZ specimens with "VERBATIM ELEVATION:" in the verbatim locality, convert as needed and copy to elevation column
for(i in 1:nrow(anolis)){  
  row <- anolis[i,]
  if(grepl("ELEVATION",row$locality)){                                            # if ELEVATION in locality, copy text after ":"
    ve.string <- strsplit(row$locality,"ELEVATION:") %>% unlist() %>% .[2]
    if(length(unlist(strsplit(ve.string,"-")))==2){                               # if a range is provided, take mean.
      range1 <- strsplit(ve.string,"-|[A-z]|'") %>% unlist() %>% .[1] %>% as.numeric()
      range2 <- strsplit(ve.string,"-|[A-z]|'") %>% unlist() %>% .[2] %>% as.numeric()
      ve <- mean(c(range1,range2))
    } else if(grepl("ca.",ve.string)){                                            
      ve <- strsplit(ve.string,"ca.|[A-Za-z]|'|-") %>% unlist() %>% .[2] %>% as.numeric()
    } else {                                        
      ve <- as.numeric(unlist(strsplit(ve.string,"[A-Za-z]|'|-"))[1])
    }
    
    if(ve > 1400 | grepl("ft|FT|feet|Feet|'",ve.string)){                         # convert to meters as needed
      ve <- ve * 0.3048
      anolis$elevation[i] <- ve
    } 
    if(grepl("M|m|meters|Meters",ve.string)){                 
      anolis$elevation[i] <- ve
    }
  }
}

#############################################################
################ georeferencing + altitude ##################
#############################################################
#add event ID for merge/split 
anolis$eventID <- paste(anolis$institutionCode,anolis$locality,anolis$year,sep=",")

#ddply summary for manual georeferencing
#loc <- ddply(anolis,.(locality,year,eventID,recordedBy,georeferenceSources,decimalLatitude,decimalLongitude,verbatimLocality),summarize,n=length(gbifID))
#write.csv(loc,"anolis_localities_27sept2016.csv",row.names = F,fileEncoding = "UTF-8")

#load georeferenced coordinates & use GPS where available
localities <- read.csv("~/Dropbox/anolis/data/anolis_localities_30sept2016.csv",stringsAsFactors = F) %>% subset(!is.na(lat)) #note fread caused encoding issues with double-quotes.
localities$long[grepl("GPS|gps",localities$georeferenceSources)] <- localities$decimalLongitude[grepl("GPS|gps",localities$georeferenceSources)]
localities$lat[grepl("GPS|gps",localities$georeferenceSources)] <- localities$decimalLatitude[grepl("GPS|gps",localities$georeferenceSources)]
localities$uncertainty[grepl("GPS|gps",localities$georeferenceSources)] <- 30

#get average elevation across uncertainty radius (21 warnings for Mona Island pts w NA elevations)
localities$alt <- data.frame(localities$long,localities$lat) %>%
                    SpatialPoints(proj4string=crs(proj4.wgs)) %>%
                      spTransform(proj4.utm) %>%
                        gBuffer(byid=T,width=localities$uncertainty) %>% 
                          extract(alt1s,.) %>%
                            lapply(.,function(e) mean(na.omit(e))) %>%
                              unlist()

#point elevation for comparison
localities$pt.alt <- data.frame(localities$long,localities$lat) %>%
                      SpatialPoints(proj4string=crs(proj4.wgs)) %>%
                        spTransform(proj4.utm) %>%
                          extract(alt1s,.)

tmp<- merge(anolis,localities,all.x=T,all.y=T)
#localities$verbatimLocality <- gsub("\"\"\"\"","\"\"",localities$verbatimLocality) #double quotes were preventing join...
#t <- join(anolis,localities,by=c("eventID","verbatimLocality"),match="first") #don't join by other vars - encoding issues...

# use verbatim elevation where available
anolis$alt[!is.na(anolis$elevation)] <- anolis$elevation[!is.na(anolis$elevation)]
anolis$pt.alt[!is.na(anolis$elevation)] <- anolis$elevation[!is.na(anolis$elevation)]

#get distance between GBIF and UW georeferencing
anolis$georefDistance <- NA
for(i in 1:nrow(anolis)){
  row <- anolis[i,]
  if(!is.na(row$decimalLatitude) & !is.na(row$lat)){
    gbif <- SpatialPoints(data.frame(row$decimalLongitude,row$decimalLatitude),proj4string=crs(proj4.wgs))
    uw <- SpatialPoints(data.frame(row$long,row$lat),proj4string=crs(proj4.wgs))
    anolis$georefDistance[i] <- spDistsN1(gbif,uw,longlat=T)
  }
}

#add age classes
anolis$ageClass[anolis$year<=1935] <- 0
anolis$ageClass[anolis$year>1935 & anolis$year<=1951] <- 1
anolis$ageClass[anolis$year>1951 & anolis$year<=1977] <- 2
anolis$ageClass[anolis$year>1977 & anolis$year<=1990] <- 3
anolis$ageClass[anolis$year>1990] <- 4

#cluster localities (deprecated - use clustering from elevation_analyses.R)
# coords <- subset(anolis,!is.na(lat))                                         #pull rows with coordinates
# coords <- ddply(coords,.(eventID,long,lat,ageClass),summarize,nobs=length(day)) #summarize by locality+ageClass
# pts <- SpatialPoints(coords=data.frame(coords$long,coords$lat),proj4string=crs(proj4.wgs))
# dist <- spDists(pts,pts) %>% as.dist()                                       #get pairwise distance matrix
# fit <- hclust(dist,method="average")                                         #cluster
# clusters <- cutree(fit,h=1)                                                  #group points within h km
# coords$localityGroup <- clusters                                             #add clusters to locality summary
# obs.per.group <- ddply(coords,.(localityGroup,ageClass),summarize,obs.per.group=sum(nobs))
# coords <- join(coords,obs.per.group,by=c("localityGroup","ageClass"))
# anolis <- join(anolis,coords,by=c("eventID","long","lat","ageClass"))

#binary columns for easy binning (deprecated - use ageClass for binning w/Helmer dataset)
# anolis$timebin <- factor(anolis$year < 1980)
# levels(anolis$timebin) <- c("1980-2015","1955-1979")
# anolis$timebin <- factor(anolis$timebin,levels=c("1955-1979","1980-2015"))

#final data filters
anolis.full <- anolis
anolis <- subset(anolis,(uncertainty <= 2000) | (is.na(uncertainty) & !is.na(elevation)))
good.species <- ddply(anolis,.(species,timebin),summarize,n=length(day)) %>% subset(.,n>200) %>% .$species
anolis <- subset(anolis,species %in% good.species)
#write.table(anolis,"~/Dropbox/anolis/data/anolis_data.csv",sep=",")

#add sight records
sight <- read.csv("~/Dropbox/anolis/data/locs/Anolis_sight_records.csv")
sight$basisOfRecord <- "sight"
sight$alt <- data.frame(sight$long,sight$lat) %>%             #get uncertainty-buffered elevation
              SpatialPoints(proj4string=crs(proj4.wgs)) %>%
                spTransform(proj4.utm) %>%
                  gBuffer(byid=T,width=sight$uncertainty) %>% 
                    extract(alt1s,.) %>%
                      lapply(.,function(e) mean(na.omit(e))) %>%
                        unlist()
sight$pt.alt <- data.frame(sight$long,sight$lat) %>%        #get point elevation
                  SpatialPoints(proj4string=crs(proj4.wgs)) %>%
                    spTransform(proj4.utm) %>%
                        extract(alt1s,.) %>%
                          lapply(.,function(e) mean(na.omit(e))) %>%
                            unlist()
sight$alt[!is.na(sight$Elev)] <- sight$Elev[!is.na(sight$Elev)]   #use verbatim elevations where available
sight$pt.alt[!is.na(sight$Elev)] <- sight$Elev[!is.na(sight$Elev)] 
#write.csv(sight,"~/Dropbox/anolis/data/locs/Anolis_sight_records_georeferenced.csv")



comb <- rbind.fill(anolis,sight)
coords <- subset(comb,!is.na(lat))                                         #pull rows with coordinates
coords <- ddply(coords,.(long,lat,ageClass),summarize,nobs=length(day)) #summarize by locality+ageClass
pts <- SpatialPoints(coords=data.frame(coords$long,coords$lat),proj4string=crs(proj4.wgs))
dist <- spDists(pts,pts) %>% as.dist()                                       #get pairwise distance matrix
fit <- hclust(dist,method="average")                                         #cluster
clusters <- cutree(fit,h=1)                                                  #group points within h km
coords$localityGroup <- clusters                                             #add clusters to locality summary
obs.per.group <- ddply(coords,.(localityGroup,ageClass),summarize,obs.per.group=sum(nobs))
coords <- join(coords,obs.per.group,by=c("localityGroup","ageClass"))
comb <- join(comb,coords,by=c("long","lat","ageClass"))


