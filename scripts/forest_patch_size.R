#experimental forest fragment size-change analyses
library(spatialEco);library(raster);library(SDMTools)
source("scripts/load_helmer_data.R")
age1 <- reclassify(age1,t(matrix(c(0,NA))))
age2 <- reclassify(age2,t(matrix(c(0,NA))))
age3 <- reclassify(age3,t(matrix(c(0,NA))))
age4 <- reclassify(age4,t(matrix(c(0,NA))))
alt2 <- projectRaster(alt,age1)

#mask for subsutteing by elevations
low <- reclassify(alt2,t(matrix(c(251,2000,NA))))
mid <- reclassify(alt2,matrix(c(-200,250,NA,
                                 501,2000,NA),nrow=2,byrow=T))
midhigh <- reclassify(alt2,matrix(c(-200,500,NA,
                                     751,2000,NA),nrow=2,byrow=T))
high <- reclassify(alt2,t(matrix(c(-200,750,NA))))


#empty raster for resampling to 500m resolution
age2low <- mask(age2,low)
age2low <- aggregate(age2,fact=20, fun=modal, expand=TRUE, na.rm=TRUE)
age2lowlab <- ConnCompLabel(age2)
age2patch <- PatchStat(age2lowlab)
age2patch <- subset(age2patch,area>1)
age2patch$ageclass <- 2

age4low <- mask(age4,low)
age4low <- aggregate(age4,fact=20, fun=modal, expand=TRUE, na.rm=TRUE)
age4lowlab <- ConnCompLabel(age4)
age4patch <- PatchStat(age4lowlab)
age4patch <- subset(age4patch,area>1)
age4patch$ageclass <- 4

tmp <- rbind(age2patch,age4patch)
ggplot(data=tmp,aes(x=log(area)))+theme_minimal()+facet_grid(~ageclass)+geom_histogram()
median(age2patch$area)
median(age4patch$area)

