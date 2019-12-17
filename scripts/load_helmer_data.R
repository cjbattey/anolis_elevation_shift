library(raster)
ext.pr <- extent(-67.6,-65.2,17.8,18.7)
ext.yunque <- extent(825000,850000,2015000,2040000)
proj4.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4.utm <- "+proj=utm +zone=19 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#1. load landcover raster (Helmer et al. 2008, JGRBiogeoscience), reclassify to combine soil/forest types as needed
alt <- raster("data/elevation/alt_1s_UTM19N.tif")
lc <- raster("data/landcover/prforestage_foresttype/IITF_JGR113_puertorico_fortype_age.img") %>% crop(extent(alt))
m.all <- read.delim("data/landcover/prforestage_foresttype/combineSoils_reclassMatrix.txt",header=F)
m.forests <- read.delim("data/landcover/prforestage_foresttype/PRlandcover_combSoils_allForest.txt",header=F)
m.humidForests <- read.delim("data/landcover/prforestage_foresttype/PRlandcover_combSoils_humidForest.txt",header=F)
m.dryForests <- read.delim("data/landcover/prforestage_foresttype/PRlandcover_combSoils_dryForest.txt",header=F)
lc.all <- reclassify(lc,m.all)
lc.forests <- reclassify(lc,m.forests)
lc.humidForests <- reclassify(lc,m.humidForests)
lc.dryForests <- reclassify(lc,m.dryForests)
coasts <- shapefile("data/ne_10m_coastline/ne_10m_coastline.shp") %>% crop(ext.pr) %>% spTransform(.,crs(lc))
#writeRaster(lc,"./data/landcover/prforestage_foresttype/PRlandcover_combSoils.tif")
#lc.all classfication key:
# 0	Bare
# 1	High Density Urban 1991
# 2	New Urban and Bare 1991-2000
# 3	Low Density Development 2000
# 4	Pasture and Agriculture
# 5	Forest Age 1-9yr
# 6	Forest Age 10-22yr
# 7	Forest Age 23-49yr
# 8	Forest Age 50-64yr
# 9 Emergent Wetland
# 10  Coastal Sand/Rock
# 11  Water

#2. reclassify by age to produce rasters by age class
age1 <- reclassify(lc.all,as.matrix(data.frame(from=c(0:11),to=c(0,0,0,0,0,0,0,0,1,0,0,0))))
age2 <- reclassify(lc.all,as.matrix(data.frame(from=c(0:11),to=c(0,0,0,0,0,0,0,1,1,0,0,0)))) 
age3 <- reclassify(lc.all,as.matrix(data.frame(from=c(0:11),to=c(0,0,0,0,0,0,1,1,1,0,0,0))))
age4 <- reclassify(lc.all,as.matrix(data.frame(from=c(0:11),to=c(0,0,0,0,0,1,1,1,1,0,0,0))))

# par(mfrow=c(2,2),oma = c(1,1,0,0) + 0.1,mar = c(2,2,2,2) + 0.1)
# plot(age1,main="1935-1951",axes=F,legend=F)
# plot(coasts,lwd=.5,add=T)
# plot(age2==1,main="1952-1977",axes=F,legend=F)
# plot(coasts,lwd=.5,add=T)
# plot(age3==1,main="1978-1990",axes=F,legend=F)
# plot(coasts,lwd=.5,add=T)
# plot(age4==1,main="1991-2000",axes=F,legend=F)
# plot(coasts,lwd=.5,add=T)
# par(mfrow=c(1,1))