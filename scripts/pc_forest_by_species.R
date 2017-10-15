#find the % occurrence reports in forested areas by species & age class. (~15min run time)
library(rgeos)
pc_forest_table <- ddply(comb,.(species),function(e) {

  sp_age2 <- subset(e,ageClass==2 & is.na(long)==F)
  sp_age2_pts <- SpatialPointsDataFrame(data.frame(sp_age2$long,sp_age2$lat),data=data.frame(sp_age2$uncertainty),
                               proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
                    spTransform(crs(age2)) %>% 
                      gBuffer(width=sp_age2$uncertainty,byid=T)
  age2_pc_forest <- extract(age2,sp_age2_pts) %>% 
                      lapply(function(e) mean(e,na.rm=T)) %>%
                        unlist() %>% 
                          mean(na.rm=T)

  sp_age4 <- subset(e,ageClass==4 & is.na(long)==F)
  sp_age4_pts <- SpatialPointsDataFrame(data.frame(sp_age4$long,sp_age4$lat),data=data.frame(sp_age4$uncertainty),
                                        proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
                  spTransform(crs(age4)) %>% 
                    gBuffer(width=sp_age4$uncertainty,byid=T)
  age4_pc_forest <- extract(age4,sp_age4_pts) %>% 
                      lapply(function(e) mean(e,na.rm=T)) %>%
                        unlist() %>% 
                          mean(na.rm=T)
  
  c(age2_pc_forest,age4_pc_forest)
})
  