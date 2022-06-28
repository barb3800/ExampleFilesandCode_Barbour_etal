## Code to Create "Pretty" Final Maps & Quantitative Indices

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

##############################################
# SET UP 
##############################################

# load packages (make sure you have installed them first!)
library(sf);library(dplyr);library(raster);library(ggplot2);library(ggspatial);library(tidyr);library(reshape2)

# set working directory
setwd("yourfilepathhere")

## load in data
load(file="./Example_Data/combined_monthly_overlap.rda")
load(file="./Example_Data/relative_overlap_matrices.rda")

##############################################
# PRETTY MAP
##############################################

# before mapping, check the distribution of values in your data:
## do this on the matrix version of your data
df<-as.data.frame(rbind(g1_s3_d[[1]],g1_s3_d[[2]],g1_s3_d[[3]],g1_s3_d[[4]],g1_s3_d[[5]],g1_s3_d[[6]],g1_s3_d[[7]],g1_s3_d[[8]],g1_s3_d[[9]],g1_s3_d[[10]],g1_s3_d[[11]],g1_s3_d[[12]]))

df2<-gather(df, drop, values, V1:V42)

boxplot(df2$values)

range(df2$values)

# NOTE: there is only one outlier value = 0.02- the rest are 0.01 or less
df2[which(df2$values >= 0.02),]

hist(df2$values,breaks=100)
## you can use these values to play with the breaks for your color bars in your map below!
## notice that much of this data is below 0.005 (lots of low values, due to not much overlap in study range)

# I really like using the "ggspatial" package for mapping
# it has open-source base maps (i use the "osm" one below)
# you can also add things like a scale bar, north arrow, and multiple layers
# for plotting rasters, I use the "layer_spatial" function
# for plotting sf objects (polygons, points), I use the "geom_sf" function

# NOTE: set the zoom option based on the extent of your data- if it's a large area, use a small zoom (e.g., 4) and if it's a small area, use a large zoom (e.g., 11)
## when in doubt, start off with a low zoom value and make bigger as needed!

# store maps for each month in a list:
g1_s3_plot<-c()
for (i in c(1:length(g1_s3_r))){
  g1_s3_plot[[i]]<-ggplot()+
    annotation_map_tile(type = 'osm', zoom = 4) + # set zoom here
    annotation_scale()+
    annotation_north_arrow(height=unit(0.5,"cm"),width=unit(0.5,"cm"),pad_y = unit(1,"cm"))+
    layer_spatial(data = g1_s3_r[[i]], aes(fill = stat(band1))) + # "stat(band1)" is used for rasters
    #shadow_spatial(data=box)+
    ylab("Latitude")+xlab("Longitude")+
    scale_fill_gradientn("Overlap Index",colours=c("lightblue","yellow","orange","red"),
                         breaks=c(0.0005,0.005,0.01,0.02), # here you can set break values for your diff colors based on the values in your data and your goal of highlighting "high overlap" areas
                         limits=c(0.0005,0.02), # you can also set the limits- you can exclude very low an doutlier high values for example to improve visualization of hotspots
                         na.value = NA)+
    theme(legend.key.height= unit(0.7, 'cm'),legend.key.width= unit(0.7, 'cm'), 
          legend.text = element_text(size=10),legend.title = element_text(size=10))
}

# print out maps (month 3 as example):
g1_s3_plot[[3]]

#######################################################
# Quantitative Indices
#######################################################

# Step 1: extract "high overlap" values and sum for each  matrix

# you can use the 75th-100th quantiles to determine "high" overlap values
## however, because our data is very skewed, quantiles may be less informative
## this is likely because we are only using 3 individuals in our example data and so there is very little spatial overlap (lots of low values for overlap)

quantile(df2$values)
## 75th-100th quantiles: > 0.0005

# to really capture "high values", you can use the maps - "high values" (yellow/red) or hotspots on the map seem to be around values > 0.005
o_g1_s3 <- lapply(g1_s3_d, function(x) {
  sum(x[which(x>0.005)])
})

# Step 2: "melt" data to turn into a dataframe
d_g1_s3 <- lapply(o_g1_s3, function(x) {
  melt(x)
})

# Step 3: # bind months together, add columns for month and geartype
geartype<-c("Drifting Longline") # note: if you had more geartypes, you could add them here!

g1_s3<-do.call("rbind",d_g1_s3)
g1_s3$Month<-c(1:nrow(g1_s3))
g1_s3$State<-rep(3,nrow(g1_s3))
g1_s3$Gear<-rep(geartype[1],nrow(g1_s3))

# you can modify this code to add in additional code for other behavior states and geartypes

# create a nicer looking data table:
names(g1_s3)[1]<-"Overlap_Index"

Month_Gear_State<-g1_s3 %>% group_by(Month,Gear,State) %>%
  summarize(Index=sum(Overlap_Index)) %>%
  data.frame()

Month_Gear_State[,c(4)]<-signif(Month_Gear_State[,c(4)],digits=2)
MGS_table<- Month_Gear_State %>% spread(Gear,Index) 

# print table
MGS_table
## the patterns seen in the table should reflect the general "hotspot" patterns we see in the maps
## here, there is especially high overlap in months: 8 and 9
## however, in general there is little overlap due to needing to incorporate more tracking data (we only used 3 ids as example)

