## Code to implement population-level, 3D Utilization Distributions for EP data
## Using the AdehabitatHR package:
## Estimate the utilization distribution (UD) in 3 dimensions: space (x,y) and time, using the
## product kernel estimator algorithm (Keating and Cherry 2009, Horne et at. 2007)
## It places a 3-D kernel function (X-Y-T) over each relocation and then sums the kernel functions
## Time can be either linear or circular (Jan 3 2001=Jan 3 2002)

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

##############################################
# SET UP 
##############################################

# load packages (make sure you have installed them first!)
library(adehabitatHR);library(raster);library(sp);library(dplyr);library(sf);library(lubridate);library(rgdal)

# set working directory
setwd("yourfilepathhere")

## load in data
# jittered tracking data (list of ids)
load(file="./Example_Data/Jitter_Turtle_Data.rda")
# bind together to makr one big dataframe
turtle_id_all<-do.call("rbind",turtle_jitter) 

#################################################
# PREP DATA
################################################

# make data into an sf (spatial) object
turtle_id_sf<-turtle_id_all %>% st_as_sf (coords = c("lon_j","lat_j"), crs = 4326)

# extract DOY (day) from datetime
turtle_id_all$doy<-lubridate::yday(turtle_id_all$date)

# split back by id
id<-split(turtle_id_all,turtle_id_all$Ref)

# NOTE: it is important that all your layers have the same spatial extent and resolution!
# make grid that is 1x1 deg res for study area (and same extent as fishing data)
base<-raster()
extent(base)<-extent(-167.99,-125.99,-5.01,40.99) # set extent to same as fishing data
crs(base)<-crs(turtle_id_sf)
res(base)<-c(1,1) # 1 deg resolution

#convert raster to a grid
# in deg
grid<-base %>% 
  as("SpatialPixels") 

########################################################
# find 3D UDs for each id, using the base raster created above

# make empty list 
ud_id<-c()
# iteratively, add sublist of uds in ud form
for (i in c(1:length(id))){
  ## create month variable
  id[[i]]$month<-as.numeric(format(id[[i]]$date,"%m"))
  # extract last date
  end<-nrow(id[[i]])
  #  data frame of lat, lon, and day of year column
  xyt<-data.frame(x=id[[i]]$lon_j,y=id[[i]]$lat_j,t=id[[i]]$doy)
  # vector of days observed for each id
  vv <-lubridate::yday(seq(as.Date(id[[i]]$date[1]), 
                           as.Date(id[[i]]$date[end]),by = "month", tz="UTC"))
  # sequence of month names
  nn <-format(seq(as.Date(id[[i]]$date[1]), 
                  as.Date(id[[i]]$date[end]),by = "month", tz="UTC"),"%m")
  # loop through time vector and find ud at each point
  ## using a smoothing parameter of 3 deg and 30 time units
  ud_id[[i]]<- lapply(1:length(vv),function(x){
    kernelkcbase(xyt, h=c(3,3,30), tcalc=vv[x], grid=grid)
  })
  # apply unique months as names to ud list for each
  names(ud_id[[i]])<-nn
}

# if you get it, don't worry about warning "In proj4string(grid) :CRS object has comment, which is lost in output"

####################################################
# DETERMINE MONTHLY UD (AVERAGE ACROSS IDS)
####################################################

# remove list structure
ud<-unlist(ud_id)

# create new lists for each month
## subset by month across all ids using month list name
nm <- names(ud)
all_months <- lapply(unique(nm), function(n) unname(unlist(ud[nm %in% n])))
names(all_months) <- unique(nm)

Jan<-all_months[["01"]]
Feb<-all_months[["02"]]
Mar<-all_months[["03"]]
Apr<-all_months[["04"]]
May<-all_months[["05"]]
Jun<-all_months[["06"]]
Jul<-all_months[["07"]]
Aug<-all_months[["08"]]
Sep<-all_months[["09"]]
Oct<-all_months[["10"]]
Nov<-all_months[["11"]]
Dec<-all_months[["12"]]

# stack ud layers
Jan_s <- stack(lapply(Jan, raster))
Feb_s <- stack(lapply(Feb, raster))
Mar_s <- stack(lapply(Mar, raster))
Apr_s <- stack(lapply(Apr, raster))
May_s <- stack(lapply(May, raster))
Jun_s <- stack(lapply(Jun, raster))
Jul_s <- stack(lapply(Jul, raster))
Aug_s <- stack(lapply(Aug, raster))
Sep_s <- stack(lapply(Sep, raster))
Oct_s <- stack(lapply(Oct, raster))
Nov_s <- stack(lapply(Nov, raster))
Dec_s <- stack(lapply(Dec, raster))

# find mean values for each month/raster stack
Jan_m<-mean(Jan_s)
Feb_m<-mean(Feb_s)
Mar_m<-mean(Mar_s)
Apr_m<-mean(Apr_s)
May_m<-mean(May_s)
Jun_m<-mean(Jun_s)
Jul_m<-mean(Jul_s)
Aug_m<-mean(Aug_s)
Sep_m<-mean(Sep_s)
Oct_m<-mean(Oct_s)
Nov_m<-mean(Nov_s)
Dec_m<-mean(Dec_s)

# create a list
UD_monthly<-list(Jan_m,Feb_m,Mar_m,Apr_m,May_m,Jun_m,Jul_m,Aug_m,Sep_m,Oct_m,Nov_m,Dec_m)
names(UD_monthly)<-c(1:12)

################################################################
# CONVERT UDS TO MATRIX FORM
################################################################

# convert to matrix
ud_matrix <- lapply(UD_monthly, function(x) {
  raster::as.matrix(x)
})

# normalize values to be 0-1
# this is because tracks had spatial bias (clustering) in some months due to dispersing from the nesting beach at the same time
# formula: xi-min(x)/ max(x)-min(x)
normalize <- function(x){(x-min(x))/(max(x)-min(x))}

ud_norm <- lapply(ud_matrix, function(x) {
  normalize(x)
})


############################################################
# PLOT UDS
############################################################

# plot matrices
library(plot.matrix)

effort_pal <- colorRampPalette(c("white", "red", "blue"), 
                               interpolate = 'linear')
# reg matrix
for (i in c(1:length(ud_matrix))){
  plot(ud_matrix[[i]],col=effort_pal,main=paste("Month",i),fmt.key="%.3f")
}

# normalized matrix
for (i in c(1:length(ud_norm))){
  plot(ud_norm[[i]],col=effort_pal,main=paste("Month",i),fmt.key="%.3f")
}

# NOTE: these surfaces don't look fantastic because I only have 3 individuals in my example data
## in practice, you need a large enough sample size to estimate "population-level" spatio-temporal distribution

##################################################
# SAVE DATA
##################################################

# save UD surfaces (in normalized, monthly matrix form) for use in matrix multiplication:
save(ud_norm,file="./Example_Data/UD_monthly_surfaces.rda")

##################################################
# CREATE RASTER VERSION FOR MAPPING
##################################################

# plot with background map

## Drop 0's and rasterize
ud_r<-c()
for (i in c(1:length(ud_norm))){
  # turn 0's into NA's (improves visualization and speeds up raster conversion)
  ud_norm[[i]][which(ud_norm[[i]]==0)]<-NA
  ud_r[[i]]<-raster(ud_norm[[i]])
  extent(ud_r[[i]])<-extent(-167.99,-125.99,-5.01,40.99)
  crs(ud_r[[i]])<-crs(turtle_id_sf)
}

for (i in c(1:length(ud_r))){
  # set margins
  par(mar = c(3, 3, 3, 3))
  plot(ud_r[[i]],zlim=c(0,1))
  ## The map:
  maps::map(xlim=c(-180,-100), ylim=c(-5,60), add=TRUE) # can change extent to see more/less of background map
  title(paste("Month",i))
}


