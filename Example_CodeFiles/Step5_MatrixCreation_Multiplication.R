## Code to Create Matrices and Perform OVerlap Analysis

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

##############################################
# SET UP 
##############################################

# load packages (make sure you have installed them first!)
library(sf);library(dplyr);library(raster);

# set working directory
setwd("yourfilepathhere")

## load in data
# 3D UD surfaces (monthly)
load(file="./Example_Data/UD_monthly_surfaces.rda")
# GFW fishing surface (monthly, drifting longline fishing)
load("./Example_Data/longlinefishing_raster.rda")
# jittered tracking data with HMM predicted states
load(file="./Example_Data/turtle_jitter_states.rda")

##################################################
# STEP 1: RASTERIZE BEHAVIOR DATA
##################################################

# currently, our behavior data is set up as having HMM predicted states for each observed day of tracking for each each id
# we want to convert this behavior data to a surface
## we want a surface for each behavior (S1-S3) and aggregated across individuals
## we will further break these surfaces up by month, to match the temporal res of our threat (fishing) and UD data surfaces

## convert states to spatial points (sf)
turtle_sf<-turtle_jitter_states %>% st_as_sf (coords = c("lon_j","lat_j"), crs = 4326)

# subset for each state
## as example, we will just focus on behavior state 3 (deeper dives, intermediate persistence velocity)
turtle_s3<-subset(turtle_sf,States=="3")

# note: if there are missing months, you need to add in some data ( we will do this below)

# drop  unused levels 
## this is needed for factors when you subset, since unused or empty levels will "stick" afterwards
turtle_s3$Month<-droplevels(turtle_s3$Month)

# split up by month
## note: if you don't have 12 months of data you won't have 12 lists as a result
turtle_s3_month<-split(turtle_s3,turtle_s3$Month)

names(turtle_s3_month)
# this data is missing month 1 (Jan)

# rasterize using states variable and counting number of points in each cell
## NOTE: it is important that all of your data have the same spatial extent and resolution!
## use extent and res same as effort raster for background raster
base<-raster()
extent(base)<-extent(g1_agg[[1]]) # set extent to same as fishing data
crs(base)<-crs(g1_agg[[1]])
res(base)<-c(1,1) # 1 deg resolution

# create list for monthly rasters of behavior:
turtle_s3_raster<-c()

# for loop for each month
# rasterize State 3 behavior for each month by counting the number of times that behavior occurs in a particular lat/long cell
for (i in c(1:length(turtle_s3_month))){
  turtle_s3_raster[[i]]<-raster::rasterize(turtle_s3_month[[i]],base,field="States",fun="count",background=0)
}

##########################################################
# STEP 2: CREATE MATRICES
##########################################################

# Note: the 3D UD surfaces are already in matrix form (list of matrices, 1 for each month)

## GFW data
g1_ma <- lapply(g1_agg, function(x) {
  raster::as.matrix(x)
})

## HMM data
turtle_s3_matrix<-c()
for (i in c(1:length(turtle_s3_raster))){
  turtle_s3_matrix[[i]]<-raster::as.matrix(turtle_s3_raster[[i]])
}

##########################################################
# STEP 3: NORMALIZE MATRICES
##########################################################

# in order to effectively perform elementwise matrix multiplication, values for each matrix need to be on the same scale
## we scale the matrices (monthly for each fishing, UD, and HMM surface) to be from 1-5 in value
## note: we use values 1-5 to prevent 0's in one matrix from excluding that layer when multiplied

library(scales)

g1_n <- lapply(g1_ma, function(x) {
  rescale(x,to=c(1,5))
})

turtle_s3_matrix_n <- lapply(turtle_s3_matrix, function(x) {
  rescale(x,to=c(1,5))
})

ud_n <- lapply(ud_norm, function(x) {
  rescale(x,to=c(1,5))
})

# ADD IN MISSING MONTH DATA
## to fill in the missing month of data (Jan), we can copy one of our current matrices, change all the values to 1, and bind to our list of scaled behavior matrices
month1<-turtle_s3_matrix_n[[1]]
month1[month1 > 1] <-1 
turtle_s3_matrix_n<-append(list(month1),turtle_s3_matrix_n)

##########################################################
# STEP 4: ELEMENT-WISE MULTIPLICATION OF MATRICES
##########################################################

# for each month, multiply the corresponding UD, GFW, and HMM matrix surfaces 
## NOTE: make sure your matrices have the same dimensions! If they don't, you won't be able to multiply them and it's likely because they don't have the same spatial extent/resolution
matrix_g1_s3<-c()

for (i in c(1:length(turtle_s3_matrix_n))){
  matrix_g1_s3[[i]] <- g1_n[[i]]  * turtle_s3_matrix_n[[i]] * ud_n[[i]]
}

##########################################################
# STEP 5: DETERMINE RELATIVE OVERLAP OF MATRICES
##########################################################

# to determine the relative overlap or interaction between these 3 matrices for each month, we divide each monthly matrix by the sum of its values

g1_s3_d<-c()
for (i in c(1:length(matrix_g1_s3))){
  g1_s3_d[[i]] <- matrix_g1_s3[[i]]/sum(matrix_g1_s3[[i]])
}

# save relative overlap matrices for later
save(g1_s3_d,file="./Example_Data/relative_overlap_matrices.rda")

##########################################################
# STEP 6: RASTERIZE FINAL MATRIX
##########################################################

# finally we convert our final matrix into a raster for ease in visualization and mapping

## raster
g1_s3_r<-c()
for (i in c(1:length(g1_s3_d))){
  g1_s3_r[[i]]<-raster(g1_s3_d[[i]])
  extent(g1_s3_r[[i]])<-extent(g1_agg[[1]])
  crs(g1_s3_r[[i]])<-crs(g1_agg[[1]])
}

#######################################################
# STEP 7: VISUALIZE AND SAVE
######################################################

# plot matrices
library(plot.matrix)

effort_pal <- colorRampPalette(c("blue", "red", "white"), 
                               interpolate = 'linear')

for (i in c(1:length(g1_s3_d))){
  plot(g1_s3_d[[i]],col=effort_pal,main=paste("Month",i),fmt.key="%.3f")
}

# plot rasters:
for (i in c(1:length(g1_s3_r))){
  plot(g1_s3_r[[i]],main=paste("Month",i))
}

## NOTE: again, because we are only using 3 tracks in this example data, the overlap matrices don't look great- there are few places of overlap due to the sparcity of behavior 3 in the HMM layers
# month 1 looks the best because there is no behavior to overlap with (it's only showing the overlap between fishing and UD distribution in Jan)

save(g1_s3_r,file="./Example_Data/combined_monthly_overlap.rda")


