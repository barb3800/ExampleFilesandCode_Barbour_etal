## EXAMPLE: EXTRACTING THREAT DATA (GLOBAL FISHING WATCH)

# This data comes from Global Fishing Watch (https://globalfishingwatch.org/) and is free to download
# The product we used was their 100th degree (0.01 deg/1 km) AIS Fishing effort, with data on effort by region and geartype
# We wanted to aggregate the data across 3 years (2018, 2019, 2020)
# We wanted the spatial resolution to be 1 deg, since our study area was quite large 

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

# NOTE: the data files you download from Global Fishing Watch can be MASSIVE and cannot be stored on Github
## The first portion of code here is just for example (lines 32-95) but does not need to be run if you aren't getting your own data

##############################################
# SET UP 
##############################################

# load libraries (note: all these packages need to be installed before use)
library(tidyverse);library(furrr);library(lubridate);library(sf);library(raster);library(maps);library(maptools);library(rgeos);library(dplyr);library(rgdal)

# set working directory
setwd("yourfilepathhere")

# read in data (.rda format)
## turtle tracking data
load("./Example_Data/Jitter_Turtle_Data.rda")

#############################################################
# STEP 1: BRING DATA INTO R
#############################################################

# You will need to make an account with GFW and download the correct product
# The data comes in as bunch of csv files (daily) for each year
# The code below should be run for each year

# NOTE: this data is massive (there is a csv file for each day of the year) and can take a LONG time to process!
## I recommend running on an external server or a computer with lots of RAM/memory

# example code below for 2020 (run also for other years of interest)

# specify data directory location (where you downloaded the data from GFW)
dir<-"D:/home/nbarbour/fleet-daily-csvs-100-v2-2020/"

# create dataframe of file name dates to filter by date range
# using 100th deg effort, daily
effort<-tibble(file=list.files(paste0(dir,"fleet-daily-csvs-100-v2-2020"),pattern=".csv",recursive=TRUE,full.names=TRUE),date=ymd(str_extract(file,pattern = '[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')))

# create vector of dates of interest
effort_dates<-seq(ymd("2020-01-01"),ymd("2020-12-31"),by="days")

# filter effort data to dates
effort_2020<-filter(effort, date %in% effort_dates)

# read in data
plan(multisession) #multiprocess for mac users
## for personal comp, may need to change memory ability
memory.limit()
#memory.limit(size=18000)

effort_df<-future_map_dfr(effort_2020$file,.f=read_csv)

#add date info
effort_df_2020<-effort_df %>% mutate(year=year(date),month=month(date))

###############################################
# STEP 2: PROCESS DATA
###############################################

# The data needs to be cropped to the boundaries of the study region & converted into an sf object

# crop to study region boundaries
## use WGS 1984 crs (same as polygons, EPSG code 4326)
effort_2020_s<-subset(effort_df_2020,cell_ll_lon > -168 & cell_ll_lon < -126 & cell_ll_lat > -5 & cell_ll_lat < 41 )
effort_2019_s<-subset(effort_df_2019,cell_ll_lon > -168 & cell_ll_lon < -126 & cell_ll_lat > -5 & cell_ll_lat < 41 )
effort_2018_s<-subset(effort_df_2018,cell_ll_lon > -168 & cell_ll_lon < -126 & cell_ll_lat > -5 & cell_ll_lat < 41 )

# make month a factor
effort_2020_s$month<-as.character(effort_2020_s$month)
effort_2019_s$month<-as.character(effort_2019_s$month)
effort_2018_s$month<-as.character(effort_2018_s$month)

# change gear type level names to look nicer (can check current with "levels(effort_2020_s$geartype)")
effort_2020_s$geartype<-as.factor(effort_2020_s$geartype)
effort_2019_s$geartype<-as.factor(effort_2019_s$geartype)
effort_2018_s$geartype<-as.factor(effort_2018_s$geartype)

levels(effort_2020_s$geartype)<-c("Drifting Longline", "Fishing", "Pole and Line","Pots and Traps","Purse Seines","Squid Jigger","Trawlers","Trollers","Tuna Purse Seines")
levels(effort_2019_s$geartype)<-c("Drifting Longline", "Fishing","Fixed Gear", "Pole and Line","Purse Seines","Set Longlines","Squid Jigger","Trawlers","Trollers","Tuna Purse Seines")
levels(effort_2018_s$geartype)<-c("Drifting Longline", "Fishing","Squid Jigger","Trawlers","Trollers","Tuna Purse Seines")

# save data
save(effort_2020_s,file="./Example_Data/2020_GFW_Data.rda")
save(effort_2019_s,file="./Example_Data/2019_GFW_Data.rda")
save(effort_2018_s,file="./Example_Data/2018_GFW_Data.rda")

##########################################################
# STEP 3: AGGREGATE ACROSS YEARS
#########################################################

# NOTE: here you can load in the intermediate data, that has been subset to the study region and can be used to run the rest of the example code
# load in data
load(file="./Example_Data/2020_GFW_Data.rda")
load(file="./Example_Data/2019_GFW_Data.rda")
load(file="./Example_Data/2018_GFW_Data.rda")

# Next we aggregate the data across years (2018-2020) by finding the total effort (sum) across years in a given lat/long cell 

# combine years and sum effort in each cell
gfw_allyrs<-bind_rows(effort_2020_s,effort_2019_s,effort_2018_s) %>%
  group_by(cell_ll_lon,cell_ll_lat, month,geartype) %>% # group by lat/lon and month
  summarize(fishing_hours = sum(fishing_hours, na.rm = T)) %>%
  ungroup() %>%
  data.frame()

# drop rows with fishing effort =0
gfw_allyrs2<-subset(gfw_allyrs,fishing_hours > 0)

# make data into an sf object
effort_all_sf<- gfw_allyrs2 %>% st_as_sf(coords=c("cell_ll_lon","cell_ll_lat"),crs=4326)

# factorize month column
effort_all_sf$month<-as.factor(effort_all_sf$month)

# reorganize levels to be month 1-12
effort_all_sf$month<-factor(effort_all_sf$month,levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))

# drop any unused levels for geartype (not all geartypes will be in your region)
effort_all_sf$geartype<-droplevels(effort_all_sf$geartype)

############################################################
# STEP 4: BREAK UP BY GEARYTYPE & MONTH
############################################################

# We are interested in understanding the overlap between individual fishing gear types, turtle behaviors, and turlte spatio-density predicitons, on a monthly resolution
## there are lots of geartypes for GFW data...since drifting longlines are impactful for leatherback survival, we use this geartype as example

## subset for drifting longline geartypes (g1)
g1<-subset(effort_all_sf,geartype=="Drifting Longline")

# split up by month (list of 12 months of fishing effort for drifting longlines)
g1_m<-split(g1,g1$month)

###########################################################
# STEP 5: RASTERIZE & AGGREGATE TO COARSER RESOLUTION
##########################################################

# next we need to convert the GFW data into a raster format
## and we need to convert it to a coarser resolution (1 deg), since 0.01 deg is too fine to capture multiple days of behavior for this population
## you should choose a resolution suitable for the scale of behavior of your data!

# rasterize using base raster set to res and correct extent
base<-raster()
extent(base)<-extent(effort_all_sf) # set extent to same as study region
crs(base)
res(base)<-c(0.01,0.01) # 0.01 deg resolution to match GFW data

# sum fishing hours in each cell
g1_r <- lapply(g1_m, function(x) {
  raster::rasterize(x=x,y=base,field="fishing_hours",fun=sum,background=0)
})

# aggregate raster to coarser resolution- 1 x 1 deg (factor of 100)
# sum fishing effort in new cells
g1_agg <- lapply(g1_r, function(x) {
  aggregate(x,fact=100,fun=sum)
})

# save
save(g1_agg,file="./Example_Data/longlinefishing_raster.rda")






