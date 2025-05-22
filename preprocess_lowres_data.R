#This code is reading in the data and putting it in the matrices format for fission-fusion analysis

#LIBRARIES
library(lubridate)
library(sf)
library(tidyverse)


#Directories
indir <-"../data/raw/gps/"
outdir <- "../data/processed/"
metadatadir <-  "../data/raw/metadata/"

setwd(indir)

#load functions
#source('C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/coati_function_library_V1.R')


#--------------------------------------------------------------------
# FORMATTING DATA FROM MOVEBANK -- don't need to run this once the txt files are made into the indir
#read in data - this is the raw CSV file from movebank
lion_all <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/lion_highres.csv") 
#when updating the code with more data, just load in the new csv file downloaded from movebank

#plotted the raw data and can see an outlier - maybe from a collar test elsewhere, so going to remove this
plot(lion_all$utm.easting, lion_all$utm.northing, col = as.character(lion_all$tag.local.identifier, type = "l"))

#removing data points where the UTM zone is 33N - in the low res, not the high res
lion_all <- lion_all[!lion_all$utm.zone == "33N",]

#plot(lion_all$utm.easting, lion_all$utm.northing, col = as.character(lion_all$tag.local.identifier, type = "l"))

#put time in posixct format
lion_all$timestamp <- as.POSIXct(lion_all$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")


#Need to read in csv and split into files for each individual and save as txt file into the indir

#setting the working directory to the raw data file in dropbox

#filter lion_all to "name", "lon", "lat", "date", "time")
lion_all_filter <- lion_all[,c(14,4,5,3)]
colnames(lion_all_filter) <- c("name", "lon", "lat", "datetime")

#remove NAs
lion_all_filter <- lion_all_filter[!is.na(lion_all_filter$lon),]

min(lion_all_filter$datetime) #for firsttime
max(lion_all_filter$datetime) #for lasttime

#split by id
split_lion <- split(lion_all_filter, lion_all_filter$name, drop = F)

#can make each individual as own dataframe in global enviro
#list2env(lapply(split_lion, as.data.frame.list), .GlobalEnv)

#save each individual as own txt file in dropbox rawdata file

allNames <- names(split_lion)

for(i in allNames){
   #get name of individual from i of list
   saveName <- paste0(i ,".txt")
   write.table(split_lion[[i]], file = saveName, sep = "%")

 }

#now each individual has a txt file with the needed info

#---------------------------------------------------------------------

#create times vector - 24h every 15 mins

firsttime <- as.POSIXct('2023-05-03 05:00', tz = 'UTC')
lasttime <-  as.POSIXct('2023-12-03 11:45', tz = 'UTC')


#first make list of 10 min intervals then remove the times overnight
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '15 min')

#making a list of all individuals in presedente
all_files <- sort(list.files())

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))


for(i in 1:length(all_files)){
  
  #filter columns to ones needed
  tagdata <- read.table(all_files[i], sep = "%")
  #getting the correct dat and time column
  tagdata <- tagdata[, c(2, 3, 4)]
  colnames(tagdata) <- c("lon", "lat", "datetime")
  #make datetime format
  tagdata$datetime <- as.POSIXct(tagdata$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  #round to the nearest 15 minutes
  tagdata$datetime <- round_date(tagdata$datetime, "15 minutes")
  
  
  #getting the last gps point per burst
  tagdata$test <- !duplicated(tagdata$datetime, fromLast = T )
  
  #removing duplicates
  tagdata <- tagdata[tagdata$test == T,]
  
  #remove rows with 0's from the lat and lon
  tagdata <- tagdata[which(tagdata$lon != 0),]
  
  #match times to get lons and lats at each time for that individual
  lon <- tagdata$lon[match(ts, tagdata$datetime)]
  lat <- tagdata$lat[match(ts, tagdata$datetime)]
  
  lats[i,] <- lat
  lons[i,] <- lon
  
  #to convert to UTM, need to use sf function and for this we need to combine the matrices rows to a dataframe
  combined_df <- data.frame(lat = lats[i,], lon = lons[i,])
  
  #convert NA's to zero for sf functions to work
  combined_df$lon[is.na(combined_df$lon)] <- 0
  combined_df$lat[is.na(combined_df$lat)] <- 0
  
  
  #convert to UTM
  #first need to give the latlon data the correct CRS so it converts to UTM correctly
  latlon_i <- st_as_sf(x=combined_df, coords=c("lon", "lat"), crs=4326)
  #convert to UTM
  utm_i <- st_transform(latlon_i, crs="+proj=utm +zone=33 +south +datum=WGS84 +units=m") #convert UTM to lat/long - coordinates are stored in a geometry column of class 'sfc_POINT'
  
  #store eastings and northings in xs and ys matrices
  xs[i,] <- unlist(map(utm_i$geometry,1))
  ys[i,] <- unlist(map(utm_i$geometry,2))
  
  
}

#removing the 0's which were converved to eastings and northings and should be replaced with NA
xs[xs < 0] <- NA
ys[ys == 10000000] <- NA


setwd(outdir)
save(list=c('xs','ys','ts'), file = paste0(outdir,'lion_xy_15min_level0.RData'))
save(list=c('lats','lons','ts'), file = paste0(outdir,'lion_latlon_15min_level0.RData'))  

#TODO metadata
setwd(metadatadir)
lion_ids <- read.csv("lion_ids.csv", header = F)
colnames(lion_ids) <- c("name", "tag_id", "age", "sex")
lion_ids$color <- '#0000FF'
lion_ids$color[which(lion_ids$age == 'Adult' & lion_ids$sex == 'Female')] <- '#FF0000'
lion_ids$color[which(lion_ids$age == 'Sub-adult' & lion_ids$sex == 'Female')] <- '#FFAA66'
lion_ids$color[which(lion_ids$age == 'Sub-adult' & lion_ids$sex == 'Male')] <- '#66AAFF'
lion_ids$color[which(lion_ids$age == 'Juvenile')] <- '#666666'
save(lion_ids, file = paste0(outdir, 'lion_ids.RData'))



