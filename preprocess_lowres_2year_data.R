#This code is reading in the Rda file Gen sent with 2 years of movement data and putting it in the matrices format for fission-fusion analysis

#LIBRARIES
library(lubridate)
library(sf)
library(tidyverse)
library(ggplot2)

#Directories
indir <-"../data/raw/gps/"
outdir <- "../../processed/"
metadatadir <-  "../data/raw/metadata/"

setwd(indir)

#load in the Rda file - dataframe called lion_df
load("C:/Users/egrout/Dropbox/lion/data/raw/lion_df_00_data_set_up.Rda")


str(lion_df)
table(lion_df$animal_id)

#plotted the raw data and can see any outliers

# Ensure animal_id is a factor or character
lion_df$animal_id <- as.character(lion_df$animal_id)

# Create unique IDs and color palette
animal_ids <- unique(lion_df$animal_id)
colors <- setNames(rainbow(length(animal_ids)), animal_ids)

# Order data by animal ID and time (if available)
lion_df <- lion_df[order(lion_df$animal_id), ]

# Plot an empty plot first
plot(NA, xlim = range(lion_df$Latitude, na.rm = T), ylim = range(lion_df$Longitude, na.rm = T), xlab = "Latitude", ylab = "Longitude")

# Loop through each individual and add their line
for (id in animal_ids) {
  inds <- lion_df[lion_df$animal_id == id, ]
  lines(inds$Latitude, inds$Longitude, col = colors[id])
}

#don't see anything that looks like an error, so won't remove any of these data
#-----------------------------------------------------------------------------

#rounding the time to the nearest 15 mins - first get it in posix format
lion_df$timestamp <- as.POSIXct(lion_df$UTC_Date_Time, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")

#plot the times when data were collected to see where we should round times too
# Extract time as decimal hours (e.g., 13.25 for 1:15 PM)
lion_df$time_of_day <- as.numeric(format(lion_df$timestamp, "%H")) +
  as.numeric(format(lion_df$timestamp, "%M")) / 60

ggplot(lion_df, aes(x = time_of_day)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  scale_x_continuous(breaks = 0:24) +
  labs(x = "Hour of Day", y = "Number of GPS Fixes",
       title = "Distribution of GPS Fixes Over the Day") +
  theme_minimal()


lion_df_filter <- lion_df[,c(28,14,13,35)]
colnames(lion_df_filter) <- c("name", "lon", "lat", "datetime")

#remove NAs
lion_df_filter <- lion_df_filter[!is.na(lion_df_filter$lon),]

min(lion_df_filter$datetime) #for firsttime
max(lion_df_filter$datetime) #for lasttime

#split by id
split_lion <- split(lion_df_filter, lion_df_filter$name, drop = F)

#save each individual as own txt file in dropbox rawdata file

allNames <- names(split_lion)

for(i in allNames){
  #get name of individual from i of list
  saveName <- paste0(i ,".txt")
  write.table(split_lion[[i]], file = saveName, sep = "%")
  
}

#now each individual has a txt file with the needed info that is saved in the gps folder in raw data


#create times vector - 24h every 15 mins

firsttime <- as.POSIXct('2023-05-04 00:00', tz = 'UTC')
lasttime <-  as.POSIXct('2025-04-12 12:00', tz = 'UTC')


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
save(list=c('xs','ys','ts'), file = 'lion_xy_15min_level0_2yr.RData')
save(list=c('lats','lons','ts'), file = 'lion_latlon_15min_level0_2yr.RData') 












