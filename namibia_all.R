
source('coati_function_library_V1.R')
source('lion_functions.R')


library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(sf)
library(tidyverse)
library(fields)
library(viridis)
library(tidygraph)

#load data downloaded from movebank 
all_nab <- read.csv("../data/raw/all_namibia.csv")

#Directories
indir <-"../data/raw/gps/all_nab/"
outdir <- "../data/processed/"
metadatadir <- "../data/raw/metadata/"

setwd(indir)


#looking at which inds have data for each month for each year


#put time in posixct format
all_nab$timestamp <- as.POSIXct(all_nab$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")
all_nab$month_year <- format(all_nab$timestamp, "%Y-%m")

presence_matrix <- all_nab %>%
  group_by(individual.local.identifier, month_year) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = month_year, values_from = present, values_fill = 0)

# Convert to long format for ggplot
presence_long <- presence_matrix %>%
  pivot_longer(-individual.local.identifier, names_to = "month_year", values_to = "present")

# Plot heatmap
gg <- ggplot(presence_long, aes(x = month_year, y = individual.local.identifier, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  theme_classic() +
  labs(x = "Month-Year", y = "Individual", fill = "Data Present") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg

ggsave(filename = "C:/Users/egrout/Dropbox/lion/results/level0/data_presence_all_namib.png", plot = gg, width = 15, height = 8, dpi = 300)

count_matrix <- all_nab %>%
  group_by(individual.local.identifier, month_year) %>%
  summarise(count = n(), .groups = "drop")

count_matrix$month_year <- factor(count_matrix$month_year, levels = sort(unique(count_matrix$month_year)))

ggplot(count_matrix, aes(x = month_year, y = individual.local.identifier, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey90") +
  theme_classic() +
  labs(
    title = "Data Point Counts per Individual per Month-Year",
    x = "Month-Year",
    y = "Individual",
    fill = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#-----------------------------------------------

#look at the proportion of time within 500m at 10am (UTC time), midday local time

#preprocess the data to make xs and ys matrices

##ONLY NEED TO RUN ONCE
# #filter lion_all to "name", "lon", "lat", "date", "time")
# nab_filt <- all_nab[,c(9,4,5,3)]
# colnames(nab_filt) <- c("name", "lon", "lat", "datetime")
# 
# #remove NAs
# nab_filt <- nab_filt[!is.na(nab_filt$lon),]
# 
# min(nab_filt$datetime) #for firsttime
# max(nab_filt$datetime) #for lasttime
# 
# 
# #split by id
# split_lion <- split(nab_filt, nab_filt$name, drop = F)
# 
# #save each individual as own txt file in dropbox rawdata file
# 
# allNames <- names(split_lion)
# 
# for(i in allNames){
#   #get name of individual from i of list
#   saveName <- paste0(i ,".txt")
#   write.table(split_lion[[i]], file = saveName, sep = "%")
#   
# }

#should remove a few days at start - so start from 6th May 
firsttime <- as.POSIXct('2020-06-07 05:00', tz = 'UTC')
lasttime  <- as.POSIXct('2025-04-03 11:00', tz = 'UTC')


#first make list of day intervals then remove the times overnight
dates <- seq.Date(from = as.Date(firsttime), to = as.Date(lasttime), by = "day")
ts <- as.POSIXct(paste(dates, "10:00:00"), tz = "UTC")

#making a list of all individuals in presedente
all_files <- sort(list.files())

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))

for(i in 1:length(all_files)){
  
  #filter columns to ones needed
  tagdata <- read.table(all_files[i], sep = "%")
  #getting the correct date and time column
  tagdata <- tagdata[, c(2, 3, 4)]
  colnames(tagdata) <- c("lon", "lat", "datetime")
  #make datetime format
  tagdata$datetime <- as.POSIXct(tagdata$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  
  #round to the nearest hour
  tagdata$datetime <- round_date(tagdata$datetime, "hour")
  
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


#make lion ids 
lion_ids <- data.frame(name = unique(all_nab$individual.local.identifier))
lion_ids$name <- as.character(lion_ids$name)

all_ids <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/metadata/namib_lion_ids.csv", header = F)



lion_ids$name[!(lion_ids$name %in% all_ids$V1)]
#missing "OPL-24" "NPL-28" "NPL-27" "NPL-40" from the GPS data


#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#subset to times when the data is every 15 mins
non_na_counts <- colSums(!is.na(xs))

R = 500
subgroup_data <- get_subgroup_data(xs, ys, R)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

split_by_year <- T
year <- unique(format(ts, "%Y"))
ff_net_year <- array(NA, dim = c(n_inds, n_inds, length(year)))
dimnames(ff_net_year)[[3]] <- year

for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    if(split_by_year){
      for(y in 1:length(seq_along(year))){
        
        y <- year[y]
        
        # Logical index for matching year-month
        ts_filter <- format(ts, "%Y") == y
        
        sub_ids_i <- subgroup_data$ind_subgroup_membership[i, ts_filter]
        sub_ids_j <- subgroup_data$ind_subgroup_membership[j, ts_filter]
        
        ff_net_year[i, j, y] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
        
      }
    }
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
order <- 1:length(n_inds)
ffnet_reorder <- ff_net[order, order]

visualize_lion_network(ff_net, lion_ids[order,])
#need to get lion_ids in right format!!








