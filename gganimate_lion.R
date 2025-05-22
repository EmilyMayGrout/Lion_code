#making animations of the lions movements



#--------PARAMS-------
plot_dir <- '../results/level0/movement_tracks/'


# libraries
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(gganimate)
library(ggmap)
library(sp)
library(ggsn)
library(ggspatial)

#csv downloaded from movebank, only gps, includes 

lion <- read.csv("../data/raw/lion_highres.csv")

lion$datetime <- as.POSIXct(x = lion$timestamp, tz = 'UTC', format = "%Y-%m-%d %H:%M:%S")
lion$local_time <- with_tz(lion$datetime, tzone = "Africa/Gaborone")

unique(as.Date(lion$datetime))

#remove NA's from lon and lat
lion <- lion[!is.na(lion$location.long),]

#convert to UTM

# Create a SpatialPoints object from latitude and longitude columns
coords <- SpatialPoints(coords = lion[, c("location.long", "location.lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))

# Define the UTM CRS for Namibia (Zone 33S)
utm_crs <- CRS("+proj=utm +zone=33 +south +datum=WGS84")  # Adjust zone as needed

# Transform coordinates to UTM
coords_utm <- spTransform(coords, utm_crs)

# Extract UTM coordinates
lion$UTM_X <- coords_utm@coords[, 1]
lion$UTM_Y <- coords_utm@coords[, 2]

#want to visualise line 9 event: 28.12.21 at 11:55

#get date column
lion$date <- as_date(lion$datetime)

# Get unique dates
unique_dates <- unique(lion$date)

#-----------------------------------------------------------------------

lion_timespan <- lion[lion$timestamp >"2023-05-08 00:00:00" & lion$timestamp < "2023-05-09 00:00:00",]
lion_timespan <- lion[lion$local_time >"2023-05-08 12:00:00" & lion$local_time < "2023-05-09 12:00:00",]
point_colors <- c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#fdbf6f','#ff7f00')  # Add more colors as needed

#add google map under the points
#register_google(key="xxx")
has_google_key()

map = get_map(location = c(lon = mean(lion_timespan$location.long), lat= mean(lion_timespan$location.lat)), zoom=11, maptype="satellite")


#first make the map
g <- ggmap(map)+ 
  geom_point(data=lion_timespan, aes(x = location.long, y=location.lat, group=individual.local.identifier, colour=as.factor(individual.local.identifier)), size=5, alpha = 0.6)+
  geom_path(data=lion_timespan, aes(x = location.long, y=location.lat, group= individual.local.identifier, colour= as.factor(individual.local.identifier)))+
  scale_color_manual(values = point_colors) +
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14)) + 
  annotation_scale(location = "bl", width_hint = 0.5) + 
  coord_sf(crs = 4326)
g

g <- g + transition_reveal(datetime) +
  labs(title = 'Time: {frame_time}') +
  view_follow(fixed_y = FALSE)+
  shadow_wake(wake_length = 0.05, alpha = FALSE)

lion_vid <- animate(g, renderer = av_renderer(), height = 1000, width =1000, fps = 10, dur = 60)
lion_vid

anim_save(paste0(plot_dir, "lion.mp4"), lion_vid)


#------------------------------------------------------------------------------------

# Create the plot for each day without the map background and following the trajectories of the individuals

for (i in seq(122, length(unique_dates), by = 1)) {
  #start_date <- unique_dates[i]
  #end_date <- unique_dates[min(i+1, length(unique_dates))]
  
  #to make the video from 12:00 to 12:00 the following day
  start_date <- as.POSIXct(paste(unique_dates[i], "10:00:00"), tz = "Africa/Gaborone")
  end_date <- as.POSIXct(paste(unique_dates[min(i + 1, length(unique_dates))], "13:00:00"), tz = "Africa/Gaborone")
  
  # Filter data for the 1-day period
  # lion_timespan <- lion %>%
  #   filter(date >= start_date & date <= end_date)
  
  lion_timespan <- lion %>%
  filter(local_time >= start_date & local_time <= end_date)
  
  
  g <- ggplot() + 
    # Overall tracks with lower alpha for background
    geom_path(data = lion, aes(x = location.long, y = location.lat, 
                               group = individual.local.identifier, 
                               colour = as.factor(individual.local.identifier)), 
              alpha = 0.1) +
    # Daily tracks with higher alpha for highlighting
    geom_line(data = lion_timespan, aes(x = location.long, y = location.lat, group = individual.local.identifier, colour = as.factor(individual.local.identifier)), linewidth = 1, alpha = 0.8) +
    geom_point(data = lion_timespan, aes(x = location.long, y = location.lat,  group = individual.local.identifier,  colour = as.factor(individual.local.identifier)), 
               size = 3, alpha = 0.7) +
    scale_color_manual(values = point_colors) +
    theme(legend.title = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) + 
    annotation_scale(location = "bl", width_hint = 0.5) + 
    coord_sf(crs = 4326) +
    theme_classic()
  
  # Add animation with transition_reveal and view_follow
  g <- g + transition_reveal(lion_timespan$local_time) +
    labs(title = 'Time: {format(frame_along, "%Y-%m-%d %H:%M")}') +
    view_follow(fixed_y = FALSE) +
    shadow_wake(wake_length = 0.05, alpha = FALSE)
  
  # Animate and save the video
  lion_vid <- animate(g, renderer = av_renderer(), height = 1000, width = 1000, fps = 10, duration = 60)
  

anim_save(filename = paste0(plot_dir, "loc_time/lion_movement_loctime_", as.Date(start_date), "_to_", as.Date(end_date), "track.mp4"), animation = lion_vid)

}

#------------------------------------------------------------
#saving movement animations every 10 days and making mp4:

#continue from i = 122

#Loop through unique dates in 10-day intervals
for (i in seq(1, length(unique_dates), by = 10)) {
  start_date <- unique_dates[i]
  end_date <- unique_dates[min(i + 9, length(unique_dates))]
  
  # Filter data for the 10-day period
  lion_timespan <- lion %>%
    filter(date >= start_date & date <= end_date)
  
  # Get the map centered around the average location for the timespan
  map <- get_map(location = c(lon = mean(lion_timespan$location.long), lat = mean(lion_timespan$location.lat)), 
                 zoom = 11, maptype = "satellite")
  
  # Create the map plot
  g <- ggmap(map) + 
    geom_point(data = lion_timespan, aes(x = location.long, y = location.lat, 
                                         group = individual.local.identifier, 
                                         colour = as.factor(individual.local.identifier)), 
               size = 5, alpha = 0.6) +
    geom_path(data = lion_timespan, aes(x = location.long, y = location.lat, 
                                        group = individual.local.identifier, 
                                        colour = as.factor(individual.local.identifier))) +
    scale_color_manual(values = point_colors) +
    theme(legend.title = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) + 
    annotation_scale(location = "bl", width_hint = 0.5) + 
    coord_sf(crs = 4326)
  
  # Add animation
  g <- g + transition_reveal(datetime, keep_last = FALSE) +
    ease_aes('linear') +
    labs(title = "Time: {format(as.POSIXct(frame_along, tz = 'UTC'), '%Y-%m-%d %H:%M')}") +
    shadow_wake(wake_length = 0.05, alpha = FALSE)
  
  # Animate and save the video
  lion_vid <- animate(g, renderer = av_renderer(), height = 1000, width = 1000, fps = 10, duration = 60)
  
  anim_save(filename = paste0(plot_dir, "/lion_movement_", start_date, "_to_", end_date, ".mp4"), animation = lion_vid)
}


#-------------------------------------------------------------------------------------
