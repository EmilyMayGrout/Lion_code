#making animations of the lions movements


#--------PARAMS-------
plot_dir <- 'C:/Users/egrout/Dropbox/lion/results/level0/'


# libraries
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(gganimate)
library(ggmap)
library(sp)
library(ggsn)

#csv downloaded from movebank, only gps, includes 

lion <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/lion_highres.csv")

lion$datetime <- as.POSIXct(x = lion$timestamp, tz = 'UTC', format = "%Y-%m-%d %H:%M:%S")

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

lion_one_day <- lion[lion$timestamp >"2023-05-08 00:00:00" & lion$timestamp < "2023-05-15 00:00:00",]

point_colors <- c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#fdbf6f','#ff7f00')  # Add more colors as needed

#add google map under the points
register_google(key="xxx")
has_google_key()


map = get_map(location = c(lon = mean(lion_one_day$location.long), lat= mean(lion_one_day$location.lat)), zoom=11, maptype="satellite")

#first make the map
g <- ggmap(map)+ 
  geom_point(data=lion_one_day, aes(x = location.long, y=location.lat, group=individual.local.identifier, colour=as.factor(individual.local.identifier)), size=5)+
  #geom_path(data=lion_one_day, aes(x = location.long, y=location.lat, group= individual.local.identifier, colour= as.factor(individual.local.identifier)))+
  scale_color_manual(values = point_colors) +
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_text(size=14), axis.title = element_text(size = 14), legend.text = element_text(size = 14))+annotation_scale(location = "bl", width_hint = 0.5) + 
  coord_sf(crs = 4326)
g


g <- g + transition_reveal(lion_one_day$datetime, keep_last = FALSE)+
  ease_aes('linear')+
  labs(title = "Time: {format(as.POSIXct(frame_along, tz = 'UTC'), '%Y-%m-%d %H:%M:%S')}")+
  shadow_wake(wake_length = 0.05, alpha = FALSE)

animate(g, height = 1000, width =1000, fps = 10, dur = 60)

anim_save(paste0(plot_dir, "test.gif"))






