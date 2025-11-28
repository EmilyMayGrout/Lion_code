#this script is making the xs, ys, and ts matrix for all namibia GPS data

#LIBRARIES
library(lubridate)
library(sf)
library(tidyverse)
library(hms)
library(patchwork)
library(reshape2)
library(ggplot2)



setwd("C:/Users/egrout/Dropbox/lion/data/raw/gps/all_nab/")
plot_dir <- 'C:/Users/egrout/Dropbox/lion/results/level0/'

#load data downloaded from movebank 
all_nab <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/all_namibia.csv")
all_ids <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/metadata/namib_lion_ids.csv", header = F)
colnames(all_ids) <- c("name", "sex", "group_type", "pride_id", "yob")

#put time in posixct format
all_nab$timestamp <- as.POSIXct(all_nab$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")
all_nab$month_year <- format(all_nab$timestamp, "%Y-%m")

all_nab$time_of_day <- as.numeric(format(all_nab$timestamp, "%H")) +
  as.numeric(format(all_nab$timestamp, "%M")) / 60

tag_summary <- data.frame(table(all_nab$tag.local.identifier, all_nab$individual.local.identifier))


g1 <- ggplot(all_nab, aes(x = as_hms(time_of_day))) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "blue") +
  scale_x_continuous(breaks = 0:24) +
  labs(x = "Hour of Day", y = "Number of GPS Fixes",
       title = "Count of GPS Fixes For each hour the Day") +
  facet_wrap(~individual.local.identifier, ncol = 3)+
  theme_classic(base_size = 12)

g1
ggsave(filename = paste0(plot_dir, "histdatacount_eachhour_perind.png"),
       plot = g1, width = 16, height = 22, dpi = 300)

#does each lion have one collar? - no, some have 2, so need to find out which collars are the ones we want to use for analysis
# make a lookup table for IDs
id_lookup <- unique(all_nab[, c("individual.local.identifier", "tag.local.identifier")])

# join tag info to dataset
all_nab2 <- all_nab %>%
  left_join(id_lookup, by = "individual.local.identifier")

# build combined label
all_nab2$facet_label <- paste0(all_nab2$individual.local.identifier, 
                               " (Tag: ", all_nab2$tag.local.identifier.x, ")")

# histogram with facet labels
g2 <- ggplot(all_nab2, aes(x = as_hms(time_of_day))) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "blue") +  # 3600 sec = 1 hr bins
  scale_x_continuous(breaks = 0:24) +
  labs(x = "Hour of Day", y = "Number of GPS Fixes",
       title = "Count of GPS Fixes by Hour of Day per Individual (with Tag ID)") +
  facet_wrap(~facet_label, ncol = 4, scales = "free_x", strip.position = "top"   # put axis labels under each facet
  ) +
  theme_classic(base_size = 12)
g2

ggsave(filename = paste0(plot_dir, "histdatacount_eachhour_perindpertag.png"),
       plot = g2, width = 20, height = 28, dpi = 300)


all_nab_filter <- all_nab[,c(9,4,5,3)]
colnames(all_nab_filter) <- c("name", "lon", "lat", "datetime")

#remove NAs
all_nab_filter <- all_nab_filter[!is.na(all_nab_filter$lon),]

min(all_nab_filter$datetime) #for firsttime
max(all_nab_filter$datetime) #for lasttime

presence_matrix <- all_nab %>%
  dplyr::group_by(individual.local.identifier, month_year) %>%
  dplyr::summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = month_year, values_from = present, values_fill = 0)

# Convert to long format for ggplot
presence_long <- presence_matrix %>%
  pivot_longer(-individual.local.identifier, names_to = "month_year", values_to = "present")

presence_long$month_year <-  as.Date(paste0(presence_long$month_year, "-01"), format = "%Y-%m-%d")

# Plot heatmap
gg <- ggplot(presence_long, aes(x = month_year, y = individual.local.identifier, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months")+
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  theme_classic() +
  labs(x = "Month-Year", y = "Individual", fill = "Data Present", title = "Raw Data Presence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg

#--------------------------------
#when is each TAG collecting data?
presence_matrix2 <- all_nab %>%
  group_by(tag.local.identifier, month_year) %>%
  dplyr::summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = month_year, values_from = present, values_fill = 0)

# Convert to long format for ggplot
presence_long2 <- presence_matrix2 %>%
  pivot_longer(-tag.local.identifier, names_to = "month_year", values_to = "present")

g3 <- ggplot(presence_long2, aes(x = month_year, y = as.factor(tag.local.identifier), fill = factor(present))) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "slateblue")) +
  theme_classic() +
  labs(x = "Month-Year", y = "Tag ID", fill = "Data Present") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 9))
  
g3

ggsave(filename = paste0(plot_dir, "datapertag.png"),
       plot = g3, width = 15, height = 8, dpi = 300)

#--------------------------------

#split by id
split_lion <- split(all_nab_filter, all_nab_filter$name, drop = F)
# Reorder list to match order in all_ids$name
split_lion <- split_lion[all_ids$name]

#save each individual as own txt file in dropbox rawdata file
allNames <- names(split_lion)

# Save each individual in the right order
for (i in seq_along(all_ids$name)) {
  ind_name <- all_ids$name[i]
  saveName <- paste0(ind_name, ".txt")
  write.table(split_lion[[ind_name]], file = saveName, sep = "%", row.names = FALSE)
}

#now each individual has a txt file with the needed info that is saved in the gps folder in raw data

#create times vector - 24h every hour - because not all on the same 2h schedule
firsttime <- as.POSIXct('2020-06-07 00:00', tz = 'UTC')
lasttime <-  as.POSIXct('2025-04-03 12:00', tz = 'UTC')
#first make list of intervals
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '2 hours')


#making a list of all individuals
all_files <- sort(list.files())

# Reorder files to match the order in all_ids$name
all_files <- all_files[match(paste0(all_ids$name, ".txt"), all_files)]

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))

for(i in 1:length(all_files)){
  
  #filter columns to ones needed
  tagdata <- read.table(all_files[i], sep = "%")
  #remove first row with titles
  tagdata <- tagdata[-1,]
  
  #getting the correct date and time column
  tagdata <- tagdata[, c(2, 3, 4)]
  colnames(tagdata) <- c("lon", "lat", "datetime")
  #make datetime format
  tagdata$datetime <- as.POSIXct(tagdata$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  
  #round to the nearest hour
  # tagdata$datetime <- round_date(tagdata$datetime, "2 hours")
  # 
  # #getting the last gps point per burst
  # tagdata$test <- !duplicated(tagdata$datetime, fromLast = T )
  # 
  # #removing duplicates
  # tagdata <- tagdata[tagdata$test == T,]
  # 
  # define target hours
  target_hours <- c(6, 10, 14, 18, 20, 22, 0, 2, 4)
  
  tagdata <- tagdata %>%
    mutate(date = as.Date(datetime))
  
  # expand all combinations of date × target hour
  grid <- expand_grid(
    date = unique(tagdata$date),
    target = target_hours
  ) %>%
    mutate(target_time = as.POSIXct(date) + hours(target))
  
  # join with actual fixes and compute time differences
  matched <- grid %>%
    left_join(tagdata, by = "date") %>%
    mutate(diff_mins = abs(as.numeric(difftime(datetime, target_time, units = "mins"))))
  
  # keep closest fix within ±10 min
  final <- matched %>%
    filter(diff_mins <= 20) %>%
    group_by(date, target) %>%
    slice_min(diff_mins, with_ties = FALSE) %>%
    ungroup()
  
  # ensure all expected rows are present (fill with NA if no fix)
  final_complete <- grid %>%
    left_join(final, by = c("date", "target", "target_time"))
  
  
  #remove rows with 0's from the lat and lon
  #tagdata <- tagdata[which(tagdata$lon != 0),]
  
  #match times to get lons and lats at each time for that individual
  #lon <- tagdata$lon[match(ts, tagdata$datetime)]
  #lat <- tagdata$lat[match(ts, tagdata$datetime)]
  
  lon <- final_complete$lon[match(ts, final_complete$target_time)]
  lat <- final_complete$lat[match(ts, final_complete$target_time)]
  
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


#removing the 0's which were converted to eastings and northings and should be replaced with NA
xs[xs < 0] <- NA
ys[ys == 10000000] <- NA

month_year <- format(ts, "%Y-%m")  # vector same length as ncol(xs)
presence_raw <- ifelse(!is.na(xs) | !is.na(ys), 1, 0)
unique_months <- sort(unique(month_year))
presence_month <- matrix(0, nrow = nrow(xs), ncol = length(unique_months))
rownames(presence_month) <- all_ids$name
colnames(presence_month) <- unique_months

for (i in seq_along(unique_months)) {
  cols_in_month <- which(month_year == unique_months[i])
  presence_month[, i] <- apply(presence_raw[, cols_in_month, drop = FALSE], 1, max)
}
presence_long <- melt(presence_month)
colnames(presence_long) <- c("individual", "month", "present")
ggplot(presence_long, aes(x = month, y = individual, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "grey85", "1" = "steelblue")) +
  theme_minimal() +
  labs(x = "Month", y = "Individual", fill = "Data present") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
save(list=c('xs','ys','ts'), file = 'allnamib_xy_2hour_level1.RData')
save(list=c('lats','lons','ts'), file = 'allnamib_latlon_2hour_level1.RData') 


#put age categories by year of birth

all_ids$sex[all_ids$sex == "F"] <- "Female"
all_ids$sex[all_ids$sex == "M"] <- "Male"

all_ids$color <- '#0000FF'
all_ids$color[which(all_ids$yob == '2012' & all_ids$sex == 'Female')] <- '#8B0000'
all_ids$color[which(all_ids$yob == '2013' & all_ids$sex == 'Female')] <- '#A31717'
all_ids$color[which(all_ids$yob == '2014' & all_ids$sex == 'Female')] <- '#B92222'
all_ids$color[which(all_ids$yob == '2015' & all_ids$sex == 'Female')] <- '#CF2D2D'
all_ids$color[which(all_ids$yob == '2016' & all_ids$sex == 'Female')] <- '#E53838'
all_ids$color[which(all_ids$yob == '2017' & all_ids$sex == 'Female')] <- '#EB6262'
all_ids$color[which(all_ids$yob == '2018' & all_ids$sex == 'Female')] <- '#F18C8C'
all_ids$color[which(all_ids$yob == '2019' & all_ids$sex == 'Female')] <- '#F7B6B6'
all_ids$color[which(all_ids$yob == '2020' & all_ids$sex == 'Female')] <- '#FDE0E0'
all_ids$color[which(all_ids$yob == '2012' & all_ids$sex == 'Male')] <- '#00441B'
all_ids$color[which(all_ids$yob == '2013' & all_ids$sex == 'Male')] <- '#1B5D31'
all_ids$color[which(all_ids$yob == '2014' & all_ids$sex == 'Male')] <- '#377548'
all_ids$color[which(all_ids$yob == '2015' & all_ids$sex == 'Male')] <- '#538F61'
all_ids$color[which(all_ids$yob == '2016' & all_ids$sex == 'Male')] <- '#70A97B'
all_ids$color[which(all_ids$yob == '2017' & all_ids$sex == 'Male')] <- '#93C399'
all_ids$color[which(all_ids$yob == '2018' & all_ids$sex == 'Male')] <- '#B7DDBA'
all_ids$color[which(all_ids$yob == '2019' & all_ids$sex == 'Male')] <- '#D5F0D7'
all_ids$color[which(all_ids$yob == '2020' & all_ids$sex == 'Male')] <- '#E5F5E0'

save(all_ids, file = 'allnamib_ids.RData')



#checking whether data collection matches original data plot
ind_names <- all_ids$name

# Convert to long format
df <- as.data.frame(xs) %>%
  mutate(individual = ind_names) %>%
  pivot_longer(
    cols = -individual,
    names_to = "time_index",
    values_to = "xs"
  ) %>%
  mutate(
    time_index = as.numeric(sub("V", "", time_index)),  # remove 'V' and convert to numeric
    datetime = ts[time_index]
  ) %>%
  select(individual, datetime, xs)


df$month_year <- floor_date(df$datetime, "month")

# Count non-NA fixes per month-year per individual
coverage_df <- df %>%
  group_by(individual, month_year) %>%
  dplyr::summarise(n_fixes = sum(!is.na(xs)), 
            .groups = "drop")

coverage_df$month_year <- as.Date(coverage_df$month_year) 

# Plot
gg1 <- ggplot(coverage_df, aes(x = month_year, y = individual, fill = n_fixes)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",    
    high = "purple"  
  )+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Month-Year",
    y = "Individual",
    fill = "Fixes",
    title = "GPS Data Coverage per Individual per Month"
  )

coverage_df <- coverage_df %>%
  mutate(present = ifelse(n_fixes > 0, 1, 0))

gg2 <- ggplot(coverage_df, aes(x = month_year, y = individual, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "purple")) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Month-Year",
    y = "Individual",
    fill = "Data Present",
    title = "Data presence after processing to matrices"
  )

gg2


comp <- gg + gg2

ggsave(filename = paste0(plot_dir, "comparing_beforeafter_processing.png"),
       plot = comp, width = 17, height = 7, dpi = 300)


#-------------------------------------------------------------
#now plotting the raw data with the tag ID data and the amount of coverage per individual

# Step 1: Summarize each tag per individual
tag_summary <- all_nab %>%
  group_by(individual.local.identifier, tag.local.identifier) %>%
  dplyr::summarise(
    start_date = min(as.Date(timestamp)),
    end_date   = max(as.Date(timestamp)),
    .groups = "drop"
  )

# Step 2: Count data points per month per tag
tag_monthly <- all_nab %>%
  mutate(month_year = floor_date(as.Date(timestamp), "month")) %>%
  group_by(individual.local.identifier, tag.local.identifier, month_year) %>%
  dplyr::summarise(n_fixes = n(), .groups = "drop")

all_ids$individual.local.identifier <- all_ids$name

#adding pride ID to the plot
tag_monthly <- tag_monthly %>%
  left_join(all_ids %>% select(individual.local.identifier, pride_id),
            by = "individual.local.identifier")

tag_summary <- tag_summary %>%
  left_join(all_ids %>% select(individual.local.identifier, pride_id),
            by = "individual.local.identifier")


# Step 3: Plot tag timelines with monthly counts as color intensity
g4 <- ggplot() +
  # Timeline bars for tag deployment
  geom_segment(data = tag_summary,
               aes(x = start_date, xend = end_date,
                   y = individual.local.identifier,
                   yend = individual.local.identifier,
                   color = as.factor(tag.local.identifier)),
               size = 17, alpha = 0.7,
               show.legend = FALSE) +
  # Tiles for monthly fix counts
  geom_tile(data = tag_monthly,
            aes(x = month_year, y = individual.local.identifier,
                fill = n_fixes, width = 25), height = 0.5) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  )+
  theme_classic(base_size = 25) +
  labs(
    x = "Date",
    y = "Individual",
    fill = "Fixes per Month",
    color = "Tag ID",
    title = "Tag Deployment and Monthly Data Coverage per Individual"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))#+
  #facet_wrap(~pride_id, scales = "free_y")


ggsave(filename = paste0(plot_dir, "rawdata_withtaginfo.png"),
       plot = g4, width = 25, height = 26, dpi = 300)


#add to this plot lines for each individual the date they died

dod <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/metadata/lion_dod.csv")

dod_ym <- dod %>%
  mutate(
    DOD_clean = trimws(DOD_clean),
    DOD_clean = ifelse(DOD_clean == "" | DOD_clean == "still alive", NA, DOD_clean)
  ) %>%
  # Keep only rows with YYYY-MM or YYYY-MM-DD pattern
  filter(grepl("^\\d{4}-\\d{2}", DOD_clean)) %>%
  mutate(
    DOD_clean = ifelse(grepl("^\\d{4}-\\d{2}$", DOD_clean),
                       paste0(DOD_clean, "-15"), DOD_clean), #make exact day to middle of month for plotting
    DOD_clean = as.Date(DOD_clean, format = "%Y-%m-%d")
  )
tag_summary <- tag_summary %>%
  left_join(dod_ym %>% select(id, DOD_clean),
            by = c("individual.local.identifier" = "id"))

tag_summary$DOD_clean <- as.Date(tag_summary$DOD_clean)


#add the inds we know died but not sure when:
dod_unk <- dod %>%
  mutate(
    DOD_clean = trimws(DOD_clean),
    DOD_clean = ifelse(DOD_clean == "" | DOD_clean == "still alive", NA, DOD_clean)
  ) %>%
  # Keep only rows with just the year (YYYY)
  filter(grepl("^\\d{4}$", DOD_clean)) %>%
  mutate(
    # Make a middle-of-year date for plotting
    DOD_clean = as.Date(paste0(DOD_clean, "-07-01"))
  )

ts_max <- as.Date("2026-01-01")
dod_unk <- dod_unk %>%
  mutate(
    plot_date = ts_max   # all points at the end of the timeline
  )

dod_unk <- dod_unk %>%
  mutate(year_label = substr(DOD_clean, 1, 4)) 


g5 <- g4 +
  geom_point(data = tag_summary %>% filter(!is.na(DOD_clean)),
                    aes(x = DOD_clean, y = individual.local.identifier),
                    shape = 21, size = 10, fill = "red", color = "black") + 
  geom_point(data = dod_unk %>% left_join(tag_summary, by = c("id" = "individual.local.identifier")),
           aes(x = plot_date, y = id),
           shape = 21, size = 5, fill = "blue", color = "black") +
  scale_fill_gradient(low = "white", high = "black")+
  geom_text(data = dod_unk %>% left_join(tag_summary, by = c("id" = "individual.local.identifier")),
            aes(x = plot_date + 15,  # shift slightly to the right
                y = id,
                label = year_label),
            hjust = 0, size = 10, color = "blue") 

ggsave(filename = paste0(plot_dir, "rawdata_tagid_dod.png"),
       plot = g5, width = 25, height = 26, dpi = 300)



#looking at data coverage per pride for the processed data

# Join pride_id to coverage_df
coverage_df <- coverage_df %>%
  left_join(all_ids %>% select(name, pride_id),
            by = c("individual" = "name"))

# Plot with faceting by pride
g5 <- ggplot(coverage_df, aes(x = month_year, y = individual, fill = n_fixes)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "purple") +
  theme_classic() +
  labs(
    x = "Month-Year",
    y = "Individual",
    fill = "Fixes",
    title = "Processed GPS Data Coverage"
  ) +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~pride_id, scales = "free_y")


ggsave(filename = paste0(plot_dir, "processeddata_bypride.png"),
       plot = g5, width = 20, height = 12, dpi = 300)



#add dod and unk dod to all_ids dataframe
dod_unk$DOD_clean <- as.Date(dod_unk$DOD_clean)
dod_unk <- dod_unk[,1:9]
dod_unk$known <- F
dod_ym$DOD_clean <- as.Date(dod_ym$DOD_clean)
dod_ym$known <- T
all_dod <- rbind(dod_ym, dod_unk)


all_ids2 <- merge(all_ids, all_dod[,c(1,9,10)], by.x = "name", by.y = "id", all = TRUE)

# Reorder files to match the order in all_ids$name
all_ids <- all_ids2[match(allNames, all_ids2$name), ]

rownames(all_ids) <- 1:nrow(all_ids)

save(all_ids, file = 'allnamib_ids_withdod.RData')


