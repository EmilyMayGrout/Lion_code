#this script is changing the sampling rates of the GPS frequencies and quantifying how much data is lost when we start to subsample and change the sampling efforts

#how often, when, and for how long prides split, and conversely when and where they come together; how far apart subgroups are when they split; and how stable group affiliation is. 
#We predict that prides will exhibit lower spatial cohesion during dry season when resources are less abundant and more widely dispersed and that larger prides will exhibit greater FF dynamics with increased frequency of fission events and greater proportion of time apart and that group and subgroup size will scale as a function of overall resource richness.
#Q. How high resolution do you need to capture consistent group level patterns
#We could simply subsample out from a fix rate of one per 15 minutes down to one per four hours. Need to decide on increments maybe: 

#15 mins
#30 mins
#60 mins
#2hrs
#4hrs
#do subsampling schedule for 2h/4h 6am-6pm every 4h, 6pm-6am every 2h 

#----------------------------------------------------------------------------------------------

#Do we want to go down to even one fix per day? Either at midday when typically resting, or midnight when typically more likely to be active
#It would be helpful (specifically for us and the ORP data) to test the Wide Horizons schedule of fixes at: 08h, 12h, 16h, 20h, 22h, 00h, 02h, 04h, 06h
#Sight collared animal(s) once every two weeks [common in study sites where animals are hard to locate without recent GPS location information or VHF]
#This could either be approached as, what are the results if you just do this, or how often would you need to sight them (and over what time period) to see all individuals
# for the 2 week sighting, could do a smaller radius as unlikely to spot all inds if they're far away?

data_dir <- "C:/Users/egrout/Dropbox/lion/data/processed/"
code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
coati_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/lion/results/level0/'
gps_file <- "lion_xy_15min_level0_2yr.RData" #level0 is when Venus is not removed
id_file <- 'lion_ids.RData'

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(ggplot2)
library(vioplot)
library(plotly)
library(patchwork)
library(ggpmisc)
library(cocomo)
library(reshape2)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')
source('lion_functions.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#removing data when sample regime went from 15 mins to 2 hours
non_na_counts <- colSums(!is.na(xs))
drops <- which(diff(non_na_counts) < -4)
xs <- xs[,1:drops[1]]
ys <- ys[,1:drops[1]]
ts <- ts[1:drops[1]]

R = 60
subgroup_data <- get_subgroup_data(xs, ys, R)

ind_subgroup_membership <- subgroup_data$ind_subgroup_membership

#----------------------------------------------------------------

#------RANDOMISE WHEN START SAMPLES ARE TAKEN--------------------



#GET SAMPLING REGIMES:

intervals <- c("00:15:00", "00:30:00", "01:00:00", "02:00:00", "04:00:00", "12:00:00")

subsampled_time_results <- subsample_multiple_intervals(ind_subgroup_membership, ts, intervals) 

day_intervals <- c(1, 2, 3, 6, 9, 12, 15, 20, 30, 60, 90)

subsampled_day_results <- subsample_multiple_day_intervals(ind_subgroup_membership, ts, day_intervals)

subsampled_results <- c(subsampled_time_results, subsampled_day_results)


run_this <- F

if(run_this == T){

#get list of subsample results by randomising the start time of extracting

# Prepare a data frame to collect split event counts
n_offsets <- 96  # 15-minute intervals in 24 hours

intervals <- c("00:15:00", "00:30:00", "01:00:00", "02:00:00", "04:00:00", "12:00:00")


# Initialize data frames to store split and merge counts
split_events_all <- data.frame(matrix(nrow = n_offsets, ncol = length(intervals)))
merge_events_all <- data.frame(matrix(nrow = n_offsets, ncol = length(intervals)))
colnames(split_events_all) <- intervals
colnames(merge_events_all) <- intervals

for (n in 1:n_offsets) {
  ind_subgroup_membership_cut <- ind_subgroup_membership[, n:ncol(ind_subgroup_membership)]
  ts_cut <- ts[n:length(ts)]
  
  subsampled_time_results_cut <- subsample_multiple_intervals(ind_subgroup_membership_cut, ts_cut, intervals)
  
  
  # Detect events for each interval
  for (interval in seq_along(subsampled_time_results_cut)){
    
  ind_sub_memb <- subsampled_time_results_cut[[interval]]$ind_subgroup_membership
    
  split_count <- nrow(detect_splits(ind_sub_memb))
  merge_count <- nrow(detect_merges(ind_sub_memb))
  
  split_events_all[n, interval] <- split_count
  merge_events_all[n, interval] <- merge_count
  
  }
  
  print(paste("Offset", n, "done"))
}


# Add identifiers
split_events_all$offset <- 1:n_offsets
merge_events_all$offset <- 1:n_offsets

# Reshape to long format
split_long <- melt(split_events_all, id.vars = "offset", variable.name = "interval", value.name = "n_events")
merge_long <- melt(merge_events_all, id.vars = "offset", variable.name = "interval", value.name = "n_events")

split_long$event_type <- "split"
merge_long$event_type <- "merge"

split_merge_time_offset <- rbind(split_long, merge_long)

}



#save(split_merge_time_offset,split_events_all, merge_events_all, file = paste0(data_dir, "time_offset_split_merge_60m.Rdata"))

load(file = paste0(data_dir, "time_offset_split_merge_60m.Rdata"))

# Plot
g0 <- ggplot(split_merge_time_offset, aes(x = offset, y = n_events, color = interval)) +
  geom_line() +
  facet_wrap(~event_type, scales = "free_y") +
  labs(title = "Variation in Detected Splits and Merges by Sampling Offset",
       x = "15-min Start Offset",
       y = "Number of Events")+
  theme_classic()

#ggsave(filename = paste0(plot_dir, 'samplingrate_offset_2yr.png'), plot = g0, width = 10, height = 5, dpi = 300)


#randomise start day for multiple day extraction - change start date across one month

run_this <- F
if(run_this == T){

day_shifts <- 0:29  # Simulate 30 different days

day_intervals <- c(1: 7, 14, 20, 30, 60, 90)


# Initialize storage for daily shift results
split_events_by_day <- data.frame(matrix(nrow = length(day_shifts), ncol = length(day_intervals)))
merge_events_by_day <- data.frame(matrix(nrow = length(day_shifts), ncol = length(day_intervals)))
colnames(split_events_by_day) <- day_intervals
colnames(merge_events_by_day) <- day_intervals

# Find all midday indices in ts
midday_indices <- which(format(ts, "%H:%M:%S") == "10:00:00")
midnight_indices <- which(format(ts, "%H:%M:%S") == "02:00:00")

for (i in seq_along(day_shifts)) {
  start_index <- midnight_indices[day_shifts[i] + 1]
  if (is.na(start_index)) next
  
  ind_subgroup_membership_day <- ind_subgroup_membership[, start_index:ncol(ind_subgroup_membership)]
  ts_day <- ts[start_index:length(ts)]
  
  subsampled_day_results <- subsample_multiple_day_intervals(ind_subgroup_membership_day, ts_day, day_intervals)
  
  #go through each subsample for each day shift
  for (j in seq_along(subsampled_day_results)) {
    subsampled_i_j <- subsampled_day_results[[j]]
    
    split_list_day <- detect_splits(subsampled_i_j$ind_subgroup_membership)
    merge_list_day <- detect_merges(subsampled_i_j$ind_subgroup_membership)
    
    split_events_by_day[i, j] <- nrow(split_list_day)
    merge_events_by_day[i, j] <- nrow(merge_list_day)
    
  }
  
  print(paste("Day shift", i, "done"))
}



# Add a column for day offset
split_events_by_day$day_shift <- 0:(nrow(split_events_by_day) - 1)
merge_events_by_day$day_shift <- 0:(nrow(merge_events_by_day) - 1)

# Convert to long format
split_day_long <- melt(split_events_by_day, id.vars = "day_shift", 
                       variable.name = "interval", value.name = "n_events")
merge_day_long <- melt(merge_events_by_day, id.vars = "day_shift", 
                       variable.name = "interval", value.name = "n_events")

# Label event type
split_day_long$event_type <- "split"
merge_day_long$event_type <- "merge"

# Combine both
split_merge_day_offset <- rbind(split_day_long, merge_day_long)

}


#save(split_merge_day_offset, split_events_by_day, merge_events_by_day, file = paste0(data_dir, "day_offset_split_merge_60m.Rdata"))

Singletons <- T

if(Singletons == T){

load(file = paste0(data_dir, "time_offset_split_merge_60m.Rdata"))
load(file = paste0(data_dir, "day_offset_split_merge_60m.Rdata"))
title <- "Variation in Detected Splits and Merges by Start Day"
title2 <- "Distribution of Split and Merge Events by Sampling Interval"
filename <- 'sampling_interval_splitsmerges_withsingles_2yr.png'

}else{
 
load(file = paste0(data_dir, "time_offset_split_merge_no_singletons.Rdata"))
load(file = paste0(data_dir, "day_offset_split_merge_no_singletons.Rdata"))
title <- "Variation in Detected Splits and Merges by Start Day without Singletons"
title2 <- "Distribution of Split and Merge Events by Sampling Interval - Without Singleton Events"
filename <- 'sampling_interval_splitsmerges_withoutsingles_2yr.png'
}


gg0 <- ggplot(split_merge_day_offset, aes(x = day_shift, y = n_events, color = interval)) +
  geom_line() +
  facet_wrap(~event_type, scales = "free_y") +
  labs(
    title = title,
    x = "Day Offset",
    y = "Number of Events"
  ) +
  theme_classic()

#ggsave(filename = paste0(plot_dir, 'samplingrate_day_offset_2yr.png'), plot = gg0, width = 10, height = 5, dpi = 300)

split_merge_day_offset$offset_type <- "day"
split_merge_time_offset$offset_type <- "time"

# Convert interval to character if it's still a factor, for better control over x-axis
split_merge_day_offset$interval <- as.numeric(levels(split_merge_day_offset$interval))[split_merge_day_offset$interval]
split_merge_day_offset$event_type <- as.factor(split_merge_day_offset$event_type)

split_merge_time_offset$interval <- as_hms(as.character(split_merge_time_offset$interval))
split_merge_time_offset$event_type <- as.factor(split_merge_time_offset$event_type)


# Create the boxplot
f1 <- ggplot(split_merge_day_offset, aes(x = interval, y = n_events, fill = event_type, group = interaction(interval, event_type))) +
  geom_boxplot(outlier.shape = NA, size = 0.2, lwd = 0.1, width = 2, alpha = 0.7, position = position_dodge(width = 3)) +
  scale_fill_manual(values=c("gold", "slateblue2")) +
  labs(
   x = "Sampling Interval (days)",
    y = "Number of Events",
    fill = "Event Type"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")


f2 <- ggplot(split_merge_time_offset, aes(x = interval, y = n_events, fill = event_type, group = interaction(interval, event_type))) +
  geom_boxplot(outlier.shape = NA, width = 1500, lwd = 0.1, alpha = 0.7, position = position_dodge(width = 1500)) + 
  scale_fill_manual(values=c("gold", "slateblue2")) +
  labs(
    title = title2,
    x = "Sampling Interval (hour)",
    y = "Number of Events",
    fill = "Event Type"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

f <- f2+f1

f

ggsave(filename = paste0(plot_dir, filename), plot = f, width = 12, height = 5, dpi = 300)

#---------------------------------------------------------------------------------------------------

#metrics to pull out at the different sampling regimes
#how often, when, and for how long prides split, and conversely when and where they come together; how far apart subgroups are when they split; and how stable group affiliation is.
#how many splits occurred
#when do splits and merges occur - not sure how to do this when there are big gaps - I guess there's just more error?
#how long are the group split
#how far do subgroups split
#how stable are subgrouping associations

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)

# Get the names of the elements in the list
element_names <- names(subsampled_results)

# Iterate over the element names and corresponding elements
for (name in element_names) {
  samplerate <- subsampled_results[[name]]
  
  for (i in 1:n_inds) {
    for (j in 1:n_inds) {
      
      # Getting subgroup id for individual i and j
      sub_ids_i <- samplerate$ind_subgroup_membership[i,]
      sub_ids_j <- samplerate$ind_subgroup_membership[j,]
      
      # Computing edge weight (fraction of time in same subgroup)
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
    }
  }
  
  diag(ff_net) <- NA
  new_order <- c(1:3, 5, 4, 6)
  ffnet_reorder <- ff_net[new_order, new_order]
  
  # Save the plot with the correct title
  png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir, name, '_subgroup_network_60m_2yr.png'))
  
  visualize_network_matrix_trago(ffnet_reorder, lion_ids[new_order,])
  mtext(name)
  
  dev.off()  # Close the PNG device
}

#-----SPLITS---------


split_list_samplingrates <- lapply(subsampled_results, function(x) {
  detect_splits(x$ind_subgroup_membership)
})

#names(split_list_samplingrates) <- names(subsampled_results)

#save list
save(split_list_samplingrates, file = paste0(data_dir, "lion_split_list.Rdata"))



#--------MERGE------------------

merge_list_samplingrates <- lapply(subsampled_results, function(x) {
  detect_merges(x$ind_subgroup_membership)
})

#names(merge_list_samplingrates) <- names(subsampled_results)

#save list
save(merge_list_samplingrates, file = paste0(data_dir, "lion_merge_list.Rdata"))



#-------------------------------------------------------------------------

#Get time difference from the last merge to the first split and from the first split to the next merge

time_diff_list <- compute_merge_split_time_differences(split_list_samplingrates, merge_list_samplingrates)

# Get the names of the elements in the list
element_names <- names(time_diff_list)

for (name in element_names) {
  time_diff <- time_diff_list[[name]]
  
  # Create the plot file name using the list element name
  plot_file_name <- paste0(plot_dir, "event_durations/", name, "_split_duration_hist.png")
  
  # Create the plot
  png(file = plot_file_name, width = 1000, height = 600, units = "px")
  par(mar = c(8, 8, 8, 8))
  hist(time_diff$diff_time_hour, breaks = 80, xlab = "Time between splits and merges (hours)", col = "lightblue3", main = "", cex.lab = 2, cex.axis = 1.5)
  dev.off()
}


#-----------------------------------------------------------------------

#plot number of splits and merges accross sampling rates

length_df <- data.frame(
  filename = names(split_list_samplingrates),
  Splits = sapply(split_list_samplingrates, nrow),
  Merges = sapply(merge_list_samplingrates, nrow),
  stringsAsFactors = FALSE
)
rownames(length_df) <- NULL  # Remove row names

#12 day period has 1 split observation, whereas there are 2 with day 15 and 30 intervals. This could be bias to the days chosen to pull from the data, so would be good to rerun these analyses with randomised start dates (to build a probability distribution) to see whether this affects the distribution of events captured (and will be more accurate representation)


#adding sampling rate to dataframe
length_df$sampling_rate <-  sapply(element_names, get_steps_skipped)

# Convert the sampling rate to hours (assuming 15 minutes = 0.25 hours)
length_df$sampling_days <- length_df$sampling_rate * 0.25 / 24
length_df$sampling_hours <- length_df$sampling_rate * 0.25


length_df_long <- pivot_longer(length_df, cols = c(Splits, Merges),
                               names_to = "Event", values_to = "Count")

# Plot
g1 <- ggplot(length_df_long, aes(x = sampling_days, y = Count, color = Event)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, max(length_df_long$sampling_days), by = 7)) +
  scale_color_manual(values = c("Splits" = "slateblue1", "Merges" = "orange1")) + # Custom colors
  labs(x = "Sampling Rate (days)", y = "Count",
       title = "Number of split and merges extracted across varying sampling rates") +
  theme_classic() +
  theme(legend.position = "none")

#zoom to less than 1 day
g2 <- ggplot(length_df_long, aes(x = sampling_hours, y = Count, color = Event)) +
  geom_line(size = 1) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, max(length_df_long$sampling_hours), by = 10)) +
  xlim(0, 24)+
  scale_color_manual(values = c("Splits" = "slateblue1", "Merges" = "orange1")) + # Custom colors
  labs(x = "Sampling Rate (hours)", y = "Count") +
  theme_classic() 

gg <- g1 + g2
gg
ggsave(filename = paste0(plot_dir, 'splitmerge_diffsamplingrates_withsingles_2yr.png'), plot = gg, width = 10, height = 5, dpi = 300)


#look at distribution of subgroup sizes observed when sampling rate changes

#extract the subgroup sizes from all data and then measure how the distribution of this changes for each subsampling regime






#look at the number of individuals that split/merge - are larger subgroups more likely to split, and singles more likely to merge?

#can't run this as split and merge events ignored single individuals... but can change it to look at this
files <- names(merge_list_samplingrates)

merge_output <- data.frame(files)
merge_output$merge_mean_n_merge <- NA
merge_output$merge_mean_n_sub1 <- NA
merge_output$merge_mean_n_sub2 <- NA
merge_output$merge_n_events <- NA

for(name in files){
  
  merge_i <- merge_list_samplingrates[[name]]
  n_events <- nrow(merge_i)
 
  # Calculate the mean of the n_merge column
  mean_n_merge <- mean(merge_i$n_merge, na.rm = TRUE)
  mean_n_sub1 <- mean(merge_i$n_sub1, na.rm = TRUE)
  mean_n_sub2 <- mean(merge_i$n_sub2, na.rm = TRUE)
  
  
  # Populate the mean_n_merge column in merge_output
  merge_output$merge_mean_n_merge[merge_output$files == name] <- mean_n_merge
  merge_output$merge_mean_n_sub1[merge_output$files == name] <- mean_n_sub1
  merge_output$merge_mean_n_sub2[merge_output$files == name] <- mean_n_sub2
  merge_output$merge_n_events[merge_output$files == name] <- n_events
} 

files <- names(split_list_samplingrates)

split_output <- data.frame(files)
split_output$split_mean_n_orig <- NA
split_output$split_mean_n_sub1 <- NA
split_output$split_mean_n_sub2 <- NA
split_output$split_n_events <- NA


for(name in files){
  
  split_i <- split_list_samplingrates[[name]]
  n_events <- nrow(split_i)
  
  # Calculate the mean of the n_split column
  mean_n_orig <- mean(split_i$n_orig, na.rm = TRUE)
  mean_n_sub1 <- mean(split_i$n_sub1, na.rm = TRUE)
  mean_n_sub2 <- mean(split_i$n_sub2, na.rm = TRUE)
  
  
  # Populate the mean_n_split column in split_output
  split_output$split_mean_n_orig[split_output$files == name] <- mean_n_orig
  split_output$split_mean_n_sub1[split_output$files == name] <- mean_n_sub1
  split_output$split_mean_n_sub2[split_output$files == name] <- mean_n_sub2
  split_output$split_n_events[split_output$files == name] <- n_events
} 


both <- merge(x = merge_output, y = split_output, by = "files", all = TRUE)
#there's more merges than splits generally, so perhaps inds join at different times, but inds leave together 
#so before a split the group is bigger and 1 or 2 inds leave
#before a merge, the group is slightly smaller and 1 or 2 inds join

#extract the number of hours
# Apply the function to the files column
both$hours <- sapply(both$files, extract_hours)

# Plot the data points and the exponential regression line
ggplot(both, aes(hours, split_n_events)) +
  geom_point() +  # Scatter plot
  geom_line() +
  theme_minimal()

ggplot(both, aes(hours, merge_n_events)) +
  geom_point() +  # Scatter plot
  geom_line() +
  theme_minimal()




#-------------------------------------------------------------------

#subsampled 2h/4h == 6am-6pm every 4h, 6pm-6am every 2h

#-------------------------------------------------------------------

subsampled_time_results_2h4h <- subsampled_time_results$subsampled_matrix_02_00_00
time_vec <- as_hms(subsampled_time_results_2h4h$subsampled_ts)

# Define time bounds
start_time <- as_hms("06:00:00") #UTC time which is 8am-8pm local
end_time <- as_hms("18:00:00")  

# Get indices between 6:00 AM and 6:00 PM UTC
day_indices <- which(time_vec >= start_time & time_vec < end_time)
night_indices <- which(!(time_vec >= start_time & time_vec < end_time))

#get indices every 4h for day indices by removing every other time index during the day
day_indices_4h <- day_indices[seq(day_indices[2], length(day_indices), by = 2)]
indices_sub <- sort(c(day_indices_4h, night_indices))
#check its subsampled correctly
t_sub <- subsampled_time_results_2h4h$subsampled_ts[indices_sub]

subsampled_2h4h <- subsampled_time_results_2h4h$ind_subgroup_membership[,indices_sub]

save(subsampled_2h4h, t_sub, indices_sub, file = "subsampled_2h4h_matrix.RData")

