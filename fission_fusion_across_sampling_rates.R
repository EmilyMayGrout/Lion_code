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

#----------FUNCTIONS------------------------------------------------------------------------

mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}


subsample_matrix <- function(matrix, ts, interval) {
  # Convert the interval to seconds
  interval_seconds <- as.numeric(difftime(as.POSIXct(interval, format="%H:%M:%S"), as.POSIXct("00:00:00", format="%H:%M:%S"), units = "secs"))
  
  # Find the indices of the time points that match the interval
  subsample_indices <- seq(1, length(ts), by = interval_seconds / (15 * 60))
  
  # Subsample the matrix and the time vector
  subsampled_matrix <- matrix[, subsample_indices, drop = FALSE]
  subsampled_ts <- ts[subsample_indices]
  
  return(list(subsampled_matrix = subsampled_matrix, subsampled_ts = subsampled_ts, subsample_indices = subsample_indices))
}


# Function to subsample the matrix at midday
subsample_matrix_midday <- function(matrix, ts) {
  # Convert the time vector to Date and Time components
  dates <- as.Date(ts)
  times <- format(ts, "%H:%M:%S")
  
  # Find the indices of the time points that match midday (12:00:00)
  midday_indices <- which(times == "12:00:00")
  
  # Subsample the matrix and the time vector
  subsampled_matrix <- matrix[, midday_indices, drop = FALSE]
  subsampled_ts <- ts[midday_indices]
  
  return(list(
    subsampled_matrix = subsampled_matrix,
    subsampled_ts = subsampled_ts,
    original_ts = ts,
    subsample_indices = midday_indices
  ))
}

# Function to subsample the matrix every 14 days at midday
subsample_matrix_14days_midday <- function(matrix, ts) {
  # Convert the time vector to Date and Time components
  dates <- as.Date(ts)
  times <- format(ts, "%H:%M:%S")
  
  # Find the indices of the time points that match midday (12:00:00)
  midday_indices <- which(times == "12:00:00")
  
  # Filter the dates to get every 7th day
  unique_dates <- unique(dates[midday_indices])
  selected_dates <- unique_dates[seq(1, length(unique_dates), by = 14)]
  
  # Find the indices that match the selected dates
  selected_indices <- which(dates %in% selected_dates & times == "12:00:00")
  
  # Subsample the matrix and the time vector
  subsampled_matrix <- matrix[, selected_indices, drop = FALSE]
  subsampled_ts <- ts[selected_indices]
  
  return(list(
    subsampled_matrix = subsampled_matrix,
    subsampled_ts = subsampled_ts,
    original_ts = ts,
    subsample_indices = selected_indices
  ))
}

# Function to subsample the matrix every N days at midday
subsample_matrix_ndays_midday <- function(matrix, ts, n_days) {
  # Convert the time vector to Date and Time components
  dates <- as.Date(ts)
  times <- format(ts, "%H:%M:%S")
  
  # Find the indices of the time points that match midday (12:00:00)
  midday_indices <- which(times == "12:00:00")
  
  # Filter the dates to get every Nth day
  unique_dates <- unique(dates[midday_indices])
  selected_dates <- unique_dates[seq(1, length(unique_dates), by = n_days)]
  
  # Find the indices that match the selected dates
  selected_indices <- which(dates %in% selected_dates & times == "12:00:00")
  
  # Subsample the matrix and the time vector
  subsampled_matrix <- matrix[, selected_indices, drop = FALSE]
  subsampled_ts <- ts[selected_indices]
  
  return(list(
    subsampled_matrix = subsampled_matrix,
    subsampled_ts = subsampled_ts,
    original_ts = ts,
    subsample_indices = selected_indices
  ))
}


# Function to extract hours from the file name
extract_hours <- function(file_name) {
  # Extract the time part using regular expressions
  time_part <- sub("subsampled_matrix_(\\d+_\\d+_\\d+).*", "\\1", file_name)
  
  # Check if the time part is in the format "HH_MM_SS"
  if (grepl("^\\d+_\\d+_\\d+$", time_part)) {
    # Split the time part into hours, minutes, and seconds
    time_parts <- as.numeric(unlist(strsplit(time_part, "_")))
    hours <- time_parts[1] + time_parts[2] / 60 + time_parts[3] / 3600
  } else if (grepl("(\\d+)days_midday", file_name)) {
    # Extract the number of days and convert to hours
    days <- as.numeric(sub(".*_(\\d+)days_midday", "\\1", file_name))
    hours <- days * 24
  } else {
    hours <- NA
  }
  
  return(hours)
}



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
gps_file <- "lion_xy_15min_level0.RData" #level0 is when Venus is not removed
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

#read in library of functions
setwd(coati_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)


#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}


#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


R = 100
subgroup_data <- get_subgroup_data(xs, ys, R)

ind_subgroup_membership <- subgroup_data$ind_subgroup_membership

#---------------------------------------------------------------------------------------------------

#GET SAMPLING REGIMES:

# Create a list to store the results
subsampled_results <- list()

# Assuming ind_subgroup_membership is your matrix and ts is your time vector
intervals <- c("00:15:00", "00:30:00", "01:00:00", "02:00:00", "04:00:00", "12:00:00")

# Subsample the matrix at different intervals
for (interval in intervals) {
  result <- subsample_matrix(ind_subgroup_membership, ts, interval)
  matrix_name <- paste0("subsampled_matrix_", gsub(":", "_", interval))
  subsampled_results[[matrix_name]] <- list(
    matrix = result$subsampled_matrix,
    subsampled_ts = result$subsampled_ts,
    original_ts = result$original_ts,
    subsample_indices = result$subsample_indices
  )
}
# Assuming ind_subgroup_membership is your matrix and ts is your time vector
day_intervals <- c(1, 2, 3, 6, 9, 12, 15, 20)

# Subsample the matrix every N days at midday
for (n_days in day_intervals) {
  result_ndays_midday <- subsample_matrix_ndays_midday(ind_subgroup_membership, ts, n_days)
  matrix_name <- paste0("subsampled_matrix_", n_days, "days_midday")
  subsampled_results[[matrix_name]] <- list(
    matrix = result_ndays_midday$subsampled_matrix,
    subsampled_ts = result_ndays_midday$subsampled_ts,
    original_ts = result_ndays_midday$original_ts,
    subsample_indices = result_ndays_midday$subsample_indices
  )
}



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
      sub_ids_i <- samplerate$matrix[i,]
      sub_ids_j <- samplerate$matrix[j,]
      
      # Computing edge weight (fraction of time in same subgroup)
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
    }
  }
  
  diag(ff_net) <- NA
  new_order <- c(1:3, 5, 4, 6)
  ffnet_reorder <- ff_net[new_order, new_order]
  
  # Save the plot with the correct title
  png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir, name, '_subgroup_network_100m.png'))
  
  visualize_network_matrix_trago(ffnet_reorder, lion_ids[new_order,])
  
  dev.off()  # Close the PNG device
}




#get splits and merges (code copied from lion_split_merge_dfs.R)

#-----SPLITS---------

# Initialize split_list_samplingrates to store all splits data frames
split_list_samplingrates <- list()

# Get the names of the elements in the list
element_names <- names(subsampled_results)

for (name in element_names) {
  
  samplerate <- subsampled_results[[name]]
  
  # Initialize splits vector
  splits <- c()
  
  n_times_cut <- ncol(samplerate$matrix)
  
  # Run through time by time and check if each time is a split
  for (t in 1:(n_times_cut - 1)) {
    
    # Get subgroup membership now and later
    subgroups_now <- samplerate$matrix[, t]
    subgroups_later <- samplerate$matrix[, t + 1]
   
    # If we have a time of completely NAs, pass
    if (sum(!is.na(subgroups_now)) == 0 | sum(!is.na(subgroups_later)) == 0) {
      next
    }
    
    # Transfer NAs from now to later and later to now
    subgroups_now[which(is.na(subgroups_later))] <- NA
    subgroups_later[which(is.na(subgroups_now))] <- NA
    
    # Get number of subgroups now and in next time step (later)
    n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
    n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
    
    # Determine if this time step is a split - include cases where singletons leave
    if ((n_subgroups_now == 1) &
        n_subgroups_later > n_subgroups_now) {
      splits <- c(splits, t)
    }
  }
  
  # Make a data frame of splits, with original group and subgroups
  splits_df <- data.frame(t = splits, orig_group = NA, sub1 = NA, sub2 = NA, sub3 = NA, sub4 = NA, sub5 = NA, subsample_indices = NA)
  
  for (i in 1:nrow(splits_df)) {
    t <- splits_df$t[i]
    subgroups_now <- samplerate$matrix[, t]
    subgroups_later <- samplerate$matrix[, t + 1]
    
    # Transfer NAs from now to later and later to now so erroneous splits are not detected
    subgroups_now[which(is.na(subgroups_later))] <- NA
    subgroups_later[which(is.na(subgroups_now))] <- NA
    
    # Identify individuals whose subgroup ID has changed
    changed_indices <- which(subgroups_now != subgroups_later)
    group_id_changers <- subgroups_later[changed_indices]
    prev_id_changers <- subgroups_now[changed_indices]
    
    # Put the IDs of the original subgroup into df
    splits_df$orig_group[i] <- list(which(subgroups_now == prev_id_changers[1]))
    
    # Find the groups where the original members went
    group_ids_later <- unique(subgroups_later[changed_indices])
    
    for (j in 1:length(group_ids_later)) {
      group_id <- group_ids_later[j]
      inds_in_group <- which(subgroups_later == group_id)
      orig_inds_in_group <- intersect(inds_in_group, changed_indices) # Only count the original group members who changed
      
      # Put lists into a data frame
      if (j == 1) {
        splits_df$sub1[i] <- list(orig_inds_in_group)
      }
      if (j == 2) {
        splits_df$sub2[i] <- list(orig_inds_in_group)
      }
      if (j == 3) {
        splits_df$sub3[i] <- list(orig_inds_in_group)
      }
      if (j == 4) {
        splits_df$sub4[i] <- list(orig_inds_in_group)
      }
      if (j == 5) {
        splits_df$sub5[i] <- list(orig_inds_in_group)
      }
    }
    # Add subsample indices
    splits_df$subsample_indices[i] <- list(samplerate$subsample_indices[t + 1]) #so splits and merges aren't at the same time 
  }
  
  # Number in each subgroup
  splits_df$n_orig <- sapply(splits_df$orig_group, function(x) { return(sum(!is.na(x))) })
  splits_df$n_sub1 <- sapply(splits_df$sub1, function(x) { return(sum(!is.na(x))) })
  splits_df$n_sub2 <- sapply(splits_df$sub2, function(x) { return(sum(!is.na(x))) })
  splits_df$n_sub3 <- sapply(splits_df$sub3, function(x) { return(sum(!is.na(x))) })
  
  # Append splits_df to split_list_samplingrates
  split_list_samplingrates[[name]] <- splits_df
  
}

#save list
save(split_list_samplingrates, file = paste0(data_dir, "lion_merge_list.Rdata"))



#--------MERGE------------------

# Initialize merge_all to store all splits data frames
merge_list_samplingrates <- list()

# Get the names of the elements in the list
element_names <- names(subsampled_results)

for (name in element_names) {
  
  samplerate <- subsampled_results[[name]]
  
  # Run through time by time and check if each time is a merge
  merge <- c()
  
  n_times_cut <- ncol(samplerate$matrix)
  
  for (t in 2:n_times_cut) {  # Start from 2 to avoid negative index
    
    # Get subgroup membership now and previous
    merge_group <- samplerate$matrix[, t]
    subgroups_previously <- samplerate$matrix[, t - 1]
    
    # If we have a time of completely NAs, pass
    if (sum(!is.na(merge_group)) == 0 | sum(!is.na(subgroups_previously)) == 0) {
      next
    }
    
    # Transfer NAs from now to later and later to now
    merge_group[which(is.na(subgroups_previously))] <- NA
    subgroups_previously[which(is.na(merge_group))] <- NA
    
    # Get number of subgroups now and in previous time step
    n_merge_group <- length(unique(merge_group[!is.na(merge_group)]))
    n_subgroups_previously <- length(unique(subgroups_previously[!is.na(subgroups_previously)]))
    
    # Determine if this time step is a merge
    if (n_subgroups_previously > n_merge_group) {
      merge <- c(merge, t)
    }
  }
  
  # Make a data frame of merges, with merged group and subgroups
  merge_df <- data.frame(t = merge, merge_group = NA, sub1 = NA, sub2 = NA, sub3 = NA, sub4 = NA, sub5 = NA, subsample_indices = NA)
  
  for (i in 1:nrow(merge_df)) {
    t <- merge_df$t[i]
    merge_group <- samplerate$matrix[, t]
    subgroups_previously <- samplerate$matrix[, t - 1]
    
    # Transfer NAs from now to later and later to now
    merge_group[which(is.na(subgroups_previously))] <- NA
    subgroups_previously[which(is.na(merge_group))] <- NA
    
    # Identify individuals whose subgroup ID has changed
    changed_indices <- which(merge_group != subgroups_previously)
    group_id_changers <- merge_group[changed_indices]
    prev_id_changers <- subgroups_previously[changed_indices]
    
    # Put the IDs of the subgroup where the merge occurs into df
    merge_df$merge_group[i] <- list(which(merge_group == group_id_changers[1]))
    
    # Find the groups where the original members were from
    group_ids_before <- unique(subgroups_previously[changed_indices])
    
    for (j in 1:length(group_ids_before)) {
      group_id <- group_ids_before[j]
      inds_in_group <- which(subgroups_previously == group_id)
      orig_inds_in_group <- intersect(inds_in_group, changed_indices) # Only count the original group members who changed
      
      # Put lists into a data frame
      if (j == 1) {
        merge_df$sub1[i] <- list(orig_inds_in_group)
      }
      if (j == 2) {
        merge_df$sub2[i] <- list(orig_inds_in_group)
      }
      if (j == 3) {
        merge_df$sub3[i] <- list(orig_inds_in_group)
      }
      if (j == 4) {
        merge_df$sub4[i] <- list(orig_inds_in_group)
      }
      if (j == 5) {
        merge_df$sub5[i] <- list(orig_inds_in_group)
      }
    }
    # Add subsample indices
    merge_df$subsample_indices[i] <- list(samplerate$subsample_indices[t])
  }
  
  
  # Number in each subgroup
  merge_df$n_merge <- sapply(merge_df$merge_group, function(x) { return(sum(!is.na(x))) })
  merge_df$n_sub1 <- sapply(merge_df$sub1, function(x) { return(sum(!is.na(x))) })
  merge_df$n_sub2 <- sapply(merge_df$sub2, function(x) { return(sum(!is.na(x))) })
  merge_df$n_sub3 <- sapply(merge_df$sub3, function(x) { return(sum(!is.na(x))) })
  
  # Append merge_df to list_all
  merge_list_samplingrates[[name]] <- merge_df
  
}

#save list
save(merge_list_samplingrates, file = paste0(data_dir, "lion_merge_list.Rdata"))


#-------------------------------------------------------------------------

# Initialize list to store time differences for each sampling rate
time_diff_list <- list()

# Get the names of the elements in the list
element_names <- names(split_list_samplingrates)

for (n in element_names) {
  
  splits_df <- split_list_samplingrates[[n]]
  merge_df <- merge_list_samplingrates[[n]]
  
  # Add event type to dataframes
  merge_df$event <- "merge"
  splits_df$event <- "split"
  
  # Subset relevant columns
  merge_df_subset <- merge_df[, c("subsample_indices", "event")]
  splits_df_subset <- splits_df[, c("subsample_indices", "event")]
  
  # Combine the data frames and ensure subsample_indices is numeric
  time_diff <- rbind(splits_df_subset, merge_df_subset) %>%
    mutate(subsample_indices = as.numeric(subsample_indices)) %>%
    arrange(subsample_indices)
  
 
  # Initialize keep columns
  time_diff$keep <- 0
  time_diff$keep2 <- 0
  
  # Mark rows to keep based on event sequence
  # duration is from the last merge to the first split
  # and from the first split to the next merge
  for (i in 1:(nrow(time_diff) - 1)) {
    if (time_diff$event[i] == "merge" & time_diff$event[i + 1] == "split") {
      time_diff$keep[i] <- 1
      time_diff$keep2[i + 1] <- 1
    }
  }
  
  # Filter rows to keep
  time_diff$keep3 <- time_diff$keep + time_diff$keep2
  time_diff <- time_diff[time_diff$keep3 == 1, c("subsample_indices", "event")]
  
  # Separate splits and merges
  splits_time_diff <- time_diff[time_diff$event == "split", ]
  merge_time_diff <- time_diff[time_diff$event == "merge", ]
  
  time_diff <- data.frame(cbind(splits_time_diff, merge_time_diff))
  time_diff <- time_diff[,-c(2,4)]
  colnames(time_diff) <- c("splits_t","merge_t")

  # Calculate time differences
  time_diff$diff <- as.numeric(time_diff$splits_t) - as.numeric(time_diff$merge_t)
  time_diff$diff_time_hour <- (time_diff$diff * 10) / 60
  time_diff$diff_time_hour <- format(round(time_diff$diff_time_hour, 1), nsmall = 1)
  time_diff$diff_time_hour <- as.numeric(time_diff$diff_time_hour)
  
  # Store the result in the list - time differences between the last merge and the first split
  time_diff_list[[n]] <- time_diff
}




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

#plot number of splits with changes in sampling effort

# Create a data frame to store the filenames and number of rows
length_sampling_df <- data.frame(
  filename = character(),
  nrow = integer(),
  stringsAsFactors = FALSE
)

# Populate the data frame
for (name in names(split_list_samplingrates)) {
  subsampled_matrix <- merge_list_samplingrates[[name]]
  nrow_matrix <- nrow(subsampled_matrix)
  
  length_sampling_df <- rbind(length_sampling_df, data.frame(
    filename = name,
    nrow = nrow_matrix,
    stringsAsFactors = FALSE
  ))
}

#adding sampling rate to dataframe
length_sampling_df$sampling_rate <- c(15, 30, 60, 120, 240, 720, 1440,2880,4320,8640,12960,17280, 21600, 28800)


# Plot the data points and the exponential regression line
ggplot(length_sampling_df, aes(sampling_rate, nrow)) +
  geom_point() +  # Scatter plot
  geom_line() +
  theme_minimal()


#look at the number of individuals that split/merge - are larger subgroups more likely to split, and singles more likely to merge?

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


