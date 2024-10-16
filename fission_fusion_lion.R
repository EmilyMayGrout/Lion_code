#Questions interested for the lion dataset:
#1. duration of splits
#2. who do they split with
#3. what are they doing when split?
#4. how far do they travel apart?
#5. how often do they come together?

#NATALIA'S TEST CHANGE 
#--------PARAMS-------
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

#----------------
#looking at the sub grouping patterns

#list of Rs
Rs <- c(10,20,30,40,50,100)
Rs <- c(100,200,300,400,500,1000)
Rs <- c(1000,2000,3000,4000,5000,10000)
Rs <- c(10000,20000,30000,40000,50000,100000)

#Number of sub groups when the radius is changed (graph put in dropbox results folder)

png(height = 2200, width = 900, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists_10000-100000m.png'))
par(mfrow=c(6,1), mar = c(8,9,12,1), mgp=c(6,1,0)) #(bottom, left, top, right)

for (i in 1:length(Rs)){
  
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "gold4", breaks = seq(.5,11,1), cex.lab = 4, cex.main = 4, cex.axis=4, freq = FALSE, ylim=c(0,.7), xlim = c(0.5, 11), yaxt = "n", axes = F)
  axis(2, at = c(0,0.3,0.6), cex.axis = 4, las = 1, pos = 0.5)
  axis(1, cex.axis = 4, las = 1, pos = 0, padj = 0.6)
}

dev.off()

#-----------------------------------

#looking at which individuals split with one another

R = 100
subgroup_data <- get_subgroup_data(xs, ys, R)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
new_order <- c(1:3,5,4,6)
ffnet_reorder <- ff_net[new_order, new_order]

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_100m.png'))

visualize_network_matrix_trago(ffnet_reorder, lion_ids[new_order,])
dev.off()

#it looks like Kiara and Sarafina split from the other group members and travelled atleast 50km from them

#how much missing data is there?

png(height = 800, width = 1400, units = 'px', filename = paste0(plot_dir, "number_tracked.png"))
par(mfrow=c(1,2), mar = c(10,9,2,1)) #c(bottom, left, top, right)
sum_tracked <- colSums(!is.na(xs))
hist(sum_tracked[ !(sum_tracked==0) ], main = "", xlab = "Number of individuals tracked", ylab = "Frequency", col = "gold2", cex.lab = 2, cex.axis = 2)

each_sum <- data.frame(sum = rowSums(!is.na(xs)))
barplot(each_sum$sum, names.arg = lion_ids$name, las=2, col = "gold3", ylab = "Number of GPS points",  cex.lab = 2, cex.axis = 2, cex.names=2, mgp=c(5,1,0))

dev.off()

#calculate the proportion of missing data
max(each_sum$sum)
min(each_sum$sum)
each_sum$missing <- max(each_sum$sum)- each_sum$sum
each_sum$prop <- (each_sum$missing/max(each_sum$sum))*100
mean(each_sum$prop)
sd(each_sum$prop)

#------------------------------------------

#what is the absolute diadic distance of the group when the group is together
#saying "together" is when they're 100m from one another
R = 10
subgroup_data <- get_subgroup_data(xs, ys, R)

#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)

#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs[,full_group_index]
subset_y <- ys[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 3)

new_order <- c(1:3,5,4,6)
ffnet_reorder <- ff_net[new_order, new_order]

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'within_10m_subgroup_network_3m.png'))
visualize_network_matrix_trago(within_group_data$proximity_net, lion_ids[new_order,])
dev.off()

#--------------------------------------------------------

#now want to work out how long they were apart for...

# this code may not be right since the function is for 1Hz data
# ff_data_200 <- detect_fissions_and_fusions(R_inner = 50, R_outer = 200, xs, ys, ts, lion_ids)
# ff_events <- ff_data_200$events_detected
# #removing wonky events
# ff_events <- ff_events[!ff_events$group_A_idxs == "0",]
# analyse_ff_event(7, events = ff_events, xs, ys, ts, max_time = 15)


#made splits_df and merge_df in lion_split_merge_dfs script

load(file = paste0(data_dir, "lion_splits_df.Rdata"))
load(file = paste0(data_dir, "lion_merge_df.Rdata"))

#rbind the dataframes to calculate the duration between fissions and fusions
merge_df$event <- "merge"
splits_df$event <- "split"

merge_df_subset <- merge_df[,c(1,12)]
splits_df_subset <- splits_df[,c(1,12)]

time_diff <- rbind(splits_df_subset, merge_df_subset)
time_diff <- time_diff %>% arrange(-desc(t))

#because the number of splits and merges are different, I'm filtering to times when its a merge then a split (so the duration of a split is from the first split to the next merge)

time_diff$keep <- NA
time_diff$keep2 <- NA

for (i in 1:nrow(time_diff)){
  
  if (time_diff$event[i] == "merge" & time_diff$event[i+1] == "split"){
    time_diff$keep[i] <- "1"
    time_diff$keep2[i+1] <- "1"
  } else {
    time_diff$keep[i] <- "0"
  } 
}


#replace NA's with 0 so I can sum the results 
time_diff$keep2[is.na(time_diff$keep2)] <- 0 
time_diff$keep3 <- as.numeric(time_diff$keep) + as.numeric(time_diff$keep2)
time_diff <- time_diff[time_diff$keep3 == 1,]
time_diff <- time_diff[,c(1,2)]
splits_time_diff <- time_diff[time_diff$event == "split",]
merge_time_diff <- time_diff[time_diff$event == "merge",]

time_diff <- data.frame(splits_t = splits_time_diff$t, merge_t = merge_time_diff$t)
time_diff$diff <- time_diff$splits_t - time_diff$merge_t
time_diff$diff_time_hour <- (time_diff$diff*10)/60
time_diff$diff_time_hour <-  format(round(time_diff$diff_time_hour, 1), nsmall = 1)
time_diff$diff_time_hour <- as.numeric(time_diff$diff_time_hour)


#time_diff_gal was made in the merge_analysis script
hist(time_diff$diff_time_hour, breaks = 100)

png(file = paste0(plot_dir, "/split_duration_hist.png"), width = 1000, height = 600, units = "px")

par(mar=c(8, 8, 8, 8))
hist(time_diff$diff_time_hour, breaks = 80, xlab = "Time between splits and merges (hours)", col = "lightblue3", main = "", cex.lab = 2,cex.axis = 1.5)

dev.off()


#looking at the raw data to see what is happening with splits and merges
colnames(merge_df) <- 1:12
colnames(splits_df) <- 1:12
split_merge_df <- rbind(splits_df, merge_df)
colnames(split_merge_df) <- c("t", "orig_group", "sub1", "sub2", "sub3", "sub4", "sub5", "n_orig", "n_sub1", "n_sub2", "n_sub3", "event")

#order the dataframe by time
split_merge_df <- split_merge_df[(order(split_merge_df$t)),]

#remove duplicates
split_merge_df<- split_merge_df[!duplicated(split_merge_df[,1:5]),]

#change index of the rows to time order
rownames(split_merge_df) <- 1:nrow(split_merge_df)

#add the time of each event
split_merge_df$datetime <- ts[split_merge_df$t]
split_merge_df$round_time <- round_date(split_merge_df$datetime, "1 hour")

split_merge_df$hour <- strftime(split_merge_df$round_time, "%H")
split_merge_df$hour <- as.numeric(split_merge_df$hour)


p1 <- ggplot(split_merge_df, aes(x=hour, fill=event)) +
  geom_histogram(alpha=0.5, position="dodge", bins=24)+
  scale_fill_manual(values = c("turquoise3", "dodgerblue4"))+
  scale_color_brewer(palette="Accent")+ 
  guides(fill="none")+
  theme_classic()

p2 <- ggplot(split_merge_df, aes(x=hour, fill=event)) +
  scale_fill_manual(values = c("turquoise1", "dodgerblue4"))+
  geom_density(alpha = 0.5)+ 
  scale_color_brewer(palette="Accent")+ 
  theme_classic()

p3 <- ggplot(time_diff, aes(x=diff_time_hour))+
  geom_histogram(bins = 150, fill = "purple4")+
  xlab("Time between splits and merges (hours)")+
  theme_classic()

#should look at the time distribution of a fission to a fusion separately to the time from a fusion to a fission

gg <- (p1 | p2) / p3
gg

ggsave(filename = paste0(plot_dir, 'duration_descriptions.png'), plot = gg, width = 10, height = 5, dpi = 300)

#want to plot the distance groups travelled against duration apart

#making new columns to add duration and max dist apart info for plotting
# Making new columns to add duration and max dist apart info for plotting
split_merge_df$dyad_dist_max_1.2 <- NA
split_merge_df$dyad_dist_max_1.3 <- NA
split_merge_df$dyad_dist_max_2.3 <- NA

split_merge_df$event_dur <- NA

# Loop through each row of the dataframe
for (i in 1:(nrow(split_merge_df) - 1)) {
  
  # Print debug information
  cat("Processing row:", i, "\n")
  
  # Check if the current index is less than the total number of rows minus 1
  if (i < nrow(split_merge_df)) {
    
    # Extract time range for the current and next event
    time_range <- split_merge_df$t[i]:split_merge_df$t[i+1]
    
    # Check if the time range is valid
    if (length(time_range) > 0 && split_merge_df$t[i] < split_merge_df$t[i+1]) {
      
      # Group 1 coordinates
      group_1_xs <- xs[split_merge_df$sub1[[i]], time_range, drop = FALSE]
      group_1_ys <- ys[split_merge_df$sub1[[i]], time_range, drop = FALSE]
      
      mean_grp1_xs <- colMeans(group_1_xs, na.rm = TRUE)
      mean_grp1_ys <- colMeans(group_1_ys, na.rm = TRUE)
      
      # Group 2 coordinates
      group_2_xs <- xs[split_merge_df$sub2[[i]], time_range, drop = FALSE]
      group_2_ys <- ys[split_merge_df$sub2[[i]], time_range, drop = FALSE]
      
      mean_grp2_xs <- colMeans(group_2_xs, na.rm = TRUE)
      mean_grp2_ys <- colMeans(group_2_ys, na.rm = TRUE)
      
      # Group 3 coordinates
      group_3_xs <- xs[split_merge_df$sub3[[i]], time_range, drop = FALSE]
      group_3_ys <- ys[split_merge_df$sub3[[i]], time_range, drop = FALSE]
      
      mean_grp3_xs <- colMeans(group_3_xs, na.rm = TRUE)
      mean_grp3_ys <- colMeans(group_3_ys, na.rm = TRUE)
      
      
      # Get distance between centroids of subgroup 1 and 2
      dyad_dist <- sqrt((mean_grp1_xs - mean_grp2_xs)^2 + (mean_grp1_ys - mean_grp2_ys)^2)
      dyad_dist_max_1.2 <- max(dyad_dist, na.rm = TRUE)
      
      # Get distance between centroids of subgroup 1 and 3
      dyad_dist <- sqrt((mean_grp1_xs - mean_grp3_xs)^2 + (mean_grp1_ys - mean_grp3_ys)^2)
      dyad_dist_max_1.3 <- max(dyad_dist, na.rm = TRUE)
      
      # Get distance between centroids of subgroup 2 and 3
      dyad_dist <- sqrt((mean_grp2_xs - mean_grp3_xs)^2 + (mean_grp2_ys - mean_grp3_ys)^2)
      dyad_dist_max_2.3 <- max(dyad_dist, na.rm = TRUE)
      
      
      # Update dataframe with max distance and event duration
      split_merge_df$dyad_dist_max_1.2[i] <- dyad_dist_max_1.2
      split_merge_df$dyad_dist_max_1.3[i] <- dyad_dist_max_1.3
      split_merge_df$dyad_dist_max_2.3[i] <- dyad_dist_max_2.3
      
      #split_merge_df$event_dur[i] <- (split_merge_df$t[i+1] - split_merge_df$t[i]) * 15
      split_merge_df$event_dur[i] <- difftime(split_merge_df$datetime[i + 1], split_merge_df$datetime[i], units = "mins")
      
    } else {
      cat("Invalid time range at row:", i, "\n")
    }
  }
}

#plotting one of the events to see the trajectories
j <- 1:ncol(group_1_xs)

xmax <- max(na.omit(c(group_1_xs, group_2_xs)))
xmin <- min(na.omit(c(group_1_xs, group_2_xs)))
ymax <- max(na.omit(c(group_1_ys, group_2_ys)))
ymin <- min(na.omit(c(group_1_ys, group_2_ys)))

plot(group_1_xs[,j], group_1_ys[,j], type = "l", col = "red", xlim = c(xmax, xmin), ylim = c(ymax,ymin))
points(group_2_xs[,j], group_2_ys[,j], type = "l", col = "blue")


#finding events where split goes to merge to plot the distance/time plot

# Initialize a logical vector to keep desired rows
keep_rows <- rep(FALSE, nrow(split_merge_df))

# Initialize a variable to track the start of a "split" sequence
split_started <- FALSE
split_index <- NA

# Loop through the dataframe to identify the first "split" and the first "merge" in each sequence
for (i in 1:(nrow(split_merge_df) - 1)) {
  if (split_merge_df$event[i] == "split" && !split_started) {
    # Mark the first split in the sequence
    split_started <- TRUE
    split_index <- i
  }
  
  if (split_started && split_merge_df$event[i + 1] == "merge") {
    # Mark the first merge after the split sequence
    keep_rows[split_index] <- TRUE
    keep_rows[i + 1] <- TRUE
    
    # Reset tracking variables for the next sequence
    split_started <- FALSE
    split_index <- NA
  }
}

# Filter the dataframe
filtered_split_merge_df <- split_merge_df[keep_rows, ]

#getting the duration in days for plotting
filtered_split_merge_df$event_dur_day <- (filtered_split_merge_df$event_dur/60)/24

#Duration groups were split
ggplot(filtered_split_merge_df[filtered_split_merge_df$event == "split",], aes(x = event_dur_day, y = dyad_dist_max_1.2)) +
  geom_point() +
  xlab("Duration (days)") +
  ylab("Distance Apart (m)") +
  theme_classic()

#durations within the day


#I'm not sure the dyadic distances are correct here, so need to make gganimate plots for each event to see what's happening.... especially for the merges which have a huge dyadic distance e.g on 2023-11-07 17:45:00
#could be GPS error! for i = 333 (same as time written here)

#plotting each dyadic distance over time


# Convert matrices to data frames for easier handling
xs_df <- as.data.frame(xs)
ys_df <- as.data.frame(ys)

colnames(xs_df) <- ts
colnames(ys_df) <- ts


# Add an 'Individual' column
xs_df$Individual <- 1:nrow(xs_df)
ys_df$Individual <- 1:nrow(ys_df)

# Convert to long format
xs_long <- xs_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "X_UTM") 
ys_long <- ys_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "Y_UTM") 

# Merge the long data frames
utm_long <- xs_long %>%
  inner_join(ys_long, by = c("Individual", "Time"))

# Function to calculate dyadic distances
calculate_dyadic_distances <- function(data) {
  distances <- data %>%
    inner_join(data, by = "Time", suffix = c(".1", ".2")) %>%
    filter(Individual.1 < Individual.2) %>%
    mutate(Distance = sqrt((X_UTM.1 - X_UTM.2)^2 + (Y_UTM.1 - Y_UTM.2)^2)) %>%
    select(Time, Individual.1, Individual.2, Distance)
  
  return(distances)
}

# Calculate dyadic distances
dyadic_distances <- calculate_dyadic_distances(utm_long)

# Function to generate pairwise data frame for a specific individual
generate_distance_df <- function(individual_id) {
  individual_data <- dyadic_distances %>%
    filter(Individual.1 == individual_id | Individual.2 == individual_id) %>%
    mutate(Other_Individual = if_else(Individual.1 == individual_id, Individual.2, Individual.1)) %>%
    select(Time, Other_Individual, Distance) %>%
    mutate(Focal_Individual = individual_id)
  
  return(individual_data)
}

# Generate distance data for all individuals
distance_data_list <- lapply(1:nrow(xs), generate_distance_df)

# Combine all distance data into one data frame
combined_distance_data <- bind_rows(distance_data_list)

# Replace ID numbers with actual names
combined_distance_data$Focal_Individual <- lion_ids$name[combined_distance_data$Focal_Individual]
combined_distance_data$Other_Individual <- lion_ids$name[combined_distance_data$Other_Individual]

# Fill missing data using linear interpolation and carry forward/backward for leading/trailing NAs
combined_distance_data <- combined_distance_data %>%
  group_by(Focal_Individual, Other_Individual) %>%
  mutate(Distance = zoo::na.locf(zoo::na.locf(zoo::na.approx(Distance, na.rm = FALSE), na.rm = FALSE), fromLast = TRUE)) %>%
  ungroup()

# Define the window size for the rolling mean
window_size <- 300

# Add a new column for the rolling mean, calculating it for each dyad separately
combined_distance_data <- combined_distance_data %>%
  group_by(Focal_Individual, Other_Individual) %>%
  mutate(Rolling_Mean_Distance = zoo::rollmean(Distance, k = window_size, fill = NA, align = "center")) %>%
  ungroup()

# Plotting function for individual dyadic distances
plot_individual_distances <- function(individual_id, data) {
  individual_plot_data <- data %>%
    filter(Focal_Individual == individual_id) %>%
    mutate(Time = as.POSIXct(Time))  # Ensure Time column is in POSIXct format
  
  # Get unique individuals for consistent colors
  unique_individuals <- unique(individual_plot_data$Other_Individual)
  # Generate colors for each unique individual
  colors <- rainbow(length(unique_individuals))
  # Create a named vector with colors for each individual
  color_mapping <- setNames(colors, unique_individuals)
  
  ggplot(individual_plot_data, aes(x = Time, y = Distance)) +
    geom_line(aes(color = Other_Individual, group = Other_Individual), size = 0.8, alpha = 0.8) +
    geom_line(aes(x = Time, y = Rolling_Mean_Distance, color = Other_Individual, group = Other_Individual), size = 6, alpha = 0.3) +
    scale_color_manual(values = color_mapping) +
    labs(title = paste("Dyadic Distances for Individual", individual_id),
         x = "Time",
         y = "Distance (meters)",
         colour = "Other Individual") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 12),
          axis.ticks.y = element_blank(),
          legend.position = "bottom") +
    scale_x_datetime(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") # Adjust the datetime scale
}

# Example usage
individual_id <- "Sarabi"
plot_individual_distances(individual_id, combined_distance_data)

# Loop through each individual and generate plots
individual_ids <- unique(combined_distance_data$Focal_Individual)
for (individual_id in individual_ids) {
  plot <- plot_individual_distances(individual_id, combined_distance_data)
  ggsave(filename = paste0("Dyadic_Distances_", individual_id, ".png"), plot = plot, width = 10, height = 6)
}

plot_individual_distances(individual_id, combined_distance_data)

# Example of how to filter and plot for one individual
individual_id <- "Sarabi"
individual_id <- "Kiara"
individual_id <- "Nala"


# Loop through each individual and save it
for (individual_id in unique(combined_distance_data$Focal_Individual)) {
  # Create the plot for the current individual
  plot <- plot_individual_distances(individual_id, combined_distance_data)
  ggsave(filename = paste0(plot_dir, "individual_", individual_id, "_distance_plot.png"), plot, width = 20, height = 6, units = "in")
  
}


#add these to doc - all data, do for 0-100, log of overall
#add speed color to plot - are they resting when 10m apart?
ggplot(combined_distance_data, aes(Distance))+
  geom_histogram()+
  facet_wrap(~Focal_Individual)+
  xlim(0, 1000)



hist(combined_distance_data$Distance, breaks = 100000, xlim = c(0, 100))

