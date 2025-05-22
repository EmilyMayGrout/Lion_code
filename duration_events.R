
#getting correct event duration distributions
#will measure it based on the time each individual is with each other individual


setwd("C:/Users/egrout/Dropbox/lion/Lion_code")

#--------PARAMS-------
data_dir <- "../data/processed/"
code_dir <- '/lion/Lion_code/'
plot_dir <- '../results/level0/'

#---------------------
#Get functions 
source('coati_function_library_V1.R')
source('lion_functions.R')

library(hms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(tidyverse)
library(stringr)
library(combinat)
library(purrr)
library(tibble)
library(patchwork)
library(grid)
library(reshape2)

#load data
load("../data/processed/lion_xy_15min_level0_2yr.RData")
load("../data/processed/lion_ids.RData")

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#subset to times when the data is every 15 mins
non_na_counts <- colSums(!is.na(xs))

drops <- which(diff(non_na_counts) < -4)

#find timepoint for first drop as this is likely where the 15 minute sampling was changed to 2hours
ts[drops[1]]

# Find date that Simba died 
# Find the index of the last non-NA value
last_valid_index <- max(which(!is.na(xs[3,])))

# Get the corresponding timestamp
ts[last_valid_index]

#insert NA to after Simba died
xs[3, last_valid_index:ncol(xs)] <- NA
ys[3, last_valid_index:ncol(xs)] <- NA

#filter gps data to when it's 15 minute sampling regime
#dates now range from 2023-05-04 to 2023-11-22 (1.5 years of 15 mins)
xs <- xs[,1:drops[1]]
ys <- ys[,1:drops[1]]
ts <- ts[1:drops[1]]


#------------------------

include_singletons <-  F

#------------------------

R = 100
subgroup_data <- get_subgroup_data(xs, ys, R)
#function code in lion_functions
splits_df <- detect_splits(subgroup_data$ind_subgroup_membership)

#because merge time occurs from before -> during and split time is from during -> after, I will change the event time of the split to one timestamp later, because we have the same times in split and merge events which isn't actually occurring at the same time
splits_df$t <- splits_df$t + 1


#adding the events when it's singletons
splits_df$singleton_event <- "TRUE"

for (i in 1:nrow(splits_df)){
  
  split_i <- splits_df[i,]
  total_size <- split_i$n_orig 
  sub_sizes <- split_i[9:13]
  
  if(sum(sub_sizes > 1, na.rm = TRUE) >= 2 ){
    
    splits_df$singleton_event[i] <- "FALSE"
  }
}

table(splits_df$singleton_event)

merges_df <- detect_merges(subgroup_data$ind_subgroup_membership)

#adding the events when it's singletons
merges_df$singleton_event <- "TRUE"

for (i in 1:nrow(merges_df)){
  
  merge_i <- merges_df[i,]
  total_size <- merge_i$n_orig 
  orig_sizes <- merge_i[9:13]
  
  if(sum(orig_sizes > 1, na.rm = TRUE) >= 2 ){
    
    merges_df$singleton_event[i] <- "FALSE"
  }
}

table(merges_df$singleton_event)


if (include_singletons == F){
  
  file <- "without_singletons"
  merges_df <- merges_df[merges_df$singleton_event == "FALSE", ]
  splits_df <- splits_df[splits_df$singleton_event == "FALSE", ]
  
}else if(include_singletons == T){
  
  file <- "with_singletons"
  merges_df <- merges_df
  splits_df <- splits_df
  
}


#cut them to events which were within the 15 mins scheduling
splits_df <- splits_df[splits_df$t < drops[1],]
merges_df <- merges_df[merges_df$t < drops[1],]


#rbind the dataframes to calculate the duration between fissions and fusions
merges_df$event <- "merge"
splits_df$event <- "split"



#looking at the raw data to see what is happening with splits and merges
colnames(merges_df) <- 1:14
colnames(splits_df) <- 1:14

split_merge_df <- rbind(splits_df, merges_df)
colnames(split_merge_df) <- c("t", "orig_group", "sub1", "sub2", "sub3", "sub4", "sub5", "n_orig", "n_sub1", "n_sub2", "n_sub3","n_sub4","n_sub5","singleton_event", "event")

#order the dataframe by time
split_merge_df <- split_merge_df[(order(split_merge_df$t)),]

#remove duplicates
split_merge_df <- split_merge_df[!duplicated(split_merge_df[,1:5]),]

#change index of the rows to time order
rownames(split_merge_df) <- 1:nrow(split_merge_df)

#add the time of each event
split_merge_df$datetime <- ts[split_merge_df$t]
split_merge_df$round_time <- round_date(split_merge_df$datetime, "1 hour")
#split_merge_df$hour <- lubridate::hour(split_merge_df$hour)



#create a list for each dyad possibility
all_dyads <- expand.grid(1:n_inds,1:n_inds)
colnames(all_dyads) <- c("other", "focal")
all_dyads <- all_dyads[all_dyads$focal != all_dyads$other,]

#go through each dyad in the list and find when they are first together in the dataframe - when they are together in merge. Then find the next time point they are in separate columns for split, then get the time point when they togeth

#extract events when these individuals are present together for merges, and separate for splits 

list_dyad <- list()

for (dyad in 1:nrow(all_dyads)){

  focal_ind <- all_dyads[dyad,2]
  other_ind <- all_dyads[dyad,1]
  
  dyad_name <- paste(focal_ind, other_ind, sep = "_")
  event_dyad <- data.frame()  # empty data frame to collect events
  
  for (event in 1:nrow(split_merge_df)) {
    
    #find all events where focal and other were in separate columns in a merge event (when they came together) and in separate columns in a split event (when they split)
    
    focal_sub <- names(which(unlist(split_merge_df[event, c("sub1", "sub2", "sub3", "sub4", "sub5")]) ==  focal_ind))
    focal_sub <- substr(focal_sub, 1,4)
    other_sub <- names(which(unlist(split_merge_df[event, c("sub1", "sub2", "sub3", "sub4", "sub5")]) ==  other_ind))
    other_sub <- substr(other_sub, 1,4)
    
    if (length(focal_sub) > 0 && length(other_sub) > 0 && focal_sub != other_sub) {
      
      #put the dataframe in the list called list_dyad
      event_dyad <- rbind(event_dyad, split_merge_df[event, ])
    } 
  }
  # Store results in the list by dyad name
  list_dyad[[dyad_name]] <- event_dyad
 }



#go through each dyads events to calculate duration between events, only including events from last split to first merge and last merge to first split 

for (dyad_name in names(list_dyad)) {
  
  events <- list_dyad[[dyad_name]]
  
  # Skip if fewer than 2 events
  if (nrow(events) < 2) {
    events$duration <- NA
    list_dyad[[dyad_name]] <- events
    next
  }
  
  # Ensure chronological order
  events <- events[order(events$t), ]
  
  # Initialize duration column
  events$duration <- NA
  
  # Loop through rows and calculate durations
  for (i in 1:(nrow(events) - 1)) {
    event_now <- events[i, ]
    event_next <- events[i + 1, ]
    
    # split → merge OR merge → split
    if ((event_now$event == "split" && event_next$event == "merge") ||
        (event_now$event == "merge" && event_next$event == "split")) {
      
      # Store duration in the *next* row
      events$duration[i + 1] <- event_next$t - event_now$t
    }
  }
  
  # Update the list
  list_dyad[[dyad_name]] <- events
}






#go through each dyad and each event to get the maximum distance they were apart between each event

for (dyad_name in names(list_dyad)) {
  
  dyad_df <- list_dyad[[dyad_name]]
  
  # Extract focal and other from dyad name
  inds <- unlist(str_split(dyad_name, "_"))
  focal_ind <- as.numeric(inds[1])
  other_ind <- as.numeric(inds[2])
  
  # Create column
  dyad_df$max_dyad_dist <- NA
  
  # Loop through each event row, starting from the second (since we need a previous row)
  for (i in 2:nrow(dyad_df)) {
    time_start <- as.numeric(dyad_df$t[i - 1])
    time_end <- as.numeric(dyad_df$t[i])
    
    # If either time point is NA, skip
    if (is.na(time_start) || is.na(time_end)) next
    
    # Get sequence of time points from t_prev to t_current (inclusive)
    time_seq <- time_start:time_end
    
    # Extract coordinates across the duration
    xs_focal <- xs[focal_ind, time_seq]
    ys_focal <- ys[focal_ind, time_seq]
    xs_other <- xs[other_ind, time_seq]
    ys_other <- ys[other_ind, time_seq]
    
    # Compute Euclidean distances
    dist_seq <- sqrt((xs_focal - xs_other)^2 + (ys_focal - ys_other)^2)
    
    # Store max distance (handle NA safely)
    dyad_df$max_dyad_dist[i] <- max(dist_seq, na.rm = TRUE)
  }
  
  list_dyad[[dyad_name]] <- dyad_df
}

#plotting duration~dyadic distance for events which are splits to merges, coloured by the dyad 

t <- list_dyad[[1]]


# Combine all dyad dataframes into one dataframe
plot_df <- bind_rows(
  lapply(names(list_dyad), function(dyad_name) {
    df <- list_dyad[[dyad_name]]
    df$dyad <- dyad_name
    return(df)
  })
)


plot_df$hours <- plot_df$duration/4

#remove durations which are 0
plot_df <- plot_df[plot_df$hours != 0,]

# Filter only rows that are merge events AND have non-NA durations/distances
plot_df <- plot_df %>%
  filter(event == "merge", !is.na(duration), !is.na(max_dyad_dist))

# Plot: duration vs. distance, coloured by dyad
g <- ggplot(plot_df, aes(x = max_dyad_dist, y = hours, color = dyad)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE)+
  labs(
    title = "Duration vs. Max Dyadic Distance (Split to Merge Events)",
    x = "Max Dyadic Distance",
    y = "Duration (hours)",
    color = "Dyad"
  ) +
  theme_classic()
g


#looking at the duration and max dyadic distance for events within 1 day
plot_df_withinday <- plot_df[plot_df$hours <= 24, ]

p1 <- ggplot(plot_df_withinday, aes(x = max_dyad_dist, y = hours)) +
  geom_point(color = "slateblue1", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = " < 24 hours") +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )

#looking at the duration and max dyadic distance for events over 1 day
plot_df_overday <- plot_df[plot_df$hours > 24, ]

p2 <- ggplot(plot_df_overday, aes(x = max_dyad_dist, y = hours)) +
  geom_point(color = "slateblue3", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = " < 24 hours") +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )
p2

p3 <- ggplot(plot_df, aes(x = max_dyad_dist, y = hours)) +
  geom_point(color = "slateblue4", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = "All events") +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )


combined_plot <- (p1 + p2) / p3  + plot_layout(heights = c(1, 1))  # Equal heights


png(file = paste0(plot_dir, "/duration_dyadicdist_2yr.png"), width = 10, height = 6, units = "in", res = 300)

grid.newpage()
grid.draw(patchworkGrob(combined_plot))
grid.text("Event Duration (hours)", x = unit(0.02, "npc"), rot = 90, gp = gpar(fontsize = 14))

dev.off()


plot_df <- plot_df %>% separate(dyad, into = c("ind1", "ind2"), sep = "_", convert = TRUE)
plot_df <- plot_df %>% mutate(id_min = pmin(ind1, ind2), id_max = pmax(ind1, ind2))

summary_df <- plot_df %>%
  group_by(id_min, id_max) %>%
  summarise(
    mean_duration = mean(duration, na.rm = TRUE),
    mean_dist = mean(max_dyad_dist, na.rm = TRUE),
    .groups = "drop"
  )

# Fill missing dyads with NA
all_ids <- sort(unique(c(summary_df$id_min, summary_df$id_max)))

# Duration matrix
duration_matrix <- summary_df %>%
  complete(id_min = all_ids, id_max = all_ids) %>%
  pivot_wider(names_from = id_max, values_from = mean_duration) %>%
  arrange(id_min)

# Distance matrix
dist_matrix <- summary_df %>%
  complete(id_min = all_ids, id_max = all_ids) %>%
  pivot_wider(names_from = id_max, values_from = mean_dist) %>%
  arrange(id_min)

# Duration heatmap
g1 <- summary_df %>%
  ggplot(aes(x = factor(id_max), y = factor(id_min), fill = mean_duration)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(title = "Mean Duration by Dyad", x = "Individual (max)", y = "Individual (min)") +
  theme_classic()

# Distance heatmap
g2 <- summary_df %>%
  ggplot(aes(x = factor(id_max), y = factor(id_min), fill = mean_dist)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", na.value = "grey90") +
  labs(title = "Mean Distance by Dyad", x = "Individual (max)", y = "Individual (min)") +
  theme_classic()

gg <- g1 + g2

ggsave(filename = paste0(plot_dir, 'mean_duration_distance_', file, '.png'), plot = gg, width = 14, height = 7, dpi = 300)



