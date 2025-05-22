#Questions interested for the lion dataset:
#1. duration of splits
#2. who do they split with
#3. what are they doing when split?
#4. how far do they travel apart?
#5. how often do they come together?

setwd("C:/Users/egrout/Dropbox/lion/Lion_code")

#--------PARAMS-------
data_dir <- "../data/processed/"
code_dir <- '/lion/Lion_code/'
plot_dir <- '../results/level0/'

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(tidygraph)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(ggplot2)
library(vioplot)
library(plotly)
library(patchwork)
library(grid)
library(reshape2)
library(igraph)
library(ggraph)
library(gganimate)
library(purrr)
library(av)

#read in library of functions
source('coati_function_library_V1.R')
source('lion_functions.R')

#load data
load("../data/processed/lion_xy_15min_level0_2yr.RData")
load("../data/processed/lion_ids.RData")


#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}


make_plots <- F #making the plots takes a bit of time for code to run
include_singletons <- T



#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#subset to times when the data is every 15 mins
non_na_counts <- colSums(!is.na(xs))

plot_samp_change <- F

if (plot_samp_change == T){

threshold <- 0.2 * nrow(xs)  # e.g., 10% of individuals
high_sampling_cols <- which(non_na_counts >= threshold)

#find consecutive indices (i.e., columns) where this is true. extracts runs of 15-min data lasting at least 1 day (96 timestamps)
runs <- rle(diff(high_sampling_cols) == 1)
run_lengths <- runs$lengths
run_starts <- cumsum(c(1, head(run_lengths, -1)))

# Find long runs
long_runs <- which(runs$values & run_lengths >= 96)
valid_indices <- unlist(lapply(long_runs, function(i) {
  start <- run_starts[i]
  seq(high_sampling_cols[start], high_sampling_cols[start] + run_lengths[i])
}))

#plot number of NAs over time to find changes in sampling regimes
plot(ts, non_na_counts, type = "l", 
     main = "Non-NA GPS Points Over Time",
     xlab = "Time", ylab = "Number of GPS Points",
     col = "steelblue")

lines(ts, zoo::rollmean(non_na_counts, k = 10, fill = NA), col = "red", lwd = 2)
legend("bottomleft", legend = c("Raw", "Smoothed"), col = c("steelblue", "red"), lty = 1)

drops <- which(diff(non_na_counts) < -4)
points(ts[drops], non_na_counts[drops], col = "yellow", pch = 16)
}

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


#before subsetting, this is the datetime range of all the data
#firsttime <- as.POSIXct('2023-05-04 00:00', tz = 'UTC')
#lasttime <-  as.POSIXct('2025-04-12 12:00', tz = 'UTC')

#filter gps data to when it's 15 minute sampling regime
#dates now range from 2023-05-04 to 2023-11-22 (1.5 years of 15 mins)
xs <- xs[,1:drops[1]]
ys <- ys[,1:drops[1]]
ts <- ts[1:drops[1]]

#----------------
#looking at the sub grouping patterns

#list of Rs
Rs1 <- c(10,20,30,40,50)
Rs2 <- c(100,200,300,400,500)
Rs3 <- c(1000,2000,3000,4000,5000)
Rs4 <- c(10000,20000,30000,40000,50000)

Rs <- list(Rs1, Rs2, Rs3, Rs4)

#choose your color scheme :) 
color <- c("gold", "gold2", "gold3", "gold4")
color <- c("slateblue1","slateblue2","slateblue3","slateblue4" )

#Number of sub groups when the radius is changed (graph put in dropbox results folder)


if(make_plots == T){
for(j in seq_along(Rs)){
  
  Rs_j <- Rs[[j]]
  Rs_name <- paste0(Rs_j[1], "-", Rs_j[5])

png(height = 1700, width = 500, units = 'px', filename = paste0(plot_dir,Rs_name,'m_2yr.png'))
par(mfrow=c(5,1), mar = c(8,9,12,1), mgp=c(6,1,0)) #(bottom, left, top, right)

for (i in 1:length(Rs_j)){
  
  R <- Rs_j[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs_j)){
    xlab <- 'Number of subgroups'
  }
  
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = color[j], breaks = seq(.5,6.5,1), cex.lab = 4, cex.main = 4, cex.axis=4, freq = FALSE, ylim=c(0,.7), xlim = c(0.5, 6.5), yaxt = "n", axes = F)
  axis(2, at = c(0,0.3,0.6), cex.axis = 4, las = 1, pos = 0.5)
  axis(1, cex.axis = 4, las = 1, pos = 0, padj = 0.6)
      }

  dev.off()
  print(paste0("done plot ", j[1]))

  }
 }

#-----------------------------------

#looking at which individuals split with one another

R = 100
subgroup_data <- get_subgroup_data(xs, ys, R)

if(make_plots == T){
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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_100m_2yr.png'))

visualize_network_matrix_trago(ffnet_reorder, lion_ids[new_order,])
dev.off()

#it looks like Kiara and Sarafina split from the other group members and travelled atleast 50km from them

#how much missing data is there?

png(height = 500, width = 1000, units = 'px', filename = paste0(plot_dir, "number_tracked_2yr.png"))
par(mfrow=c(1,2), mar = c(10,9,2,1)) #c(bottom, left, top, right)
sum_tracked <- colSums(!is.na(xs))
hist(sum_tracked[ !(sum_tracked==0) ], main = "", xlab = "Number of individuals tracked", ylab = "Frequency", col = "gold2", cex.lab = 2, cex.axis = 2)

each_sum <- data.frame(sum = rowSums(!is.na(xs)))
barplot(each_sum$sum, names.arg = lion_ids$name, las=2, col = "gold3", ylab = "Number of GPS points",  cex.lab = 2, cex.axis = 2, cex.names=2, mgp=c(5,0.2,0))

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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'within_10m_subgroup_network_3m_2yr.png'))
visualize_network_matrix_trago(within_group_data$proximity_net, lion_ids[new_order,])
dev.off()

}

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------


#now want to work out how long they were apart for...
# Detect events

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


#distribution of subgroup sizes for splits and merges
#hist(cbind(splits_df$n_sub1, splits_df$n_sub2, splits_df$n_sub3), breaks = 4)
#hist(cbind(merges_df$n_sub1, merges_df$n_sub2, merges_df$n_sub3), breaks = 4)

#cut them to events which were within the 15 mins scheduling
splits_df <- splits_df[splits_df$t < drops[1],]
merges_df <- merges_df[merges_df$t < drops[1],]


#rbind the dataframes to calculate the duration between fissions and fusions
merges_df$event <- "merge"
splits_df$event <- "split"

merges_df_subset <- merges_df[,c(1,15)]
splits_df_subset <- splits_df[,c(1,15)]

time_diff <- rbind(splits_df_subset, merges_df_subset)
time_diff <- time_diff %>% arrange(-desc(t))

#because the number of splits and merges are different, I'm filtering to times when its a split then a merge (so the duration of a split is from the first split to the next merge), or the times when its a merge then a split (as these distributions of times are different depending on the ordering here)

# Loop through rows (except the last to avoid out-of-bounds)
#make copy to store merge the split data
time_diff$keep <- rep("0", nrow(time_diff))
time_diff$keep2 <- rep("0", nrow(time_diff))
#time_diff will store split them merge
split_then_merge <- time_diff
merge_then_split <- time_diff

for (i in 1:(nrow(split_then_merge) - 1)) {
  if (split_then_merge$event[i] == "split" & split_then_merge$event[i + 1] == "merge") {
    split_then_merge$keep[i] <- "1"
    split_then_merge$keep2[i + 1] <- "1"
  }
}

for (i in 1:(nrow(merge_then_split) - 1)) {
  if (merge_then_split$event[i] == "merge" & merge_then_split$event[i + 1] == "split") {
    merge_then_split$keep[i] <- "1"
    merge_then_split$keep2[i + 1] <- "1"
    }
  }
  
time_diff_list <- list(split_then_merge, merge_then_split)

# Process both sets
split_merge_change_order <- lapply(time_diff_list, function(td) {
  td$keep3 <- as.numeric(td$keep) + as.numeric(td$keep2)
  td <- td[td$keep3 == 1, ]
  td <- td[, c("event", "t")]  # Keep relevant columns
  
  splits_td <- td[td$event == "split", ]
  merges_td <- td[td$event == "merge", ]
  
  # Align by row
  td_out <- data.frame(
    splits_t = splits_td$t,
    merges_t = merges_td$t
  )
  
  # Compute differences
  td_out$diff <- abs(td_out$merges_t - td_out$splits_t)
  td_out$diff_time_hour <- round((td_out$diff * 10) / 60, 1)
  
  
  return(td_out)
})

# Optional: name them for clarity
names(split_merge_change_order) <- c("split_then_merge", "merge_then_split")


if(make_plots == T){

#time_diff_gal was made in the merge_analysis script
#hist(time_diff$diff_time_hour, breaks = 100)

png(file = paste0(plot_dir, "/split_duration_hist_2yr_", file, ".png"), width = 1000, height = 600, units = "px")

par(mar=c(8, 8, 8, 8))
hist(split_merge_results[[2]]$diff_time_hour, breaks = 80, xlab = "Time between splits and merges (hours)", col = "lightblue3", main = "", cex.lab = 2, cex.axis = 1.5)

dev.off()

}

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

split_merge_df$hour <- strftime(split_merge_df$round_time, "%H")
split_merge_df$hour <- as.numeric(split_merge_df$hour)

if(make_plots == T){

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

p3 <- ggplot(split_merge_change_order[[1]], aes(x=diff_time_hour))+
  geom_histogram(bins = 80, fill = "slateblue")+
  xlab("Time from last split to next merge (hours)")+
  theme_classic()

p4 <- ggplot(split_merge_change_order[[2]], aes(x=diff_time_hour))+
  geom_histogram(bins = 80, fill = "slateblue4")+
  xlab("Time from last merge to next split (hours)")+
  theme_classic()


#should look at the time distribution of a fission to a fusion separately to the time from a fusion to a fission

gg <- (p1 | p2) / p3 / p4
gg

ggsave(filename = paste0(plot_dir, 'duration_descriptions_2yr_', file, '.png'), plot = gg, width = 10, height = 7, dpi = 300)

}

#plot the before and after of each split and merge event to save in a folder to get an idea of what occurred

# Define a semi-transparent color for the arrows
arrow_col <- rgb(red = 0.9, green = 0.6, blue = 0.1,  alpha = 0.3)  # RGB values + alpha (0 to 1)
arrow_col2 <- rgb(red = 0.1, green = 0.3, blue = 0.5,  alpha = 0.3)  # RGB values + alpha (0 to 1)


#this code takes a while to run
if(make_plots == T){

#change the time of the split event back to its original time index
  split_merge_df$t[split_merge_df$event == "split"] <- split_merge_df$t[split_merge_df$event == "split"] -1
  
  
for(i in 1:nrow(split_merge_df)){
  
  event_type <- split_merge_df$event[i]
  event_time <- split_merge_df$t[i]
  
  inds_involved <- unlist(split_merge_df$orig_group[i])
  
  xs_event <- xs[inds_involved, event_time]
  xs_bef_event <- xs[inds_involved, event_time - 1]
  xs_aft_event <- xs[inds_involved, event_time + 1]
  
  ys_event <- ys[inds_involved, event_time]
  ys_bef_event <- ys[inds_involved, event_time - 1]
  ys_aft_event <- ys[inds_involved, event_time + 1]
  
  xlim <- range(c(xs_event, xs_bef_event, xs_aft_event), na.rm = TRUE)
  ylim <- range(c(ys_event, ys_bef_event, ys_aft_event), na.rm = TRUE)
  
  png(file = paste0(plot_dir, "movement_tracks/splitmerge_events/", event_type, "plot_event", event_time, ".png"),width = 6, height = 6, units = "in", res = 300)
  
  par(mar=c(4.5, 4.1, 3.5, 7.1), xpd=TRUE) #bottom, left, top, right
  
  plot(xs_event, ys_event, xlim = xlim, ylim = ylim, col = "yellow2", pch = 16, cex = 2,
       main = paste(event_type, "event at t =", event_time),
       xlab = "X", ylab = "Y", asp = 1)
  
  points(xs_bef_event, ys_bef_event, col = "orange1", pch = 16, cex = 2)
  points(xs_aft_event, ys_aft_event, col = "slateblue", pch = 16, cex = 2)
  
  # Arrows: before → event
  arrows(xs_bef_event, ys_bef_event, xs_event, ys_event, col = arrow_col, length = 0.15, angle = 20, lwd = 2)
  
  # Arrows: event → after
  arrows(xs_event, ys_event, xs_aft_event, ys_aft_event, col = arrow_col2, length = 0.15, angle = 20, lty = 1, lwd = 2)
  
  # Add legend
  legend("topright",
         inset=c(-0.3,0),
         legend = c("Before", "During", "After"),
         col = c("orange1", "yellow2", "slateblue"),
         pch = 16,
         title = "Time")
  
  
  # ---- Add scale bar ----
  # Define scale bar length in same units as xs/ys (e.g., 100 meters)
  scale_length <- 10
  
  # Place scale bar 5% from bottom and 20% from right
  x_pad <- 0.2 * diff(xlim)   # 20% from right
  y_pad <- 0.05 * diff(ylim)  # 5% from bottom
  
  x_end <- xlim[2] - x_pad
  x_start <- x_end - scale_length
  y_pos <- ylim[1] + y_pad
  
  # Draw scale bar
  segments(x_start, y_pos, x_end, y_pos, lwd = 2, col = "black")
  text(x = (x_start + x_end)/2, y = y_pos + 0.02 * diff(ylim),
       labels = paste0(scale_length, " m"), cex = 0.8)
    # ------------------------
  
  
  dev.off()
 }
}


#which individuals often leave one another - are there individuals who often break association? Here I will go through each split event, and give a score for each focal to each other individual, 1 if they are in the same subgroup and 0 if they are in a different subgroup, run through all events and tally up these scores

event_type <- "merge"

# Get all split events
event_all <- split_merge_df[split_merge_df$event == event_type,]
n_events <- nrow(event_all)
n_inds <- nrow(xs)
#all_inds <- as.character(1:n_inds)
all_inds <- lion_ids$tag_id

# Create a 3D array to store one matrix per event
event_association_array <- array(NA, dim = c(n_inds, n_inds, n_events),
                                 dimnames = list(all_inds, all_inds, NULL))

for (i in 1:n_events) {
  event <- event_all[i, ]
  orig_inds <- event$orig_group[[1]]
  
  # Get all subgroups, remove NULLs or NAs
  subgroups <- list(event$sub1[[1]], event$sub2[[1]], event$sub3[[1]], 
                    event$sub4[[1]], event$sub5[[1]])
  subgroups <- subgroups[!sapply(subgroups, function(g) is.null(g) || all(is.na(g)))]
  
  # Initialize matrix with NA
  assoc_matrix <- matrix(NA, nrow = n_inds, ncol = n_inds, dimnames = list(all_inds, all_inds))
  
  # For each pair in the original group, score whether they ended up together
  for (ind1 in orig_inds) {
    for (ind2 in orig_inds) {
      if (ind1 == ind2) {
        assoc_matrix[ind1, ind2] <- NA
      } else {
        together <- any(sapply(subgroups, function(g) all(c(ind1, ind2) %in% g)))
        assoc_matrix[ind1, ind2] <- ifelse(together, 1, 0)
      }
    }
  }
  event_association_array[,,i] <- assoc_matrix
}


# Sum the array across events, ignoring NAs
overall_event_association <- apply(event_association_array, c(1,2), function(x) sum(x, na.rm = TRUE))

diag(overall_event_association) <- NA

# Convert to long format
assoc_df <- melt(overall_event_association, varnames = c("Ind1", "Ind2"), value.name = "Association")
# remove NA entries (e.g., self-pairs)
assoc_df <- assoc_df[!is.na(assoc_df$Association), ]

if(make_plots == T){

#overall plot showing how often each dyad (pair of individuals) ended up together across splits
g1 <- ggplot(assoc_df, aes(x = as.factor(Ind1), y = as.factor(Ind2), fill = Association)) +
  geom_tile() +
  scale_fill_gradient(low = "slateblue4", high = "yellow", na.value = "white") +
  theme_minimal() +
  labs(title = paste0("Overall ", event_type, " Association Matrix"),
       x = "Individual", y = "Individual", fill = "Association") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()


# (Optional) Normalize by number of events each dyad was involved in
dyad_counts <- apply(!is.na(event_association_array), c(1,2), sum)
normalized_event_association <- overall_event_association / dyad_counts
diag(normalized_event_association) <- NA
# Convert to long format
norm_assoc_df <- melt(normalized_event_association, varnames = c("Ind1", "Ind2"), value.name = "Association")
# remove NA entries (e.g., self-pairs)
norm_assoc_df <- norm_assoc_df[!is.na(norm_assoc_df$Association), ]

# normalised plot
g2 <- ggplot(norm_assoc_df, aes(x = as.factor(Ind1), y = as.factor(Ind2), fill = Association)) +
  geom_tile() +
  scale_fill_gradient(low = "slateblue4", high = "yellow", na.value = "white") +
  theme_minimal() +
  labs(title = paste0("Normalized ", event_type, " Association Matrix"),
       x = "Individual", y = "Individual", fill = "Association") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

gg <- g1 + g2
gg

ggsave(filename = paste0(plot_dir, event_type, '_association_', file, '.png'), plot = gg, width = 12, height = 5, dpi = 300)

}

run_consistency <- F
if(run_consistency == T && make_plots == T){

#function that takes in a splits dataframe and outputs a p_dyad_together matrix (probability of splitting together for each dyad)
p_dyad_together <- get_p_dyad_together(splits_df_local = splits_all, n_inds_local = n_inds)

#Compute consistency based on consistency metric (see coati_function_library)
consistency_data <- get_consistency(p_dyad_together)

#randomize splits and recompute consistency metric
set.seed(24)
n_rands <- 1000
rando_consistencies <- rep(NA, n_rands)
for (i in 1:n_rands){
  
  rando <- randomise_splits(splits_all)
  rando_dyad <- get_p_dyad_together(splits_df_local = rando, n_inds_local = n_inds)
  consist <- get_consistency(rando_dyad)
  rando_consistencies[i] <- consist
  
}


# Create data frame for the histogram
df <- data.frame(rando_consistencies = rando_consistencies)

if(include_singletons == F){
  y_start <- 0
  y_end <- 180
  breaks = seq(0.2, 0.35, 0.001)
  
}else{
  y_start <- 0
  y_end <- 250
  breaks = seq(0.05, 0.25, 0.002)
}


gg <- ggplot(df, aes(x = rando_consistencies)) +
  geom_histogram(breaks = breaks,
                 fill = "slateblue", color = "black") +
  
  # Use annotate() to add a vertical line without triggering warnings
  annotate("segment", 
           x = consistency_data, xend = consistency_data,
           y = y_start, yend = y_end,
           color = "orange", linewidth = 1.5) +
  
  labs(x = "Consistency of random sub-group allocations", y = "Count") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.margin = margin(t = 20, r = 30, b = 40, l = 40)
  )
gg

ggsave(filename = paste0(plot_dir, 'consistency_subgrouping_', file, '.png'), plot = gg, width = 8, height = 5, dpi = 300)

}




#plot the distance groups traveled against duration apart

#making new columns to add duration and max dist apart info for plotting

split_merge_df$dyad_dist_max_1.2 <- NA
split_merge_df$dyad_dist_max_1.3 <- NA
split_merge_df$dyad_dist_max_2.3 <- NA

split_merge_df$event_dur <- NA

###
###-----This needs work - need to calculate duration of events on a dyadic level as this current code gives inaccurate values for distance apart against duration 
###
###

# Loop through each row of the dataframe
for (i in 1:(nrow(split_merge_df) - 1)) {
  
  # Print debug information
  #cat("Processing row:", i, "\n")
  
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
      #cat("Invalid time range at row:", i, "\n")
    }
  }
}


#clean the environment
#rm(list = ls()[grep("group_", ls())])
rm(group_1_xs, group_1_ys, group_2_xs, group_2_ys, group_3_xs, group_3_ys, split_i, merge_i, assoc_df, assoc_matrix, event, merges_df_subset, splits_df_subset, sub_sizes, subgroups, time_diff, time_diff_list, Rs, orig_sizes, merge_then_split, split_then_merge, split_merge_change_order, event_all, overall_event_association)

#Dyadic distance ~ event duration

#for plotting need to reshape to long format
long_df <- split_merge_df %>%
  pivot_longer(
    cols = starts_with("dyad_dist_max_"),
    names_to = "dyad_pair",
    values_to = "dyad_distance"
  )
#removing -Inf values
long_df_clean <- long_df[is.finite(long_df$dyad_distance), ]

#convert duration to hours
long_df_clean$event_dur_hours <- long_df_clean$event_dur / 60

if(make_plots == T){

p4 <- ggplot(long_df_clean, aes(x = dyad_distance, y = event_dur_hours)) +
  geom_point(color = "slateblue4", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = "All durations", x = "Dyad distance between subgroups (meters)", y = NULL) +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )


#looking at the duration and max dyadic distance for events within 1 day
one_day_df <- long_df_clean[long_df_clean$event_dur_hours <= 24, ]

p5 <- ggplot(one_day_df, aes(x = dyad_distance, y = event_dur_hours)) +
  geom_point(color = "slateblue1", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = " < 24 hours", x = NULL, y = NULL) +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )

#this plot suggests a subgroup were 2km apart for 5 hours so would have needed to travel on average 400m an hour to reach them - I assume this is not GPS error

#looking at event durations over 24 hours long
over_one_day_df <- long_df_clean[long_df_clean$event_dur_hours > 24, ]

p6 <- ggplot(over_one_day_df, aes(x = dyad_distance, y = event_dur_hours)) +
  geom_point(color = "slateblue1", size = 2, alpha = 0.6) +  # Colored points by dyad
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # One overall regression line
  labs(title = " > 24 hours", x = NULL, y = NULL) +
  theme_classic()+
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # top, right, bottom, left (in 'pt')
  )


combined_plot <- (p5 + p6) / p4 + plot_layout(heights = c(1, 1))  # Equal heights

png(file = paste0(plot_dir, "/duration_dyadicdist_2yr_", file, ".png"), width = 10, height = 6, units = "in", res = 300)
    
grid.newpage()
grid.draw(patchworkGrob(combined_plot))
grid.text("Event Duration (hours)", x = unit(0.02, "npc"), rot = 90, gp = gpar(fontsize = 14))

dev.off()

}




#which individuals often leave/merge on their own?

if(include_singletons == T){

  event_type <- "merge"
  event_all <- split_merge_df[split_merge_df$event == event_type,]  
  
  # Set up count vectors
  singleton_count <- setNames(rep(0, n_inds), as.character(1:n_inds))
  subgroup_count <-  setNames(rep(0, n_inds), as.character(1:n_inds))
  
  for (i in 1:nrow(event_all)) {
    row_i <- event_all[i, ]
    
    # Extract all subgroups
    subgroups <- list(row_i$sub1[[1]], row_i$sub2[[1]], row_i$sub3[[1]], 
                      row_i$sub4[[1]], row_i$sub5[[1]])
    
    for (subgroup in subgroups) {
      if (!is.null(subgroup) && !all(is.na(subgroup))) {
        ids <- as.character(subgroup)
        
        if (length(ids) == 1) {
          # Singleton individual
          singleton_count[ids] <- singleton_count[ids] + 1
        } else if (length(ids) > 1) {
          # Group split
          subgroup_count[ids] <- subgroup_count[ids] + 1
        }
      }
    }
  }
}

if(include_singletons == T && make_plots == T){
  
  singleton_count <- data.frame(id = lion_ids$tag_id, count = singleton_count)
  subgroup_count <- data.frame(id = lion_ids$tag_id, count = subgroup_count)
  
  # Rename columns for clarity
  colnames(singleton_count)[2] <- "singleton"
  colnames(subgroup_count)[2] <- "group"
  
  # Combine into one data frame
  combined_df <- merge(singleton_count, subgroup_count, by = "id")
  
  # Reshape to long format for stacking
  long_df <- pivot_longer(combined_df, cols = c("singleton", "group"),
                          names_to = "event_type", values_to = "count")
  proportions_df <- long_df %>%
    pivot_wider(names_from = event_type, values_from = count, values_fill = 0) %>%
    mutate(
      total = singleton + group,
      proportion_singleton = singleton / total
    )
  
  # Plot
stacked_plot <- ggplot(long_df, aes(x = id, y = count, fill = event_type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("singleton" = "slateblue3", "group" = "slateblue1")) +
    theme_classic() +
    labs(title = paste0("Number of ", event_type, " events per individual"),
         x = NULL, y = "Count", fill = "Event Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ 
  geom_text(
      data = proportions_df,
      aes(label = paste0(round(proportion_singleton * 100, 1), "%"), y = total + 60),  # adjust `+ 50` as needed
      color = "black", size = 3)

ggsave(filename = paste0(plot_dir, 'singleton_and_subgroup_count_perind_per', event_type, '.png'), plot = stacked_plot, width = 5, height = 4, dpi = 300)



}


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#-----DIADIC DISTANCES---------------------------------

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

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

# Calculate dyadic distances
dyadic_distances <- calculate_dyadic_distances(utm_long)

# Generate distance data for all individuals
distance_data_list <- lapply(1:nrow(xs), generate_distance_df)

# Combine all distance data into one data frame
combined_distance_data <- bind_rows(distance_data_list)

# Replace ID numbers with actual names
combined_distance_data$Focal_Individual <- lion_ids$name[combined_distance_data$Focal_Individual]
combined_distance_data$Other_Individual <- lion_ids$name[combined_distance_data$Other_Individual]

# remove the dyadic data for Simba after his death

combined_distance_data <- combined_distance_data[
  !(combined_distance_data$Time > as.character(ts[last_valid_index]) &
      (combined_distance_data$Focal_Individual == "Simba" |
         combined_distance_data$Other_Individual == "Simba")),
]

# Fill missing data using linear interpolation and carry forward/backward
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

# Example usage
individual_id <- "Sarabi"
#plot_individual_distances(individual_id, combined_distance_data)

if(make_plots == T){

# Loop through each individual and generate plots
individual_ids <- unique(combined_distance_data$Focal_Individual)
for (individual_id in individual_ids) {
  plot <- plot_individual_distances(individual_id, combined_distance_data) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0(plot_dir, "Dyadic_Distances_2yr_", individual_id, ".png"), plot = plot, width = 9, height = 5, dpi = 300)
  }


#add these to doc - all data, do for 0-100, log of overall
#add speed color to plot - are they resting when 10m apart?
ggplot(combined_distance_data, aes(Distance))+
  geom_histogram()+
  facet_wrap(~Focal_Individual)+
  xlim(0, 1000)

hist(combined_distance_data$Distance, breaks = 100000, xlim = c(0, 100))

}

rm(distance_data_list, dyadic_distances, long_df, long_df_clean, utm_long, xs_long, ys_df, xs_df, ys_long)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------


#------------SUBGROUP SIZE ~ TIME OF DAY-----------------------------


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#only done for all the data (with singleton events)
if(make_plots == T){

hist(subgroup_data$n_subgroups)

n_subs <- data.frame(n_subgroups = subgroup_data$n_subgroups, mean_subsize = colMeans(subgroup_data$subgroup_counts, na.rm = T), datetime = ts)

n_subs$time_hms <- as_hms(n_subs$datetime)

p0 <- ggplot(n_subs, aes(x = time_hms, y = n_subgroups)) +
  geom_smooth(se = FALSE, method = "loess", color = "steelblue2", size = 1.2) +
  theme_classic() +
  xlab("Time of Day") +
  ylab("Number of Subgroups") +
  labs(title = "Smoothed Number of Subgroups Throughout the Day")

p1 <- ggplot(n_subs, aes(x = time_hms, y = mean_subsize)) +
  geom_smooth(se = FALSE, method = "loess", color = "steelblue3", size = 1.2) +
  theme_classic() +
  xlab("Time of Day") +
  ylab("Subgroups Size") +
  labs(title = "Smoothed Subgroup Size Throughout the Day")


n_subs$time_bin <- floor_date(n_subs$datetime, unit = "1 hour") %>% format("%H:%M")

# Summarize mean and standard error
agg_df <- n_subs %>%
  group_by(time_bin) %>%
  summarise(
    mean_n_subs = mean(n_subgroups, na.rm = TRUE),
    mean_size = mean(mean_subsize, na.rm = TRUE),
    se_n_subs = sd(n_subgroups, na.rm = TRUE) / sqrt(n()),
    se_mean_size = sd(mean_subsize, na.rm = TRUE) / sqrt(n())
  )

# Plot with error bars
p2 <- ggplot(agg_df, aes(x = time_bin, y = mean_n_subs)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_n_subs - se_n_subs,
                    ymax = mean_n_subs + se_n_subs),
                width = 0.3, color = "black") +
  theme_classic() +
  xlab("Time of Day (UTC)") +
  ylab("Mean Number of Subgroups") +
  labs(title = "Mean Number of Subgroups per Hour with SE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

p3 <- ggplot(agg_df, aes(x = time_bin, y = mean_size)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_size - se_mean_size,
                    ymax = mean_size + se_mean_size),
                width = 0.3, color = "black") +
  theme_classic() +
  xlab("Time of Day (UTC)") +
  ylab("Mean Subgroup Size") +
  labs(title = "Mean Subgroup Size per Hour with SE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


gg <- (p0 + p1) / (p2 + p3)
gg

ggsave(filename = paste0(plot_dir, "subgroupsize_time.png"), plot = gg, width = 10, height = 8, dpi = 300)


}


rm(agg_df, event_all, n_subs, p0, p1, p2, p3, gg, row_i)



#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#Each individuals proportion time alone vs any other situation - each month
#Proportion of time group is together vs any other situation - each month

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

test <- matrix(ncol = n_times, nrow = n_inds)
for (i in 1:ncol(subgroup_data$ind_subgroup_membership)){
  
  for (j in 1:n_inds){
    
    group_i_j <- subgroup_data$ind_subgroup_membership[j,i]
    group_i <- subgroup_data$ind_subgroup_membership[,i]
    sub_counts <- table(group_i)
    
    
    if(group_i_j %in% group_i[-j]){
      test[j,i] <- 0
      
      
    } else{ 
      test[j,i] <- 1
      }
  }
}





#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#Investigate how group dynamics shift before and after Simba's death

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#split subgrouping data to before and after Simba died first to see whether this affects the results 

bef <- subgroup_data$subgroup_counts[,1:last_valid_index]
aft <- subgroup_data$subgroup_counts[,last_valid_index:ncol(subgroup_data$subgroup_counts)]

bef_n_groups <- subgroup_data$n_subgroups[1:last_valid_index]
aft_n_groups <- subgroup_data$n_subgroups[last_valid_index:ncol(subgroup_data$subgroup_counts)]

# Convert to long format with group labels
bef_long <- data.frame(value = as.vector(as.matrix(bef)), period = "Before")
aft_long <- data.frame(value = as.vector(as.matrix(aft)), period = "After")
bef_long_n_groups <- data.frame(value = as.vector(as.matrix(bef_n_groups)), period = "Before")
aft_long_n_groups <- data.frame(value = as.vector(as.matrix(aft_n_groups)), period = "After")

# Combine both
combined_long <- bind_rows(bef_long, aft_long)
combined_long_n_groups <- bind_rows(bef_long_n_groups, aft_long_n_groups)

combined_long$period <- factor(combined_long$period, levels = c("Before", "After"))
combined_long_n_groups$period <- factor(combined_long_n_groups$period, levels = c("Before", "After"))

if(make_plots == T){

# Plot with ggplot2
p1 <- ggplot(combined_long, aes(x = value)) +
  geom_histogram(bins = 6, fill = "slateblue1", color = "slateblue3", alpha = 0.5) +
  facet_wrap(~ period, ncol = 2) +
  ylim(0, 20000) +
  labs(x = "Subgroup Count", y = "Frequency", title = "Distribution of subgroup sizes before and after Simba's death") +
  theme_classic()


# Plot with ggplot2
p2 <- ggplot(combined_long_n_groups, aes(x = value)) +
  geom_histogram(bins = 6, fill = "gold", color = "gold4", alpha = 0.7) +
  facet_wrap(~ period, ncol = 2) +
  ylim(0, 12000) +
  labs(x = "Number of subgroups", y = "Frequency", title = "Distribution of the number of subgroups before and after Simba's death") +
  theme_classic()

gg <- p1 / p2

ggsave(filename = paste0(plot_dir, 'subgrp_dynamics_bef_aft_simba_2yr.png'), plot = gg, width = 7, height = 6, dpi = 300)
}

rm(aft, aft_long, aft_long_n_groups, bef, bef_long, bef_long_n_groups, combined_long_n_groups, combined_long, gg, p1, p2)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#--------------PLOT DYADIC DISTANCES IN NETWORK------------------------

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

combined_distance_data$date <- as.POSIXct(combined_distance_data$Time, 
                                              format = "%Y-%m-%d", tz = "UTC")

#get the mean, median, and sd for diadic distances
dyad_df <- compute_dyadic_metrics(combined_distance_data)

# Convert the 3D array to long format
dyad_df <- melt(dyad_df, varnames = c("Focal", "Other", "Statistic"), value.name = "Distance")

# Now filter to just the 'mean' distances and exclude self-pairs
dyad_edges <- dyad_df %>%
  filter(Statistic == "mean", Focal != Other)

# Convert distance to weight (inverse for layout)
dyad_edges$weight <- 1 / dyad_edges$Distance

# Create igraph object with weights
g <- graph_from_data_frame(dyad_edges, directed = FALSE, vertices = lion_ids)

# Use Kamada-Kawai layout (preserves edge length relative to weights)
layout_kk <- layout_with_kk(g, weights = E(g)$weight)

if(make_plots == T){
# Plot
gg <- ggraph(g, layout = "manual", x = layout_kk[,1], y = layout_kk[,2]) +
  geom_edge_link(aes(edge_width = weight), color = "lightgrey", show.legend = FALSE) +
  geom_node_point(aes(color = sex, shape = age), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Male" = "cyan4", "Female" = "orange1")) +
  scale_shape_manual(values = c("Adult" = 18, "Sub-adult" = 20)) +
  scale_edge_width_continuous(range = c(0.2, 2))+
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Mean Dyadic Distance Network")

ggsave(filename = paste0(plot_dir, 'dyadic_distance_all.png'), plot = gg, width = 6, height = 5, dpi = 300)

}

#compare before and after Simba died:

#split the data into before and after SImba's date of death
combined_distance_data_list <- list(combined_distance_data[combined_distance_data$date < as.Date(ts[last_valid_index]),], combined_distance_data[combined_distance_data$date >= as.Date(ts[last_valid_index]),])

dyad_before <- compute_dyadic_metrics(combined_distance_data_list[[1]])
dyad_after  <- compute_dyadic_metrics(combined_distance_data_list[[2]])

df_before <- melt(dyad_before, varnames = c("Focal", "Other", "Statistic"), value.name = "Distance")
df_before$Period <- "Before"
df_after <- melt(dyad_after, varnames = c("Focal", "Other", "Statistic"), value.name = "Distance")
df_after$Period <- "After"
#remove SImba to calculations for after
df_after$Distance[df_after$Focal == "Simba" | df_after$Other == "Simba"] <- NA

dyad_all <- rbind(df_before, df_after)

dyad_all_wide <- dyad_all %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Distance
  )

dyad_all_wide <- subset(dyad_all_wide, Focal != Other)
dyad_all_wide$Period <- factor(dyad_all_wide$Period, levels = c("Before", "After"))

if(make_plots == T){

plot <- ggplot(dyad_all_wide, aes(x = interaction(Focal, Other), y = mean, fill = Period)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Before" = "slateblue4", "After" = "slateblue1")) +
  geom_errorbar(
    aes(ymin = pmax(mean - sd, 0),
      ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.7, color = "lightgrey"
  ) +
  labs(x = "Dyad", y = "Mean Distance (± SD)", title = "Dyadic Distances Before and After") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(filename = paste0(plot_dir, 'dyadic_distance_bef_aft_barplot.png'), plot = plot, width = 9, height = 4, dpi = 300)

}


#compare network before and after Simba's death

# Convert distance to weight (inverse for layout)
dyad_all_wide$weight <- 1 / dyad_all_wide$mean  


# Combine before and after (excluding Simba from "after")
dyad_avg <- dyad_all_wide %>%
  filter(Focal != Other) %>%
  group_by(Focal, Other) %>%
  summarise(mean_weight = mean(1 / mean, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.infinite(mean_weight), !is.na(mean_weight)) %>%
  select(from = Focal, to = Other, weight = mean_weight)

# Create a layout graph
g_layout <- graph_from_data_frame(dyad_avg, directed = FALSE, vertices = lion_ids)

# Use KK layout for consistency
layout_fixed <- layout_with_kk(g_layout, weights = E(g_layout)$weight)

dyad_all_wide$weight <- 1 / dyad_all_wide$mean  

df_before_edges <- dyad_all_wide %>%
  filter(Period == "Before", Focal != Other) %>%
  select(from = Focal, to = Other, weight = weight)

df_after_edges <- dyad_all_wide %>%
  filter(Period == "After", Focal != Other) %>%
  select(from = Focal, to = Other, weight = mean)
#remove Simba
df_after_edges <- df_after_edges[!(df_after_edges$from == "Simba") & !(df_after_edges$to == "Simba"),]

#Create graph objects
g_before <- graph_from_data_frame(df_before_edges, directed = FALSE, vertices = lion_ids)
g_after <- graph_from_data_frame(df_after_edges, directed = FALSE, vertices = lion_ids[-3,])

# Reuse same layout coordinates for before and after
V(g_before)$x <- layout_fixed[, 1]
V(g_before)$y <- layout_fixed[, 2]

V(g_after)$x <- layout_fixed[, 1][match(V(g_after)$name, V(g_before)$name)]
V(g_after)$y <- layout_fixed[, 2][match(V(g_after)$name, V(g_before)$name)]

if(make_plots == T){

p_before <- ggraph(g_before, layout = "manual", x = V(g_before)$x, y = V(g_before)$y) +
  geom_edge_link(aes(edge_width = weight), color = "lightgrey", show.legend = FALSE) +
  geom_node_point(aes(color = sex, shape = age), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Male" = "cyan4", "Female" = "orange1")) +
  scale_shape_manual(values = c("Adult" = 18, "Sub-adult" = 20)) +
  scale_edge_width_continuous(range = c(0.2, 2)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "left") +
  labs(title = "Mean dyadic distance before Simba")

p_after <- ggraph(g_after, layout = "manual", x = V(g_after)$x, y = V(g_after)$y) +
  geom_edge_link(aes(edge_width = weight), color = "lightgrey", show.legend = FALSE) +
  geom_node_point(aes(color = sex, shape = age), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Male" = "cyan4", "Female" = "orange1")) +
  scale_shape_manual(values = c("Adult" = 18, "Sub-adult" = 20)) +
  scale_edge_width_continuous(range = c(0.2, 2)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "none") +
  labs(title = "Mean dyadic distance after Simba")


gg <- p_before + p_after
gg

ggsave(filename = paste0(plot_dir, 'dyadic_distance_bef_aft.png'), plot = gg, width = 9, height = 5, dpi = 300)

}




#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

#-------------------------ANIMATION TIME----------------------------

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------


#now want to make the network into an animation over time for before Simba died (as can't get it to work when he's not in the group)

# Step 1: Prepare daily dyadic distances
daily_edges <- combined_distance_data %>%
  filter(Focal_Individual != Other_Individual) %>%
  group_by(date, Focal_Individual, Other_Individual) %>%
  summarise(mean_distance = mean(Distance, na.rm = TRUE), .groups = "drop") %>%
  rename(from = Focal_Individual, to = Other_Individual, weight = mean_distance)

# Optional: reorder columns for clarity
daily_edges <- daily_edges[, c("from", "to", "weight", "date")]

# Step 2: Prepare node (vertex) data
all_nodes <- unique(c(daily_edges$from, daily_edges$to))

nodes <- lion_ids %>%
  filter(name %in% all_nodes) %>%
  rename(name = name)  # ensures consistency with 'from'/'to' columns

# Step 3: Create a static layout (fixed across all days)
g_base <- graph_from_data_frame(
  daily_edges %>% filter(date == min(date)),
  vertices = nodes,
  directed = FALSE
)

layout_static <- create_layout(g_base, layout = "kk")  # Kamada-Kawai layout

# Step 4: Save layout as data frame for merging
layout_static_df <- layout_static %>%
  as.data.frame() %>%
  select(name, x, y)  # x/y = static node positions

# Step 5: Add coordinates to edges (for animation)
edges_for_plot <- daily_edges %>%
  left_join(layout_static_df, by = c("from" = "name"), suffix = c("", ".from")) %>%
  left_join(layout_static_df, by = c("to" = "name"), suffix = c("", ".to")) %>%
  rename(
    xend = x.to,
    yend = y.to
  )

network_plot <- ggplot() +
  # Edges: from `edges_for_plot`, using x/xend/y/yend
  geom_segment(data = edges_for_plot,
               aes(x = x, y = y, xend = xend, yend = yend, 
                   size = 1 / weight, alpha = 1 / weight),
               color = "grey40") +
  
  # Nodes: from layout_static_df (no xend/yend here)
  geom_point(data = layout_static_df, aes(x = x, y = y, color = name), size = 6) +
  geom_text(data = layout_static_df, aes(x = x, y = y, label = name), vjust = 1.5, size = 3) +
  scale_size(range = c(0.2, 5)) +
  scale_alpha(range = c(0.2, 1)) +
  theme_void() +
  labs(title = "Dyadic Distance Network - {closest_state}") +
  transition_states(edges_for_plot$date, transition_length = 2, state_length = 1) +
  ease_aes("linear")



# Step 7: Render or save animation
animate(network_plot, width = 800, height = 600, fps = 10, duration = 20, renderer = gifski_renderer("../results/level0/dyadic_network.gif"))

#-------------------------------------------------------------------

#changing the node locations per day

#-------------------------------------------------------------------


all_nodes <- unique(c(daily_edges$from, daily_edges$to))

nodes <- lion_ids %>%
  filter(name %in% all_nodes) %>%
  distinct(name, .keep_all = TRUE)

# List of layouts per day
layouts_by_day <- daily_edges %>%
  group_by(date) %>%
  group_split() %>%
  map_dfr(function(day_edges) {
    g <- graph_from_data_frame(day_edges, directed = FALSE, vertices = nodes)
    
    # Convert distances to similarity (shorter distance = stronger tie)
    E(g)$weight <- 1 / day_edges$weight
    
    layout <- layout_with_kk(g, weights = E(g)$weight)
    layout_df <- as_tibble(layout)
    layout_df$name <- V(g)$name
    layout_df$date <- unique(day_edges$date)
    
    return(layout_df)
  })

# Merge layout with edge data per day
edges_for_plot <- daily_edges %>%
  left_join(layouts_by_day, by = c("from" = "name", "date")) %>%
  rename(x = V1, y = V2) %>%
  left_join(layouts_by_day, by = c("to" = "name", "date")) %>%
  rename(xend = V1, yend = V2)


layouts_by_day <- layouts_by_day %>%
  left_join(lion_ids, by = "name")


#removing Simba
layouts_by_day_clean <- layouts_by_day[!(layouts_by_day$name == "Simba" & layouts_by_day$date > as.Date(ts[last_valid_index])), ]
edges_for_plot_clean <- edges_for_plot[!((edges_for_plot$from == "Simba" | edges_for_plot$to == "Simba") & edges_for_plot$date > as.Date(ts[last_valid_index])), ]


network_plot <- ggplot() +
  # Edge layer
  geom_segment(
    data = edges_for_plot_clean,
    aes(x = x, y = y, xend = xend, yend = yend, size = 1 / weight, alpha = 1 / weight),
    color = "grey40", show.legend = FALSE
  ) +
  # Node points with sex and age mappings
  geom_point(
    data = layouts_by_day_clean,
    aes(x = V1, y = V2, color = sex, shape = age),
    size = 5
  ) +
  # Node labels
  geom_text(
    data = layouts_by_day_clean,
    aes(x = V1, y = V2, label = name),
    vjust = 1.5,
    size = 3
  ) +
  scale_color_manual(values = c("Male" = "cyan4", "Female" = "orange1")) +
  scale_shape_manual(values = c("Adult" = 18, "Sub-adult" = 20)) +
  scale_size(range = c(0.2, 2)) +
  scale_alpha(range = c(0.2, 1)) +
  theme_void() +
  labs(title = 'Mean Dyadic Distance Network — {frame_time}') +
  transition_time(date) +
  ease_aes("linear")

animate(network_plot, width = 8, height = 6, units = "in", fps = 30, res = 300, duration = 120, renderer = av_renderer("../results/level0/dyadic_network_clean.mp4"))


