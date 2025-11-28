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
#got the subsampled_2h_4h matrix and t_sub which are the times matching the matrix in fission_fusion_across_sampling_rates
load("../data/processed/subsampled_2h4h_matrix.RData")

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


#----------------------------------------------------------

#running the consistency metric for the 2h/4h to compare with 15 minute sampling
#get indices for 2h/4h 
indices_in_ts <- match(t_sub, ts)
indices_in_ts <- indices_in_ts[!is.na(indices_in_ts)]

sample_rate <- "15min"

if(sample_rate == "15min"){
  
  xs <- xs
  ys <- ys
  
  }else if(sample_rate == "2h4h"){
    
  xs <- xs[,indices_in_ts]
  ys <- ys[,indices_in_ts]  
  ts <- t_sub
  
  }


R = 60

subgroup_data <- get_subgroup_data(xs, ys, R)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_', R, 'm_', sample_rate, '.png'))

visualize_lion_network(ffnet_reorder, lion_ids[new_order,])

dev.off()
  
#-------------------------------------------------------------------

#now compare on a monthly scale:
R = 600
subgroup_data <- get_subgroup_data(xs, ys, R)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
split_by_month <- F
year_months <- unique(format(ts, "%Y-%m"))
ff_net_month <- array(NA, dim = c(n_inds, n_inds, length(year_months)))
dimnames(ff_net_month)[[3]] <- year_months
  
#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
      
      if(split_by_month){
        for(m in 1:length(seq_along(year_months))){
          
          ym <- year_months[m]
          # Logical index for matching year-month
          ts_filter <- format(ts, "%Y-%m") == ym
          
          sub_ids_i <- subgroup_data$ind_subgroup_membership[i, ts_filter]
          sub_ids_j <- subgroup_data$ind_subgroup_membership[j, ts_filter]
          
          ff_net_month[i, j, m] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
          #print(ff_net_month[i, j, m])
        }}
      
      #getting subgroup id for individual i and j
      sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
      sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
      
      #computing edge weight (fraction of time in same subgroup)
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
    }}
  
for (i in 1:dim(ff_net_month)[3]) {
    diag(ff_net_month[,,i]) <- NA
  }
    
new_order <- c(1:3,5,4,6)
ffnet_reorder <- ff_net[new_order, new_order]
    
png(height = 2200, width = 3200, units = 'px', filename = paste0(plot_dir,'monthly_subgroup_network_', R, 'm_', sample_rate, '.png'))
par(mfrow=c(4,5)) 
    
for(i in 1:dim(ff_net_month)[3]){
  ffnet_reorder <- ff_net_month[new_order, new_order, i]
  visualize_lion_network(ffnet_reorder, lion_ids[new_order,])
  mtext(line=1, side=3, year_months[i] , outer=F)
    }
    
dev.off()
  



#--------------------------------------------------------------

run_consistency <- F
if(run_consistency == T){

#function code in lion_functions
splits_df <- detect_splits(subgroup_data$ind_subgroup_membership)

#function that takes in a splits dataframe and outputs a p_dyad_together matrix (probability of splitting together for each dyad)
p_dyad_together <- get_p_dyad_together(splits_df_local = splits_df, n_inds_local = n_inds)

#Compute consistency based on consistency metric (see coati_function_library)
consistency_data <- get_consistency(p_dyad_together)

#randomize splits and recompute consistency metric
set.seed(24)
n_rands <- 1000
rando_consistencies <- rep(NA, n_rands)
for (i in 1:n_rands){
  
  rando <- randomise_splits(splits_df)
  rando_dyad <- get_p_dyad_together(splits_df_local = rando, n_inds_local = n_inds)
  consist <- get_consistency(rando_dyad)
  rando_consistencies[i] <- consist
  
}


# Create data frame for the histogram
df <- data.frame(rando_consistencies = rando_consistencies)

y_start <- 0
y_end <- 250
breaks = seq(0.05, 0.25, 0.002)

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

ggsave(filename = paste0(plot_dir, 'consistency_subgrouping_', sample_rate, "_", R, 'm.png'), plot = gg, width = 8, height = 5, dpi = 300)


# Detect events
splits <- detect_splits(subgroup_data$ind_subgroup_membership)
merges <- detect_merges(subgroup_data$ind_subgroup_membership)
colnames(merges)[c(2,8)] <- c("orig_group", "n_orig")
splits$event <- "split"
merges$event <- "merge"

split_merge_df <- rbind(splits, merges)
colnames(split_merge_df) <- c("t", "orig_group", "sub1", "sub2", "sub3", "sub4", "sub5", "n_orig", "n_sub1", "n_sub2", "n_sub3","n_sub4","n_sub5", "event")

#order the dataframe by time
split_merge_df <- split_merge_df[(order(split_merge_df$t)),]

#remove duplicates
split_merge_df <- split_merge_df[!duplicated(split_merge_df[,1:5]),]

#change index of the rows to time order
rownames(split_merge_df) <- 1:nrow(split_merge_df)

#add the time of each event
split_merge_df$datetime <- t_sub[split_merge_df$t]
split_merge_df$round_time <- round_date(split_merge_df$datetime, "1 hour")

split_merge_df$hour <- strftime(split_merge_df$round_time, "%H")
split_merge_df$hour <- as.numeric(split_merge_df$hour)

gg <- ggplot(split_merge_df, aes(x=hour, fill=event)) +
  scale_fill_manual(values = c("turquoise1", "dodgerblue4"))+
  geom_density(alpha = 0.5)+ 
  scale_color_brewer(palette="Accent")+ 
  theme_classic()
gg

ggsave(filename = paste0(plot_dir, 'splitmerge_overtime_', sample_rate, "_", R, 'm.png'), plot = gg, width = 8, height = 5, dpi = 300)



if(sample_rate == "2h4h"){
  split_merge_df_2h4h <- split_merge_df
}else if(sample_rate == "15min"){
  split_merge_df_15min <- split_merge_df
}



splits_15min <- split_merge_df_15min[split_merge_df_15min$event == "split",]
splits_2h4h <- split_merge_df_2h4h[split_merge_df_2h4h$event == "split",]
splits_15min$sampling <- "15min"
splits_2h4h$sampling <- "2h4h"
combined_splits <- rbind(splits_15min, splits_2h4h)

merges_15min <- split_merge_df_15min[split_merge_df_15min$event == "merge",]
merges_2h4h <- split_merge_df_2h4h[split_merge_df_2h4h$event == "merge",]
merges_15min$sampling <- "15min"
merges_2h4h$sampling <- "2h4h"
combined_merges <- rbind(merges_15min, merges_2h4h)

g1 <- ggplot(combined_splits, aes(x = hour, fill = sampling)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("15min" = "red", "2h4h" = "turquoise")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Split events") 

g2 <- ggplot(combined_merges, aes(x = hour, fill = sampling)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("15min" = "red", "2h4h" = "turquoise")) +
  theme_classic() +
  labs(title = "Merge events", fill = "Sampling Interval") 

gg <- g1 + g2
gg

ggsave(filename = paste0(plot_dir, 'splitmerge_overtime_', R, 'm.png'), plot = gg, width = 10, height = 5, dpi = 300)

}


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

#2h/4h vs 15 min prop time each is alone

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


n_inds <- nrow(xs)
n_times <- ncol(xs)

singles <- matrix(ncol = n_times, nrow = n_inds)

for (i in 1:ncol(subgroup_data$ind_subgroup_membership)) {
  for (j in 1:n_inds) {
    
    group_i_j <- subgroup_data$ind_subgroup_membership[j, i]
    group_i <- subgroup_data$ind_subgroup_membership[, i]
    
    if (is.na(group_i_j)) {
      singles[j, i] <- NA
    } else {
      # Count how many individuals are in the same subgroup
      same_group_members <- sum(group_i == group_i_j, na.rm = TRUE)
      
      # If only this individual is in that group, mark as single
      if (same_group_members == 1) {
        singles[j, i] <- 1
      } else {
        singles[j, i] <- 0
      }
    }
  }
}


# Get month index for each time point
month_indx <- format(ts, "%m")
months <- sort(unique(month_indx))  

# Initialize matrices
alone_count <- matrix(NA, nrow = n_inds, ncol = length(unique(month_indx)))
alone_prop <- matrix(NA, nrow = n_inds, ncol = length(unique(month_indx)))

# Column names as months
months <- sort(unique(month_indx))
colnames(alone_count) <- months
colnames(alone_prop) <- months

# Loop through each month
for (i in seq_along(months)) {
  ts_month_i <- which(month_indx == months[i])
  #print(ts[min(ts_month_i)])
  #print(ts[max(ts_month_i)])
  
  singles_month_i <- singles[, ts_month_i, drop = FALSE]
  
  # Count how often individual is alone (value == 1)
  alone_count[, i] <- rowSums(singles_month_i == 1, na.rm = TRUE)
  
  # Count number of valid time points (not NA)
  valid_counts <- rowSums(!is.na(singles_month_i))
  
  # Proportion of time alone = alone / valid observations
  alone_prop[, i] <- alone_count[, i] / valid_counts
}

# Add individual IDs as a column for reshaping
alone_count_df <- as.data.frame(alone_count)
alone_count_df$id <- 1:nrow(alone_count_df)
alone_prop_df <- as.data.frame(alone_prop)
alone_prop_df$id <- 1:nrow(alone_prop_df)

alone_count_long <- alone_count_df %>%
  pivot_longer(cols = -id, names_to = "month", values_to = "alone_count")
alone_prop_long <- alone_prop_df %>%
  pivot_longer(cols = -id, names_to = "month", values_to = "alone_prop")

alone_long_combined <- left_join(alone_count_long, alone_prop_long, by = c("id", "month"))

alone_long_combined <- alone_long_combined %>%
  mutate(name = lion_ids$name[id])


gg <- ggplot(alone_long_combined, aes(x=month, y = alone_prop, colour = name))+
  geom_point(size = 3)+
  geom_line(aes(group = name), alpha = 0.5, linewidth = 2)+
  labs(
    x = "Month",
    y = "Proportion of time alone",
    colour = "Individual ID"
  )+
  theme_classic()
gg

ggsave(filename = paste0(plot_dir, 'prop_timealone_', sample_rate, '.png'), plot = gg, width = 8, height = 5, dpi = 300)




#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#Proportion of time group is together vs any other situation - each month
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------



together <- matrix(ncol = n_times, nrow = 1)

for (i in 1:ncol(subgroup_data$ind_subgroup_membership)) {
  
  #if all individuals are in the same subgroup and there are 1 or no inds who are NA, the full group is together (NA rule because of Simba)
  if(length(unique(na.omit(subgroup_data$ind_subgroup_membership[,i]))) == 1 & sum(is.na(subgroup_data$ind_subgroup_membership[,i])) <= 1) {
    together[i] <- 1
  } else{
    together[i] <- 0
    
  }
}

together_prop <- matrix(NA, nrow = 1, ncol = length(unique(month_indx)))
colnames(together_prop) <- months

# code copied from above
for (i in seq_along(months)) {
  ts_month_i <- which(month_indx == months[i])
  together_month_i <- together[ts_month_i, drop = FALSE]
  
  # Count how often individual is alone (value == 1)
  together_count <- sum(together_month_i == 1, na.rm = TRUE)
  
  # Count number of valid time points (not NA)
  valid_counts <- sum(!is.na(together_month_i))
  
  # Proportion of time together = together / valid observations
  together_prop[i] <- together_count / valid_counts
}


together_long <- as.data.frame(t(together_prop))
colnames(together_long) <- "prop_together"
together_long$month <- sort(unique(month_indx))  

gg <- ggplot(together_long, aes(x=month, y = prop_together, group = 1))+
  geom_point(size = 3)+
  geom_line()+
  labs(
    x = "Month",
    y = "Proportion of time together"
  )+
  theme_classic()

gg

ggsave(filename = paste0(plot_dir, 'prop_timetogether', sample_rate, '.png'), plot = gg, width = 8, height = 5, dpi = 300)

