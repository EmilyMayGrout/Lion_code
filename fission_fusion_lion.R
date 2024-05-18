#Questions interested for the lion dataset:
#1. duration of splits
#2. who do they split with
#3. what are they doing when split?
#4. how far do they travel apart?
#5. how often do they come together?


#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/lion/data/processed/"
code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
coati_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/lion/results/level0/'
gps_file <- "lion_xy_15min_level0.RData" #level0 is when Venus is not removed
id_file <- 'lion_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
Rs <- c(100,200,300,400,500,1000)
Rs <- c(1000,2000,3000,4000,5000,10000)
Rs <- c(10000,20000,30000,40000,50000,100000)


R <- 50

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
source('coati_function_library_V1.R')

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

R = 10
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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_10m.png'))

visualize_network_matrix_galaxy(ffnet_reorder, lion_ids[new_order,])
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
R = 1000
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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'within_1000m_subgroup_network_3m.png'))
visualize_network_matrix_galaxy(within_group_data$proximity_net, lion_ids[new_order,])
dev.off()

#--------------------------------------------------------

#now want to work out how long they were apart for...

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

#add the time of each event
split_merge_df$datetime <- ts[split_merge_df$t]
split_merge_df$round_time <- round_date(split_merge_df$datetime, "1 hour")

split_merge_df$hour <- strftime(split_merge_df$round_time, "%H")
split_merge_df$hour <- as.numeric(split_merge_df$hour)


p1 <- ggplot(split_merge_df, aes(x=hour, fill=event)) +
  geom_histogram(alpha=0.5, position="dodge", aes(y=..density..), bins=24)+
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


