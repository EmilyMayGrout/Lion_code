#looking at the subgrouping dynamics across all the prides and coalitions
#look at seasonal changes

library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)
library(fields)
library(patchwork)
library(hms)

code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
plot_dir <-  'C:/Users/egrout/Dropbox/lion/results/level0/'

setwd(code_dir)
source('coati_function_library_V1.R')
source('lion_functions.R')

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
load('allnamib_xy_2hour_level0.RData')
load('allnamib_latlon_2hour_level0.RData') 
load('allnamib_ids.RData')

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


#get dyadic distances between all individuals 

# Convert matrices to data frames for easier handling
xs_df <- as.data.frame(xs)
ys_df <- as.data.frame(ys)
colnames(xs_df) <- ts
colnames(ys_df) <- ts
# Add an 'Individual' column
xs_df$Individual <- all_ids$name
ys_df$Individual <- all_ids$name

# Convert to long format
xs_df <- xs_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "X_UTM") 
ys_df <- ys_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "Y_UTM") 

# Merge the long data frames
utm_long <- xs_df %>%
  inner_join(ys_df, by = c("Individual", "Time"))

# Calculate dyadic distances
dyadic_distances <- calculate_dyadic_distances(utm_long)

#get time in posixct format - need to add time when its midnight, otherwise gives NA
dyadic_distances$Time <- as.POSIXct(
  ifelse(nchar(dyadic_distances$Time) == 10,
         paste0(dyadic_distances$Time, " 00:00:00"),
         dyadic_distances$Time),
  format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
)

# Get unique individuals for consistent colors
unique_individuals <- unique(dyadic_distances$Individual.2)
# Generate colors for each unique individual
colors <- rainbow(length(unique_individuals))
# Create a named vector with colors for each individual
color_mapping <- setNames(colors, unique_individuals)


for(i in unique(dyadic_distances$Individual.1)){

ind_i <- i

dyadic_distances_i <- dyadic_distances[dyadic_distances$Individual.1 == ind_i,]

#remove NA's in distances
dyadic_distances_i <- dyadic_distances_i[!is.na(dyadic_distances_i$Distance),]
#hist(dyadic_distances_i$Distance)

#filter to individuals who have been within 600m of one another
dyadic_distances_i <- dyadic_distances_i %>%
  group_by(Individual.2) %>%
  filter(min(Distance, na.rm = TRUE) <= 600) %>%
  ungroup()
#hist(dyadic_distances_i$Distance)

#get rolling mean for plotting a smoother line
dyadic_distances_i <- dyadic_distances_i %>%
  group_by(Individual.1, Individual.2) %>%
  arrange(Time) %>%
  mutate(
    RollingMean = rollmean(Distance, k = 12, fill = NA, align = "right")
  ) %>%
  ungroup()


plot <- ggplot(dyadic_distances_i, aes(x = Time)) +
  geom_line(aes(y = Distance, color = Individual.2, group = Individual.2), linewidth = 0.8, alpha = 0.8) +
  geom_line(aes(y = RollingMean, color = Individual.2, group = Individual.2), size = 3, alpha = 0.3) +
  scale_color_manual(values = color_mapping) +
  labs(title = paste("Dyadic Distances for Individual", ind_i),
       x = "Time",
       y = "Distance (meters)",
       colour = "Other Individual") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom") +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%m-%d")
plot

ggsave(filename = paste0(plot_dir, "allnamib_dyadicdist_", ind_i, ".png"), plot = plot, width = 12, height = 5, dpi = 300)

}


#-----------------------------------------------------------
#make matrix plots for prop time together for dyads that at somepoint were within 600m of one another 

R = 600

mean_prox_ind_monthly <- list()

year_month <- format(ts, "%Y-%m")

for(ym in unique(year_month)){
  
  # Subset time indices
  ts_ym <- which(year_month == ym)
  xs_ym <- xs[, ts_ym]
  ys_ym <- ys[, ts_ym]
  
  # Keep individuals with any data
  valid_inds <- rowSums(!is.na(xs_ym) | !is.na(ys_ym)) > 0
  xs_ym <- xs_ym[valid_inds, , drop = FALSE]
  ys_ym <- ys_ym[valid_inds, , drop = FALSE]
  ids_ym <- all_ids[valid_inds, , drop = FALSE]
  # Skip if proximity network is NULL/NA
  proxim_data <- get_proximity_data(xs_ym, ys_ym, R)
  if (is.null(proxim_data$proximity_net) || all(is.na(proxim_data$proximity_net))) next

  # Compute mean proximity per individual (row-wise)
  mean_per_ind <- rowMeans(proxim_data$proximity_net, na.rm = TRUE)
  # Combine with individual IDs
  df_month <- data.frame(
    YearMonth = ym,
    Individual = ids_ym$name,
    MeanProximity = mean_per_ind
  )
  mean_prox_ind_monthly[[ym]] <- df_month
  
  
  #plot the proximity matrices for the individuals which did have data for that month
  png(height = 800, width = 850, units = 'px', filename = paste0(plot_dir,'within_600m_', ym, '.png')) 
  # Visualize network with consistent positions
  visualize_lion_network(proxim_data$proximity_net, ids_ym)
  mtext(ym, side = 3, line = 0, cex = 1.5, font = 2)
  dev.off()

}


# Combine into a single dataframe
mean_prox_ind_monthly_df <- do.call(rbind, mean_prox_ind_monthly)

#adding 01 to date, so its the first of each month
mean_prox_ind_monthly_df$YearMonth <- as.Date(paste0(mean_prox_ind_monthly_df$YearMonth, "-01"))

 mean_prox_ind_monthly_df[mean_prox_ind_monthly_df$MeanProximity == 0,] <- NA

color_mapping <- setNames(all_ids$color, all_ids$name)


ggplot(mean_prox_ind_monthly_df, aes(x = YearMonth, y = MeanProximity, color = Individual, group = Individual)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = color_mapping) + 
  scale_x_date(date_labels = "%y-%m", date_breaks = "1 year") +
  labs(
    title = "Mean Proximity to other individuals within 600m",
    x = "Year-Month",
    y = "Mean Proximity (m)",
    color = "Individual"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

#-------------------------------------------------------------------

#when inds are within 600m, what is the proportion of time they are within 60m


#don't run!! just load file saved below
# subgroup_data_600 <- get_subgroup_data(xs, ys, 600)
# subgroup_membership_600 <- subgroup_data_600$ind_subgroup_membership
# subgroup_data_60 <- get_subgroup_data(xs, ys, 60)
# subgroup_membership_60 <- subgroup_data_60$ind_subgroup_membership
# 
# ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
# split_by_month <- T
# year_months <- unique(format(ts, "%Y-%m"))
# ff_net_month <- array(NA, dim = c(n_inds, n_inds, length(year_months)))
# dimnames(ff_net_month)[[3]] <- year_months
# 
# for(i in 1:n_inds){
#   for(j in 1:n_inds){
#     
#     if(split_by_month){
#       for(m in 1:length(seq_along(year_months))){
#         
#         ym <- year_months[m]
#         
#         # Logical index for matching year-month
#         ts_filter <- format(ts, "%Y-%m") == ym
#         
#         sub_ids_i <- subgroup_data_60$ind_subgroup_membership[i, ts_filter]
#         sub_ids_j <- subgroup_data_60$ind_subgroup_membership[j, ts_filter]
#         
#         ff_net_month[i, j, m] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
#         #print(ff_net_month[i, j, m])
#       }
#     }
#     
#     #getting subgroup id for individual i and j
#     sub_ids_i <- subgroup_data_60$ind_subgroup_membership[i,]
#     sub_ids_j <- subgroup_data_60$ind_subgroup_membership[j,]
#     
#     #computing edge weight (fraction of time in same subgroup)
#     ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
#   }
# }
# 
# diag(ff_net) <- NA

#save(ff_net, ff_net_month, dyadic_distances, all_ids, subgroup_data_60, file = "allnamib_60m_ffnet_monthly.RData")

load("allnamib_60m_ffnet_monthly.RData")

#look at one of the months plots 
visualize_lion_network(ff_net_month[,,15], all_ids)


#-------------------------------------------------------
#plot subgroup sizes across the year
hist(subgroup_data_60$n_subgroups)

#make empty dataframe to store mean subsize for each datapoint 
mat <- subgroup_data_60$ind_subgroup_membership

mean_subs <- data.frame(
  time = ts,
  n_tracked = colSums(!is.na(mat)),
  n_subs = subgroup_data_60$n_subgroups
)

mean_subs$mean_subs <- mean_subs$n_subs/mean_subs$n_tracked

ggplot(data = mean_subs, aes(x= time, y = n_subs))+
  geom_line(aes(x= time, y = n_subs))+
  geom_line(aes(x= time, y = n_tracked), color = "red")

#as there is a lot of data, getting weekly mean sub size and mean number of inds tracked

# Extract week start dates
week_start <- as.Date(cut(mean_subs$time, "week"))

# Compute weekly means for subgroup size and n_tracked
weekly_summary <- aggregate(cbind(mean_subs, n_tracked) ~ week_start,
                            data = mean_subs, FUN = mean, na.rm = TRUE)

ggplot(weekly_summary, aes(x = as.Date(week_start))) +
  geom_line(aes(y = mean_subs, color = "Mean subgroup size"), size = 1) +
  geom_line(aes(y = n_tracked, color = "Mean individuals tracked"), size = 1, linetype = "dashed") +
  scale_x_date(date_labels = "%m-%y", date_breaks = "1 year") +
  labs(title = "Weekly Mean Subgroup Size ~ Number of Individuals Tracked",
       x = "Week", y = "Value") +
  scale_color_manual(values = c("Mean subgroup size" = "slateblue",
                                "Mean individuals tracked" = "orangered")) +
  theme_classic()

hist(weekly_summary$n_tracked, breaks = 20)

#filter subgroup counts to when there are at least 15 tracked to look at annual patterns in subgrouping associations

which(mean_subs$n_tracked == 15)[1]

#filter to the dates when there are first atleast 15 inds collecting data
mean_subs_filt <- mean_subs[which(mean_subs$n_tracked == 15)[1]:nrow(mean_subs),]
#get weekly summaries again with the subsetted data
week_start_filt <- as.Date(cut(mean_subs_filt$time, "week"))

weekly_summary_filt <- aggregate(cbind(mean_subs, n_tracked) ~ week_start_filt, data = mean_subs_filt, FUN = mean, na.rm = TRUE)


ggplot(weekly_summary_filt, aes(x = as.Date(week_start_filt))) +
  geom_line(aes(y = mean_subs, color = "Mean subgroup size"), size = 1) +
  geom_line(aes(y = n_tracked, color = "Mean individuals tracked"), size = 1, linetype = "dashed") +
  scale_x_date(date_labels = "%m-%y", date_breaks = "1 year") +
  labs(title = "Weekly Mean Subgroup Size ~ Number of Individuals Tracked",
       x = "Week", y = "") +
  scale_color_manual(values = c("Mean subgroup size" = "slateblue",
                                "Mean individuals tracked" = "orangered")) +
  theme_classic()

#don't think this gives much info, I wanted to see whether there are seasonal shifts in subgrouping patterns, but this might be confounded by the collars being on multiple prides and coalitions, I think to do this properly, I need to run this analysis per group to look at within grouping patterns over time
#also could look at amount of between group interactions over seasons

# Convert matrices to data frames for easier handling
R <- 600

split_subgroups_by_group_id <- list()

for(group_i in unique(all_ids$pride_id)){
  
  #group_i <- "southern_anabeb"
  
  group_index <- which(all_ids$pride_id == group_i)
  
  if (length(group_index) <= 1) next
    
  xs_group_i <- xs[group_index,]
  ys_group_i <- ys[group_index,]
  
  # For each column, check if all individuals have data (no NA)
  valid_cols <- colSums(!is.na(xs_group_i)) > 1
  
  if (all(!valid_cols)) next
  
  # First and last column indices where more than one has data
  first_col <- which(valid_cols)[1]
  last_col  <- tail(which(valid_cols), 1)
  
  # Subset to that range
  xs_group_i <- xs_group_i[, first_col:last_col]
  ys_group_i <- ys_group_i[, first_col:last_col]
  ts_i <- ts[first_col:last_col]
  
    #skip if only one individual in the group as can't do subgrouping patterns on them
  if(sum(table(all_ids$pride_id)[group_i]) == 1 || is.na(table(all_ids$pride_id)[group_i])) next
  
  subgroups <- get_subgroup_data(xs_group_i, ys_group_i, R)
  
  subgroups_df <- as.data.frame(subgroups$ind_subgroup_membership)
  colnames(subgroups_df) <- ts_i
  subgroups_df$Individual <- all_ids$name[group_index]
  
# Convert to long format
  subgroups_df <- subgroups_df %>%
    pivot_longer(cols = -Individual, names_to = "Time", values_to = "ind_subgroup_membership")

  n_subs_df <- data.frame(n_subgroups = subgroups$n_subgroups, mean_subsize = colMeans(subgroups$subgroup_counts, na.rm = T), Time = ts_i)
  
  #merge subgroups_df with n_subs_df
  #*****!!only issue with this is that it adds data to individuals who don't have individual subgroup membership!! - need to look into this
  df_combined <- merge(subgroups_df, n_subs_df, by = "Time")
  
  #put each groups subgrouping patterns in separate elements of a list
  split_subgroups_by_group_id[[group_i]] <- df_combined
  
}


#to get the mode subgroup ID per day
mode_subgroup <- function(x) {
  x <- x[!is.na(x)]           # ignore NAs
  if(length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#now can look at each groups subgrouping patterns separately!
for (group_name in names(split_subgroups_by_group_id)){

  group_i <- split_subgroups_by_group_id[[group_name]] 
  
#group_i <- split_subgroups_by_group_id$Okavariona
  group_i$Time <- as.POSIXct(group_i$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

#add a jitter so easier to see each ind in the plot
  group_i <- group_i %>%
  mutate(
    Individual_num = as.numeric(factor(Individual)),
    y_offset = (Individual_num - 1) * 0.02,  # adjust 0.01 as needed
    subgroup_jittered = as.numeric(ind_subgroup_membership), 
    x_offset = Time + (Individual_num - 1) * 3600 
  )

  g1 <- ggplot(data = group_i, aes(x= Time, y = subgroup_jittered, group = Individual, color = Individual))+
  geom_line(alpha = 0.6)+
  theme_classic()+
  labs(
    x = "Week",
    y = "Subgroup Membership",
    title = paste0("Raw Subgrouping Patterns at ", R, "m")
  ) +
  scale_x_datetime(date_labels = "%y-%m", date_breaks = "6 months")

  group_mode <- group_i %>%
  mutate(subgroup_num = as.numeric(ind_subgroup_membership),
         Week = floor_date(Time, unit = "1 day")) %>%        
  group_by(Individual, Week) %>%
    summarise(weekly_mode = mode_subgroup(subgroup_num), .groups = "drop")

#add a jitter so easier to see each ind in the plot
  group_mode <- group_mode %>%
  mutate(
    Individual_num = as.numeric(factor(Individual)),
    y_offset = (Individual_num - 1) * 0.01,  # adjust 0.01 as needed
    mode_jittered = weekly_mode + y_offset,
    x_offset = Week + (Individual_num - 1) * 3600 
  )

  g2 <- ggplot(group_mode, aes(x = x_offset, y = mode_jittered, color = Individual)) +
  geom_line(size = 1, alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Week",
    y = "Mode Subgroup Membership",
    title = "Daily Smoothing"
  ) +
  scale_x_datetime(date_labels = "%y-%m", date_breaks = "6 months")

  gg <- g1/g2


ggsave(filename = paste0(plot_dir, "subgrouping_patterns_", group_name, "_", R, ".png"), plot = gg, width = 12, height = 7, dpi = 1200)

}



#now looking at the number of subgroups over time for each group
for (group_name in names(split_subgroups_by_group_id)) {
  
  #group_name <- "Okavariona"
  group_i <- split_subgroups_by_group_id[[group_name]]
 
  group_i$Time <- as.POSIXct(group_i$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  group_i$only_time <- as_hms(group_i$Time)
  
  #if only_time is NA, fill with 00:00:00 for midnight
  group_i$only_time[is.na(group_i$only_time)] <- "00:00:00"
  group_i$only_time <- as_hms(group_i$only_time)
  
  #summarise the number of subgroups and number of inds tracked per day (mode)
  summary_i <- group_i %>%
    mutate(Date = floor_date(Time, unit = "day")) %>%
    group_by(Date) %>%
    summarise(
      n_tracked = mode_subgroup(ind_subgroup_membership),  # mode of subgroup membership
      n_subgroups = mode_subgroup(n_distinct(ind_subgroup_membership[!is.na(ind_subgroup_membership)])), # mode of subgroup counts
      .groups = "drop"
    )
  
  # Add 7-day rolling mean
  summary_i <- summary_i %>%
    arrange(Date) %>%
    mutate(
      n_tracked_7d = rollmean(n_tracked, k = 7, fill = NA, align = "right"),
      n_subgroups_7d = rollmean(n_subgroups, k = 7, fill = NA, align = "right")
    )
  #remove NA 
  summary_i <- summary_i[!is.na(summary_i$Date),]
  
  p <- ggplot(summary_i, aes(x = Date)) +
    geom_line(aes(y = n_tracked, color = "No. inds tracked"), alpha = 0.4) +
    geom_line(aes(y = n_tracked_7d, color = "No. inds tracked"), size = 1) +
    geom_line(aes(y = n_subgroups, color = "Mode subgroup size"), alpha = 0.4) +
    geom_line(aes(y = n_subgroups_7d, color = "Mode subgroup size"), size = 1) +
    scale_color_manual(name = "Metric", values = c("No. inds tracked" = "slateblue", "Mode subgroup size" = "gold2")) +
    theme_classic() +
    labs(y = "Count",
      title = paste0("Subgroup size when R is ", R, "m for ", group_name))
  # Save plot
  ggsave(filename = paste0(plot_dir, "subgroupsummary_", group_name, ".png"),
         plot = p, width = 12, height = 4, dpi = 1200)
  
  #look at mean subgroup size and number of subgroups across time of day
  
  #because the times that data are collected are reduced to every 4h during the day, need to think about how to filter the data to remove these points
  #for now, filtering data to only 4h during the day (so recording at 6am,10am, 2pm, 6pm) - but need to think about this more as we are removing data
  
  # group_i <- group_i[!group_i$only_time == as_hms("08:00:00"),]
  # group_i <- group_i[!group_i$only_time == as_hms("12:00:00"),]
  # group_i <- group_i[!group_i$only_time == as_hms("16:00:00"),]
  # group_i <- group_i[!group_i$only_time == as_hms("04:00:00"),]
  # group_i <- group_i[!group_i$only_time == as_hms("20:00:00"),]
  
  #this filtering didn't fix the issue, I think that the gps points during the day need to be shifted slightly so that they are all on the same schedule, or that GPS points could be interpreted to every 2h during the day before subgroup measures are extracted
  
  plot(group_i$mean_subsize, group_i$n_subgroups)
  
  p0 <- ggplot(group_i, aes(x = only_time, y = n_subgroups)) +
    geom_smooth(se = FALSE, method = "loess", color = "steelblue2", size = 1.2) +
    #geom_point()+
    theme_classic() +
    xlab("Time of Day") +
    ylab("Number of Subgroups") +
    labs(title = "Smoothed Number of Subgroups Throughout the Day")
  
  p1 <- ggplot(group_i, aes(x = only_time, y = mean_subsize)) +
    geom_smooth(se = FALSE, method = "loess", color = "steelblue3", size = 1.2) +
    #geom_point()+
    theme_classic() +
    xlab("Time of Day") +
    ylab("Subgroups Size") +
    labs(title = "Smoothed Subgroup Size Throughout the Day")
  
  
  # Summarize mean and standard error
  agg_df <- group_i %>%
    group_by(only_time) %>%
    summarise(
      mean_n_subs = mean(n_subgroups, na.rm = TRUE),
      mean_size = mean(mean_subsize, na.rm = TRUE),
      se_n_subs = sd(n_subgroups, na.rm = TRUE) / sqrt(n()),
      se_mean_size = sd(mean_subsize, na.rm = TRUE) / sqrt(n())
    )
  
  # Plot with error bars
  p2 <- ggplot(agg_df, aes(x = only_time, y = mean_n_subs)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_n_subs - se_n_subs,
                      ymax = mean_n_subs + se_n_subs),
                  width = 0.3, color = "black") +
    theme_classic() +
    xlab("Time of Day (UTC)") +
    ylab("Mean Number of Subgroups") +
    labs(title = "Mean Number of Subgroups per Hour with SE") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  p3 <- ggplot(agg_df, aes(x = only_time, y = mean_size)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_size - se_mean_size,
                      ymax = mean_size + se_mean_size),
                  width = 0.3, color = "black") +
    theme_classic() +
    xlab("Time of Day (UTC)") +
    ylab("Mean Subgroup Size") +
    labs(title = "Mean Subgroup Size per Hour with SE") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  
  #gg <- (p0 + p1) / (p2 + p3)
  #these plots don't show the pattern expected whereby when the group size is larger there should be fewer subgroups, I think the issue lies in the 4h resolution during the day
  
}

#---------------------------------------------------------------

#daily changes in subgroup size per group
na_count <- list()

#look at which times there is data for each individual
for (group_name in names(split_subgroups_by_group_id)){

  group_i <- split_subgroups_by_group_id[[group_name]]
  group_i$Time <- as.POSIXct(group_i$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  group_i$only_time <- as_hms(group_i$Time)
  #if only_time is NA, fill with 00:00:00 for midnight
  group_i$only_time[is.na(group_i$only_time)] <- "00:00:00"
  group_i$only_time <- as_hms(group_i$only_time)
  
  for (ind in unique(group_i$Individual)){
    
    ind_i <- group_i[group_i$Individual == ind,]
    ind_i_summary <- ind_i %>%
      group_by(only_time) %>%
      summarise(
       sum_nas = length(which(is.na(ind_subgroup_membership))),
       sum_data = length(which(!is.na(ind_subgroup_membership)))
    )
      
    #ggplot(ind_i_summary, aes(x = only_time, y = sum_data))+
    #  geom_count()
    
    na_count[[ind]] <- ind_i_summary
  }
}


#plot the number of data points collected for each time point for each individual
na_count_df <- NA
for (ind in names(na_count)){
 na_count_df_ind <- data.frame(na_count[[ind]])
 na_count_df_ind$ind <- ind
 na_count_df <- rbind(na_count_df, na_count_df_ind)
}
#remove NA row
na_count_df <- na_count_df[!is.na(na_count_df$only_time),]

#add pride ID to change colors in plot to pride rather than ind
na_count_df <- merge(all_ids, na_count_df, by.x = "name", by.y = "ind")

timeplot <- ggplot(na_count_df, aes(x = only_time, y= sum_data, group = pride_id, fill = pride_id))+
  geom_bar(stat = "identity", position = "stack")+
  xlab("Time of Day (UTC)") +
  ylab("Sum of GPS data points") +
  theme_classic()

ggsave(filename = paste0(plot_dir, "countofdatapergroup_eachhour.png"),
       plot = timeplot, width = 11, height = 8, dpi = 1200)

timeplot2 <- ggplot(na_count_df, aes(x = only_time, y= sum_data, fill = pride_id))+
  geom_bar(stat = "identity", position = "stack")+
  xlab("Time of Day (UTC)") +
  ylab("Sum of GPS data points") +
  facet_wrap(~name, ncol = 3)+
  theme_classic()

ggsave(filename = paste0(plot_dir, "countofdatapergroup_eachhour_perind.png"),
       plot = timeplot2, width = 14, height = 14, dpi = 1200)


#why is NPL-32 missing data?? - have I got the ID order right?
#this ind has data but I think when I clean it, it somehow gets lost

