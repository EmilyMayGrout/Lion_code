#this script is using the cleaned 2h/4h sampling regime to investigate temporal patterns in subgrouping for periods of time when all individuals are tracked per pride/coalition


library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(zoo)
library(fields)
library(patchwork)
library(hms)
library(cocomo)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
plot_dir <-  'C:/Users/egrout/Dropbox/lion/results/level1/'

setwd(code_dir)
source('coati_function_library_V1.R')
source('lion_functions.R')

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
load('allnamib_xy_2hour_level1.RData') #level 1 is filtering to gps at the correct 2h/4h interval
load('allnamib_ids.RData')
all_nab <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/all_namibia.csv")

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


#----------------------------------------------------------------------------
#look at Omukutu 2h/4h subsampling - why are midday points not there when all are tracked?
#----------------------------------------------------------------------------

Omukutu_xs <- xs[which(all_ids$pride_id == "Omukutu"),]
Omukutu_ys <- ys[which(all_ids$pride_id == "Omukutu"),]

Omukutu_ids <- all_ids[which(all_ids$pride_id == "Omukutu"),]
individual_labels <- Omukutu_ids$name

# --- 1. Convert to data frame and add timestamps correctly ---
Omukutu_df <- as.data.frame(Omukutu_xs)
ts_fixed <- format(ts, "%Y-%m-%d %H:%M:%S")
colnames(Omukutu_df) <- ts_fixed

# --- 2. Pivot longer (each individual × timestamp) ---
Omukutu_long <- Omukutu_df %>%
  mutate(individual = rownames(Omukutu_df)) %>%
  pivot_longer(
    cols = -individual,
    names_to = "ts",
    values_to = "x_value"
  )

Omukutu_long$ts <- as.POSIXct(Omukutu_long$ts, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")

# --- 4. Extract hour and month ---
Omukutu_long <- Omukutu_long %>%
  mutate(
    hour = as.numeric(strftime(ts, "%H", tz = "UTC")),
    year_month = format(ts, "%Y-%m"),
    has_data = !is.na(x_value)
  )

# --- 5. Count number of data points per hour per individual per month ---
Omukutu_long <- Omukutu_long %>%
  filter(year_month >= "2022_05")

hourly_density <- Omukutu_long %>%
  group_by(year_month, individual, hour) %>%
  dplyr::summarise(
    n_points = sum(has_data, na.rm = TRUE),
    .groups = "drop"
  )

# --- 6. Plot ---
ggplot(hourly_density, aes(x = factor(hour), y = n_points, fill = individual)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ year_month, ncol = 3) +
  theme_classic() +
  labs(
    title = "Omukutu Pride – Data Density by Hour (May 2022 and After)",
    x = "Hour (UTC)",
    y = "Number of Data Points",
    fill = "Individual ID"   # This sets the legend title
  )+
  scale_fill_manual(
    name = "Lion Name",
    values = brewer.pal(n = length(individual_labels), name = "Set1"),  # optional: assign colors
    labels = individual_labels)

#----------------------------------------------------------------------------




R <- 60

#how many days are all group members from each pride simultaneously tracked?

#make an empty dataframe
pride_summary <- data.frame(
  pride_id = character(),
  n_days_tracked = numeric(),
  first_date = as.Date(character()),
  last_date = as.Date(character()),
  mean_pct_tracked = numeric(),
  n_inds = numeric(),
  stringsAsFactors = FALSE
)


#i <- "Omukutu"

for(i in unique(all_ids$pride_id)){

xs_pride_i <- xs[which(all_ids$pride_id == i), , drop = FALSE] #this keeps it as a matrix even if its 1 individual
ys_pride_i <- ys[which(all_ids$pride_id == i), , drop = FALSE]
n_inds_i <- nrow(xs_pride_i)

#skip if only one individual
if (n_inds_i <= 1) next

n_tracked_i <- colSums(!is.na(xs_pride_i))

#if there are more than 6 individuals in the group, I want the number of times all but two are tracked
if (n_inds_i <= 6) {
  # Require all individuals tracked
  all_tracked_idxs_i <- which(n_tracked_i == n_inds_i)
} else {
  all_tracked_idxs_i <- which(n_tracked_i >= (n_inds_i - 2))
}

n_days_tracked <- length(unique(as.Date(ts[all_tracked_idxs_i])))
mean_pct_tracked <- mean(n_tracked_i[all_tracked_idxs_i] / n_inds_i) * 100
first_date <- min(unique(as.Date(ts[all_tracked_idxs_i])))
last_date <- max(unique(as.Date(ts[all_tracked_idxs_i])))

pride_summary <- rbind(
  pride_summary,
  data.frame(pride_id = i, n_inds = n_inds_i, mean_pct_tracked = mean_pct_tracked, first_date = first_date,
    last_date = last_date, n_days_tracked = n_days_tracked, stringsAsFactors = FALSE))
}

rm(xs_pride_i, ys_pride_i)


#for each day, how many inds are tracked
n_tracked_df <- data.frame(datetime = ts, n_tracked = colSums(!is.na(xs)))
ggplot(n_tracked_df, aes(x = datetime, y = n_tracked)) +
  geom_line(color = "steelblue") +
  labs(x = "Time", y = "Number of Individuals Tracked Simultaneously", title = "Tracking Coverage Over Time") +
  theme_classic()



#subset times when all inds are tracked to look at subgrouping patterns


subgroup_data_perpride <- list()
xs_perpride <- list()

for(i in unique(pride_summary$pride_id)){
  
  #filter to date range where all individuals are tracked
  first_ts_indx <- which(ts == as.Date(pride_summary$first_date[pride_summary$pride_id == i]))
  last_ts_indx  <- which(ts == as.Date(pride_summary$last_date[pride_summary$pride_id == i]))
  
  if(length(first_ts_indx) == 0) next
  
  #index xs and ys for each prides individuals
  xs_pride_i <- xs[which(all_ids$pride_id == i),first_ts_indx:last_ts_indx , drop = FALSE]#this keeps it as a matrix even if its 1 individual
  ys_pride_i <- ys[which(all_ids$pride_id == i),first_ts_indx:last_ts_indx , drop = FALSE]
  ts_pride_i <- ts[first_ts_indx:last_ts_indx]

  subgroup_data_i <- get_subgroup_data(xs_pride_i, ys_pride_i, R)
  #add time to list
  subgroup_data_i$Time <- ts_pride_i
  subgroup_data_i$orig_ts <- first_ts_indx:last_ts_indx
  
  subgroup_data_perpride[[i]] <- subgroup_data_i
  
  colnames(xs_pride_i) <- first_ts_indx:last_ts_indx
  xs_perpride[[i]] <- xs_pride_i
  
  
}



#now looking at the number of subgroups over time for each group
for (group_name in names(subgroup_data_perpride)) {
  
  #group_name <- "Uniab"
  group_i <- subgroup_data_perpride[[group_name]]
  
  group_i$Time <- as.POSIXct(group_i$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  group_i$only_time <- as_hms(group_i$Time)
  
  #if only_time is NA, fill with 00:00:00 for midnight
  group_i$only_time[is.na(group_i$only_time)] <- "00:00:00"
  group_i$only_time <- as_hms(group_i$only_time)
  
  group_i_df <- data.frame(Time = group_i$Time, only_time = group_i$only_time, n_subgroups = group_i$n_subgroups)
  group_i_df$mean_subsize <- colMeans(group_i$subgroup_counts, na.rm = T)
  
  plot1 <- ggplot(group_i_df, aes(x = only_time)) +
    geom_smooth(aes(y = n_subgroups, color = "Number of Subgroups"),
                se = F, method = "loess", linewidth = 1.2, span = 0.4) +
    geom_smooth(aes(y = mean_subsize, color = "Mean Subgroup Size"),
                se = F, method = "loess", linewidth = 1.2, span = 0.4) +
    scale_y_continuous(name = "Count") +
    scale_color_manual(values = c("Number of Subgroups" = "steelblue2",
                                  "Mean Subgroup Size" = "orangered")) +
    theme_classic() +
    xlab("Time of Day") +
    labs(title = paste0("Smoothed Number of Subgroups and Subgroup Size for ", group_name),
         color = "Metric")
  
  #subgroup dynamics over time all inds are collared
  plot2 <- ggplot(group_i_df, aes(x = Time))+
    geom_smooth(aes(y = mean_subsize, color = "Mean Subgroup Size"))+
    geom_smooth(aes(y = n_subgroups, color = "Number of Subgroups"))+
    scale_color_manual(values = c("Number of Subgroups" = "steelblue2",
                                  "Mean Subgroup Size" = "orangered"))+
    scale_y_continuous(name = "Count") +
    theme_classic()+
    xlab("Date") +
    labs(title = paste0("Smoothed Number of Subgroups and Subgroup Size for ", group_name),
         color = "Metric")
  
  plot <- plot1 + plot2
  ggsave(filename = paste0(plot_dir, "subgroupsize_", group_name, "_", R, ".png"), plot = plot, width = 18, height = 7, dpi = 300)
  
  
  #also plotting the prop time together per month
  n_inds <- nrow(group_i$ind_subgroup_membership)
  ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
  split_by_month <- T
  year_months <- unique(format(group_i$Time, "%Y-%m"))
  ff_net_month <- array(NA, dim = c(n_inds, n_inds, length(year_months)))
  dimnames(ff_net_month)[[3]] <- year_months
  
  #going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      
      if(split_by_month){
        for(m in 1:length(seq_along(year_months))){
          
          ym <- year_months[m]
          
          # Logical index for matching year-month
          ts_filter <- format(group_i$Time, "%Y-%m") == ym
          
          sub_ids_i <- group_i$ind_subgroup_membership[i, ts_filter]
          sub_ids_j <- group_i$ind_subgroup_membership[j, ts_filter]
          
          ff_net_month[i, j, m] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
          #print(ff_net_month[i, j, m])
        }
      }
      
      #getting subgroup id for individual i and j
      sub_ids_i <- group_i$ind_subgroup_membership[i,]
      sub_ids_j <- group_i$ind_subgroup_membership[j,]
      
      #computing edge weight (fraction of time in same subgroup)
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
    }
  }
  
  diag(ff_net) <- NA
  
  ids <- all_ids[which(all_ids$pride_id == group_name),]
  
  for (i in 1:dim(ff_net_month)[3]) {
    diag(ff_net_month[,,i]) <- NA
  }
  n_months <- dim(ff_net_month)[3]
  
  n_row <- if (n_months < 4) {
    1
  } else {
    2
  }
  
  png(height = 400*n_row, width = 500*n_months/n_row, units = 'px', filename = paste0(plot_dir,'subgroup_network', R, group_name, '_monthly.png'))
  
  par(mfrow=c(n_row,(n_months/n_row))) 
  
  for(i in 1:dim(ff_net_month)[3]){
    ffnet_m <- ff_net_month[, , i]
    visualize_lion_network(ffnet_m, ids)
    mtext(line=1, side=3, paste(year_months[i], group_name), outer=F)
  }
  dev.off()
}



#-----------------------------------------------------------------


#where are the distributions of the prides? - how much overlap is there in space use, who interacts, when, and where?


# group_name <- names(subgroup_data_perpride)[1]
# xs_group_i <- xs[which(all_ids$pride_id == group_name),]
# ys_group_i <- ys[which(all_ids$pride_id == group_name),]

all_sub_data <- get_subgroup_data(xs, ys, 600)

#get the travel distance for each tracked individual over time

xs_1 <- xs[1,]

n_inds <- nrow(xs)
first_ts_indx <- which(ts == as.Date(pride_summary$first_date[pride_summary$pride_id == i]))
last_ts_indx  <- which(ts == as.Date(pride_summary$last_date[pride_summary$pride_id == i]))




plot(NULL, xlim=c(min(xs[,],na.rm=T), max(xs[,],na.rm=T)), 
     ylim=c(min(ys[,],na.rm=T), max(ys[,],na.rm=T)), asp=1) 

for(i in seq(1, n_inds)){	
  points(xs[i,], ys[i,], 
         col=rainbow(n_inds)[i], pch=19, cex=0.1, asp=1)
}



#alternative plotting method: 

#plot the locations of all inds onto the map
namibia <- ne_countries(country = "Namibia", returnclass = "sf")

lon_min <- quantile(all_nab$location.long, 0.001, na.rm = TRUE)
lon_max <- quantile(all_nab$location.long, 0.999, na.rm = TRUE)
lat_min <- quantile(all_nab$location.lat, 0.001, na.rm = TRUE)
lat_max <- quantile(all_nab$location.lat, 0.999, na.rm = TRUE)

xlim <- c(lon_min, lon_max) + c(0, 1)
ylim <- c(lat_min, lat_max)

all_nab <- all_nab %>%
  left_join(all_ids %>% select(name, pride_id), 
            by = c("individual.local.identifier" = "name"))

unique_prides <- unique(all_nab$pride_id)
pride_colors <- rainbow(length(unique_prides))
names(pride_colors) <- unique_prides
alpha_cols <- adjustcolor(pride_colors, alpha.f = 0.05)
names(alpha_cols) <- unique_prides   # <- assign names here



png(width = 7, height = 6, units = "in", res = 1200, filename = paste0(plot_dir,'grouplocations.png'))


plot(st_geometry(namibia), xlim = xlim, ylim = ylim, col = "gray95", border = "black",
     asp = 1, xlab = "Longitude", ylab = "Latitude", main = "")

for(pride in unique_prides){
  points(
    all_nab$location.long[all_nab$pride_id == pride],
    all_nab$location.lat[all_nab$pride_id == pride],
    col = alpha_cols[pride],
    pch = 20,
    cex = 0.8
  )
}

legend("topright", legend = unique_prides, col = pride_colors, pch = 19, cex = 1, bg = "white")

dev.off()


table(all_ids$group_type, all_ids$pride_id)




#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#
#when do splits and merges occurs per pride?
#
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#subgroup_data_perpride is the subset times when all inds are tracked 
splits_list <-  vector("list", length(subgroup_data_perpride))
names(splits_list) <- names(subgroup_data_perpride)
merges_list <-  vector("list", length(subgroup_data_perpride))
names(merges_list) <- names(subgroup_data_perpride)

for (i in seq_along(subgroup_data_perpride)) {
  
  subgroup_memb <- subgroup_data_perpride[[i]]$ind_subgroup_membership
  
  splits_df <- detect_splits(subgroup_memb)
  splits_df$Time <-  subgroup_data_perpride[[i]]$Time[splits_df$t]
  splits_df$hour <- as.numeric(strftime(splits_df$Time, "%H", tz = "UTC"))
  splits_list[[i]] <- splits_df
  
  merges_df <- detect_merges(subgroup_memb)
  merges_df$Time <-  subgroup_data_perpride[[i]]$Time[merges_df$t]
  merges_df$hour <- as.numeric(strftime(merges_df$Time, "%H", tz = "UTC"))
  merges_list[[i]] <- merges_df
  
  
}

all_pride_splits <- ldply(splits_list, data.frame)
all_pride_splits$event <- "split"
all_pride_splits_cut <- all_pride_splits[,c(1,16,17)]
all_pride_merges <- ldply(merges_list, data.frame)
all_pride_merges$event <- "merge"
all_pride_merges_cut <- all_pride_merges[,c(1,16,17)]

all_pride_split_merge <- rbind(all_pride_splits_cut, all_pride_merges_cut)

gg <- ggplot(all_pride_split_merge, aes(x = hour, fill = .id)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  scale_x_continuous(breaks = 0:23) +
  theme_classic(base_size = 18) +
  labs(
    title = "Density of Events by Hour",
    x = "Hour of Day (UTC)",
    y = "Density"
  ) + 
  scale_fill_discrete(name = "Pride/Coalition ID")+
  facet_wrap(~event)
gg

ggsave(filename = paste0(plot_dir, "splitmerge_density", "_", R, ".png"), plot = gg, width = 18, height = 7, dpi = 300)



#-----------------------------------------------------------------------------------------
#How much of this pattern is driven by the amount of data available??
#-----------------------------------------------------------------------------------------

#go through the xs_perpride I made in an earlier forloop to count which hours have data for each pride

all_pride_long_list <- list()


for(i in seq_along(xs_perpride)){
  
xs_pride_i <- xs_perpride[[i]]

col_ts_idx <- as.numeric(colnames(xs_pride_i))
hour_row <- as.numeric(strftime(ts[col_ts_idx], "%H", tz = "UTC"))

#in case I rerun this code, it stops adding rows 
if ("hour" %in% rownames(xs_pride_i)) {
  xs_pride_i["hour", ] <- hour_row
} else {
  xs_pride_i <- rbind(xs_pride_i, hour_row)
  rownames(xs_pride_i)[nrow(xs_pride_i)] <- "hour"
}
xs_perpride[[i]] <- xs_pride_i

xs_pride_i_long <- as.data.frame(t(xs_pride_i)) 

longer_df <- xs_pride_i_long %>%
  pivot_longer(cols = starts_with("V"), names_to = "individual", values_to = "value")

all_pride_long_list[[i]] <- longer_df


}

#now plotting the amount of data collected per hour for each ind in each pride - but first need to get the list of xs into a long dataframe
names(all_pride_long_list) <- names(xs_perpride)

all_pride_long_df <- ldply(all_pride_long_list, data.frame)

ggplot(all_pride_long_df, aes(x = hour, y = value, fill = individual)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_classic()+
  facet_wrap(~.id)+
  theme(legend.position="none")


#how many NAs are there per hour?

#counting the number of NAs per ind per hour
all_pride_long_df$hour <- as.integer(all_pride_long_df$hour)
prides <- unique(all_pride_long_df$.id)
hours <- 0:23

#create all combinations of pride × individual × hour
all_combinations <- expand.grid(
  pride = prides,
  individual = unique(all_pride_long_df$individual),
  hour = hours,
  stringsAsFactors = FALSE)

all_combinations$n_na <- 0
all_combinations$n_total <- 0
all_combinations$prop_na <- 0

# Loop through each row and count NAs
for (i in seq_len(nrow(all_combinations))) {
  p <- all_combinations$pride[i]
  ind <- all_combinations$individual[i]
  h <- all_combinations$hour[i]

  subset_data <- all_pride_long_df[
    all_pride_long_df$.id == p &
      all_pride_long_df$individual == ind &
      all_pride_long_df$hour == h,]

  n_na <- sum(is.na(subset_data$value))
  n_total <- nrow(subset_data)
  prop_na <- if (n_total > 0) n_na / n_total else NA
  
  all_combinations$n_na[i] <- n_na
  all_combinations$n_total[i] <- n_total #number of obs available for that hour, ind, pride
  all_combinations$prop_na[i] <- prop_na
}

#restructure df to long format to show the total data vs NA data
plot_long <- plot_data %>%
  pivot_longer(
    cols = c(n_total, n_na),
    names_to = "type",
    values_to = "count"
  ) %>%
  mutate(type = factor(type, levels = c("n_total", "n_na"), labels = c("Available Data", "NA")))


#only show the hours which are at the 2/4h schedule
hours_to_show <- c(6, 10, 14, 18, 20, 22, 0, 2, 4)
plot_data <- subset(all_combinations, hour %in% hours_to_show)

#NAs on a pride level
g1 <- ggplot(plot_long, aes(x = factor(hour), y = count, fill = type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ pride) +
  scale_fill_manual(
    name = "Data Type",
    values = c("Available Data" = "skyblue", "NA" = "red")
  ) +
  theme_classic() +
  labs(
    x = "Hour (UTC)",
    y = "Number of Data Points",
    title = "Data Availability per Pride and Hour"
  )

ggsave(filename = paste0(plot_dir, "data_available_perpride_", R, ".png"), plot = g1, width = 25, height = 15, dpi = 300)


#NA's on an individual level

g2 <- ggplot(plot_long, aes(x = factor(hour), y = count, fill = type, group = individual)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity", show.legend = TRUE) +
  facet_wrap(~ pride) +
  scale_fill_manual(
    name = "Data Type",
    values = c("Available Data" = "skyblue", "NA" = "red")
  ) +
  theme_classic(base_size = 25) +
  labs(
    x = "Hour (UTC)",
    y = "Number of Data Points",
    title = "Data Availability per Individual"
  ) +
  guides(fill = guide_legend(title = "Data Type"))  # ensures only NA vs Available is shown

ggsave(filename = paste0(plot_dir, "data_available_perind_", R, ".png"), plot = g2, width = 25, height = 15, dpi = 300)

