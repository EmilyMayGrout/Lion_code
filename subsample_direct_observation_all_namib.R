#this script is using all the namibia groups to run the subsampling analysis to simulate traditional field methods
#compare accuracy of results when visiting at 6am, 10am, and 2pm

library(hms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fields)
library(asnipe)
library(fields)

code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
plot_dir <-  'C:/Users/egrout/Dropbox/lion/results/level2/'

setwd(code_dir)
source('coati_function_library_V1.R')
source('lion_functions.R')

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
#load('allnamib_xy_2hour_level2.RData') #level 1 is filtering to gps at the correct 2h/4h interval, getting the closest fix within 20 mins of the desired hours
load('allnamib_omukutu_xy_2hour_level3.RData')

load('allnamib_ids_withdod_V2.RData')


#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

R <- 60

# add 0s for xs and ys after lions died - so they're not missing but are not with the group, so getting times when all tracked would not be affected by an animal that died 
#if the dod is known, add 0's after that, if unknown, add 0's when there is no more collar data
add_0s <- F

if(add_0s == T){

for (id in seq_len(nrow(all_ids))) {
  
  known <- all_ids$known[id]
  if (is.na(known)) { #skip if known is NA
    next
  }
  xs_vec <- xs[id, ]
  ys_vec <- ys[id, ]
  
  if (known) {
    dod <- all_ids$DOD_clean[id]
    
    if(dod > max(ts)){
      next
    }
    
    #get ts for date of death
    dod_idx <- which(ts == dod)
    
    # set everything AFTER DOD to 0
    xs[id, dod_idx:length(ts)] <- 0
    ys[id, dod_idx:length(ts)] <- 0
    
  } else { #when the dod isn't known
    #get last NON-NA value in xs or ys
    last_fix_idx <- max(which(!is.na(xs_vec) | !is.na(ys_vec)))
    
    # Skip if last fix is the end of the ts
    if (last_fix_idx == length(ts)) {
      next
    } 
    
    # If all NA, skip to avoid errors
    if (is.finite(last_fix_idx)) {
      xs[id, (last_fix_idx+1):length(ts)] <- 0
      ys[id, (last_fix_idx+1):length(ts)] <- 0
    }
  }
}
}


#look at general patterns of association across groups
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
new_order <- c(1:42)
ffnet_reorder <- ff_net[new_order, new_order]

make_plot <- T
if(make_plot == T){
png(height = 1100, width = 1200, units = 'px', filename = paste0(plot_dir,'allnamib_network_60m_2yr_withdod_V1.png'))

visualize_lion_network(ffnet_reorder, all_ids[new_order,])

dev.off()
}

#-------------------------------------------------------------------------------------
#splitting data to each pride/coalition association to simulate traditional field observations 
#-------------------------------------------------------------------------------------

unique_groups <- unique(unlist(strsplit(all_ids$pride_id, "-")))

group_matrices <- lapply(unique_groups, function(p) {
  inds <- which(grepl(p, all_ids$pride_id))
  list(
    group_id = p,
    xs = xs[inds, , drop = FALSE],
    ys = ys[inds, , drop = FALSE],
    ids = all_ids$name[inds],
    id_info = all_ids[inds,]
  )
})

names(group_matrices) <- unique_groups

#group_i <- "Omukutu"


#making plots for data presence per month (whether there was any data collected in that month for each group, and add grey in plot when the animal died - if the dod is unknown, assumed dod was when there were no more data). Aim to find the amount of time where most/all inds in a group were tracked, excluding the inds which died
for(group_i in unique(unique_groups)){

group <- group_matrices[[group_i]]

#excluding the individuals in this calculation after they have died by assigning whether they have data/no data/dead in a seperate matrix:
result <- mask_after_death_with_presence(
  xs = group$xs,
  ys = group$ys,
  ts = ts,
  id_info = group$id_info
)

presence <- result$presence
n_tracked <- result$n_tracked

df <- data.frame(n_tracked, ts)

#ggplot(df, aes(x = ts, y = n_tracked))+
#  geom_point()


#converting presence matrix to long format for plotting
presence_long <- as.data.frame(result$presence)
rownames(presence_long) <- group$ids

presence_long$individual <- rownames(presence_long)

# Convert to long format
presence_long <- result$presence %>%
  as.data.frame() %>%
  mutate(individual = as.numeric(rownames(.))) %>%
  pivot_longer(
    cols = -individual,
    names_to = "time_idx",
    values_to = "presence"
  ) %>%
  mutate(
    time_idx = as.integer(gsub("[^0-9]", "", time_idx)), #as the colnames were V1, V2, need to extract the numbers
    ts = ts[time_idx],  #assign timestamp from ts vector
    month_year = as.Date(paste0(format(ts, "%Y-%m"), "-01"))   
  )

presence_long <- presence_long %>%
  mutate(
    individual = group$ids[individual],  # ensures correct row order
    status = case_when(
      is.na(presence) ~ "After death",
      presence == 1 ~ "Data",
      presence == 0 ~ "Missing"
    ),
    status = factor(status, levels = c("Data", "Missing", "After death"))
  )

#make a plot showing if any data were present in that month

data_presence <- ggplot(presence_long, aes(x = month_year, y = individual, fill = status)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Data" = "darkolivegreen3",
      "Missing" = "lightblue1",
      "After death" = "grey"
    )
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Month-Year",
    y = "Individual",
    fill = "Status",
    title = paste0("Data Presence for Pride: ", group$group_id)
  )


#generally, how much data is there per month?

#calculating the proportion of data points per individual per month
monthly_data <- presence_long %>%
  group_by(individual, month_year) %>%
  summarise(
    data_points = sum(status == "Data", na.rm = TRUE),
    total_points = n(),
    dead = any(status == "After death"),
    .groups = "drop"
  ) %>%
  mutate(
    prop_data = ifelse(dead, NA, data_points / total_points)  # NA for dead months
  )

# Plot with darkness representing amount of data
data_prop <- ggplot(monthly_data, aes(x = month_year, y = individual)) +
  geom_tile(aes(fill = prop_data), color = "white") +
  geom_tile(data = monthly_data %>% filter(dead),
            fill = "grey", color = "white") +
  scale_fill_gradient(low = "white", high = "darkolivegreen3", na.value = "grey", name = "Data Proportion") +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Month-Year",
    y = "Individual",
    title = paste0("Monthly Data Coverage with Death for Pride: ", group$group_id)
  )


#ggsave(filename = paste0(plot_dir, "datapresence2", "_", group_i, ".png"), plot = data_presence, width = 10, height = 7, dpi = 300)

#ggsave(filename = paste0(plot_dir, "dataprop2", "_", group_i, ".png"), plot = data_prop, width = 10, height = 7, dpi = 300)

}


#next todo is find the timestamps for each group where all are tracked at the same time - how much data coverage does this give? Do I need to change the time difference from the correct time stamps to 60mins rather than 20mins to get more data?

#make empty df to get the dates where all inds which are alive are tracked
dur_df <- data.frame(group_name = names(group_matrices))

for(group_i in names(group_matrices)){

  group <- group_matrices[[group_i]]

  if(nrow(group$xs) == 1) next

  result <- mask_after_death_with_presence(
    xs = group$xs,
    ys = group$ys,
    ts = ts,
    id_info = group$id_info
  )
  
# 1 = data present, 0 = missing, NA = after death
presence <- result$presence
n_inds_i <- nrow(presence)
n_times <- ncol(presence)

#find times where every alive individual has data in the presence matrix, going row by row)
all_alive_tracked <- sapply(1:n_times, function(t) {
  alive <- !is.na(presence[, t]) #alive individuals (regardless of whether they have data)
  if (!any(alive)) return(FALSE) #skip if no one alive
  all(presence[alive, t] == 1)   #all alive individuals must have presence == 1
})

# extract time indices
all_alive_tracked_idx <- which(all_alive_tracked)

#timestamps where everyone alive was tracked
times_all_tracked <- ts[all_alive_tracked]

# number of such time points
num_timepoints_all_alive_tracked <- length(all_alive_tracked)

#number of days where there is a GPS point during the day
hours_all_tracked <- hour(times_all_tracked)

# keep only 6, 10, or 14
day_hours <- c(6, 10, 14)
night_hours <- c(18, 20, 22, 0, 2, 4)

day_hours_filt <- times_all_tracked[hours_all_tracked %in% day_hours]
night_hours_filt <- times_all_tracked[hours_all_tracked %in% night_hours]

# Extract the dates for day and night fixes
day_dates <- unique(as.Date(day_hours_filt))
night_dates <- unique(as.Date(night_hours_filt))

# Find dates where there is both day and night data
dates_with_day_and_night <- as.Date(intersect(day_dates, night_dates), origin = "1970-01-01")
dates_with_day_and_night

dur <- as.numeric(max(dates_with_day_and_night) - min(dates_with_day_and_night)) 
min_date <- min(dates_with_day_and_night)
max_date <- max(dates_with_day_and_night)

dur_df$dates_day_night[dur_df$group_name == group_i] <- dur
dur_df$min_date_daynight[dur_df$group_name == group_i] <- as.character(min_date)
dur_df$max_date_daynight[dur_df$group_name == group_i] <- as.character(max_date)
dur_df$all_dates[dur_df$group_name == group_i] <- paste(as.character(dates_with_day_and_night), collapse = ", ")


}

#put NA for inf values
dur_df[sapply(dur_df, is.infinite)] <- NA


rm('monthly_data', 'presence', 'presence_long', 'df', 'data_prop', 'data_presence', 'ff_net', 'ffnet_reorder', 'result', 'subgroup_data')













#---------------------------------------------------------
#now we have the dates where all inds had day and night data for each group separately, going through each group, extracting the 2 week subsampled data for all inds (then will do it as though one was collared after this) - only for dates with day and night data


#making large for loop to run the same analyses across the different day subsampling
#assume they are observed at 10am UTC time (midday 12)
hour <- " 10:00:00"

# Make empty dataframe for the subsampled data
sub_df <- data.frame(group_i = character(),
                     focal_ind = character(),
                     focal_name = character(),
                     datetime = as.POSIXct(character()),
                     sub_com = I(list()),
                     stringsAsFactors = FALSE)


for(group_i in names(group_matrices)){
  
  group <- group_matrices[[group_i]]
  
  n_inds_i <- nrow(group$xs)
  
  if(n_inds_i == 1) next

  if(is.na(dur_df$dates_day_night[which(dur_df$group_name == group_i)])) next
  
  dates_str <- dur_df$all_dates[which(dur_df$group_name == group_i)]
  #convert string back to a date format to get the 2 week times when all that had data were tracked
  dates <- as.Date(strsplit(dates_str, ",\\s*")[[1]])
  
  # #go though each start date option
  # first_date <- dates[1:14]
  # 
  # for (date in seq_along(first_date)) {
  #   
  # 
    
  #dates_2week <- dates[seq(date, length(dates), 14)]
  sub_times <- as.POSIXct(paste0(dates, hour), tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  
  #get the xs and ys for the 10am subset each day
  sub_idx <- which(ts %in% sub_times)
  xs_10am_groupi <- group$xs[,sub_idx]
  ys_10am_groupi <- group$ys[,sub_idx]
  
  subgroup_data_group_i <- get_subgroup_data(xs_10am_groupi, ys_10am_groupi, R)
  
  ff_net_all <- matrix(NA, nrow = n_inds_i, ncol = n_inds_i)
  
  #going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
  for(i in 1:n_inds_i){
    for(j in 1:n_inds_i){
      
      #getting subgroup id for individual i and j
      sub_ids_i <- subgroup_data_group_i$ind_subgroup_membership[i,]
      sub_ids_j <- subgroup_data_group_i$ind_subgroup_membership[j,]
      
      #computing edge weight (fraction of time in same subgroup)
      ff_net_all[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
    }
  }
  
  diag(ff_net_all) <- NA
  visualize_lion_network(ff_net_all, group$id_info)



#comparing results when one ind in group is tracked
  for (t in 1:length(sub_idx)) {
    
    for (ind in 1:n_inds_i) {
  
    #get subgroup info for each individual at time t (subgroup id at that time for that ind)
    subgroup_ind <- subgroup_data_group_i$ind_subgroup_membership[ind, t]
  
    if (is.na(subgroup_ind)) next
  
      grp_members <- which(subgroup_data_group_i$ind_subgroup_membership[,t] == subgroup_ind)
      grp_members <- grp_members[grp_members != ind]
  
      # Append row to sub_df
      sub_df <- rbind(sub_df, data.frame(
        group_i = group_i,
        focal_ind = ind,
        focal_name = group$ids[ind],
        datetime = sub_times[t],
        sub_com = I(list(grp_members)),
        stringsAsFactors = FALSE
     ))
   }
  }
} 



# sub_df$sub_count <- sapply(sub_df$sub_com, length)+1 #to get total subgroup size, need to dd one to add the focal individual
# sub_df$date <- as.Date(sub_df$datetime)
# sub_df$month <- format(sub_df$datetime, "%m")
# 
# 
# #extract the mean and SE of subgroup size for each individual for each month for plotting
# summary_df <- sub_df %>%
#   mutate(month = lubridate::month(datetime, label = TRUE)) %>%
#   group_by(group_i, month) %>%
#   summarise(
#     mean_sub = mean(sub_count),
#     se_sub = sd(sub_count) / sqrt(n()),
#     .groups = "drop"
#   )
# 

# ggplot(summary_df, aes(x = month, y = mean_sub, color = factor(group_i))) +
#   geom_point(position = position_dodge(width = 0.5), size = 2) +
#   geom_errorbar(aes(ymin = mean_sub - se_sub, ymax = mean_sub + se_sub),
#                 width = 0.2,
#                 position = position_dodge(width = 0.5)) +
#   theme_classic() +
#   #facet_wrap(~group_i)+
#   labs(
#     title = "Monthly Mean Subgroup Sizes",
#     x = "Month",
#     y = "Mean Subgroup Size",
#     color = "Group"
#   ) +
#   theme(axis.text.x = element_text(hjust = 1))






#-------------------------------------------------------------

#----------MAKE PROP TIME TOGETHER PER FOCAL------------------

#-------------------------------------------------------------

#go through the sub_df dataframe for each focal, and recreate the matrix used for prop time together

#group_i <- "Etendeka"

#create list to store one matrix per focal individual
#for each subset start date every 14 days

subgroup_matrices <- list()

for(group_i in unique(sub_df$group_i)){
  
  if(is.na(dur_df$dates_day_night[dur_df$group_name == group_i])) next
  
  all_dates <- unique(sub_df$datetime[sub_df$group_i == group_i])
  first_date <- all_dates[1:14]
  
  if(is.null(subgroup_matrices[[group_i]])) subgroup_matrices[[group_i]] <- list()
  
  for (date in seq_along(first_date)) {
    
    date_i <- as.character(all_dates[date])
    dates_2week <- all_dates[seq(date, length(all_dates), 14)]
    
    subgroup_matrices[[group_i]][[date_i]] <- list()  
    
    sub_df_date_i <- sub_df[sub_df$group_i == group_i & sub_df$datetime %in% dates_2week, ]
    sub_df_date_i$time <- as_hms(sub_df_date_i$datetime)
    n_inds_i <- nrow(group_matrices[[group_i]]$xs) 
    names_inds_i <- group_matrices[[group_i]]$ids
    
    #loop over each focal individual
    for(i in unique(sub_df_date_i$focal_name)) {
      
      ind_i <- as.numeric(unique(sub_df_date_i$focal_ind[which(sub_df_date_i$focal_name == i)]))
      sub_df_focal_i <- sub_df_date_i[sub_df_date_i$focal_name == i, ]
      n_timepoints <- nrow(sub_df_focal_i)
      
      # Initialize matrix
      matrix_i <- matrix(0, 
                         nrow = n_inds_i, 
                         ncol = n_timepoints,
                         dimnames = list(names_inds_i, paste0("t", 1:n_timepoints )))
      
      for(j in 1:n_timepoints) {
        matrix_i[ind_i, j] <- 1
        sub_comp <- unlist(sub_df_focal_i$sub_com[j])
        if (length(sub_comp) > 0 && !all(is.na(sub_comp))) {
          matrix_i[sub_comp, j] <- 1
        }
      }
      
      # Store in nested list
      subgroup_matrices[[group_i]][[date_i]][[as.character(i)]] <- matrix_i 
    }
  }
}


#--------------------------------------------------------------------------
  
  
  
# Create a list to store focal matrices
focal_matrices <- list()

#group_i <- "Aub"

for(group_i in names(subgroup_matrices)){
  
  #get group for correct plotting of IDs in matrix later
  group <- group_matrices[[group_i]]
  
  sub_df_i <- sub_df[sub_df$group_i == group_i,]
  sub_df_i$time <- as_hms(sub_df_i$datetime)
  n_inds_i <- nrow(group$xs) #cant do length as sometimes there is an ind missing data
  names_inds_i <- group$ids
  
  subgroup_matrices_i <- subgroup_matrices[[group_i]]
  
  if(n_inds_i == 1) next

# Matrix to store results
ff_net_matrix <- matrix(NA, nrow = n_inds_i, ncol = n_inds_i)

#go through each dyad from the focal inds perspective to get prop time all inds are with one another
for(focal in names(subgroup_matrices_i)){
  
  #get iteration number for bar plot code when removing focal ind in calc
  iter_num <- match(focal, group$id_info$name)
  
  
  focal_i_matrix <- subgroup_matrices_i[[focal]]
  focal_obs <- !is.na(focal_i_matrix[focal, ])
  
  for (i in 1:n_inds_i) {
    for (j in 1:n_inds_i) {
      
      # Timepoints where all 3 (focal, i, j) are observed
      valid_t <- focal_obs & !is.na(focal_i_matrix[i, ]) & !is.na(focal_i_matrix[j, ])
      
      if (sum(valid_t) > 0) {
        # Proportion of time i and j were in same subgroup, when with focal
        ff_net_matrix[i, j] <- mean(focal_i_matrix[i, valid_t] == focal_i_matrix[j, valid_t])
      } else {
        ff_net_matrix[i, j] <- NA
      }
    }
  }
  
  diag(ff_net_matrix) <- NA
  #new_order <- 1:n_inds_i
  #ffnet_reorder <- ff_net_matrix[new_order, new_order]
  
  # Store in list, named by focal individual
  focal_matrices[[group_i]][[focal]] <- ff_net_matrix
  
  png(height = 400, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_network_focal_', focal, '-', group_i, ' ', days, '-days.png'))
  
  par(mfrow = c(1, 2))
  visualize_lion_network(ff_net_matrix, group$id_info)
  mtext(group_i,  cex = 2)
  
  # Remove the focal individual from comparison
  values <- ff_net_matrix[iter_num, -iter_num]
  names <- group$id_info$name[-iter_num]
  
  bar_mids <- barplot(values, col = "slateblue",
                      main = paste("Proportion of Time in Subgroup with", focal),
                      ylim = c(0, 1), xaxt = "n", cex.axis = 1.5)
  
  # Add correct axis labels at the bar midpoints
  axis(1, at = bar_mids, labels = names, las = 2, cex.axis = 1.5)
  
  dev.off()
}
}





  

#----------------------------------------------------------------------------

#Full group association matrices from total dataset (to compare with subsampled matrices)

#----------------------------------------------------------------------------



#compare which focals have closest result to the total dataset

for(group_i in names(group_matrices)){
  
  group <- group_matrices[[group_i]]
  n_inds_i <- nrow(group$xs)
  
  if(n_inds_i == 1) next
  subgroup_data <- get_subgroup_data(group$xs, group$ys, R)
  ff_net_all <- matrix(NA, nrow = n_inds_i, ncol = n_inds_i)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds_i){
  for(j in 1:n_inds_i){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net_all[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net_all) <- NA
total_matrix <- ff_net_all

#--------compare rank of sub sampled matrices to full matrix-----------------

# Function to extract vectorized upper triangle of a matrix (excluding diagonal)
get_upper_triangle <- function(mat) {
  mat[lower.tri(mat)] <- NA  # keep upper triangle only
  return(as.vector(mat[upper.tri(mat)]))  # extract as vector
}

# Vector from real data
real_vector <- get_upper_triangle(total_matrix)

# Initialize result vector for each ind
ids <- group_matrices[[group_i]]$ids
spearman_scores <- setNames(numeric(length(ids)), ids)
spearman_scores[] <- NA


#going though all inds in each group to see whether is focal data for the subset matrices, if not, replace with NA
for (ind in group_matrices[[group_i]]$ids) {
  
  # check if focal matrix exists for this individual
  if (is.null(focal_matrices[[group_i]][[ind]])) {
    
    # create empty NA matrix
    focal_matrices[[group_i]][[ind]] <- matrix(
      NA,
      nrow = n_inds_i,
      ncol = n_inds_i,
      dimnames = list(group_matrices[[group_i]]$ids,
                      group_matrices[[group_i]]$ids)
    )
    
    focal_vec <- rep(NA, n_inds_i * (n_inds_i - 1) / 2)
    
  } else {
  
  
  focal_vec <- get_upper_triangle(focal_matrices[[group_i]][[ind]])
  }
  
  # Ensure comparison only for dyads with valid values in both
  valid_idx <- !is.na(real_vector) & !is.na(focal_vec)
  
  # Spearman correlation
  if (sum(valid_idx) >= 2) {
    spearman_scores[ind] <- cor(real_vector[valid_idx], focal_vec[valid_idx], method = "spearman")
  
  } else {
    spearman_scores[ind] <- NA  # not enough valid points
  }
}

spearman_long <- rbind(
  spearman_long,
  data.frame(
    start_date = date,
    individual = names(spearman_scores),
    spearman   = spearman_scores
  ))

if (nrow(spearman_long) == 0) {
  message("No Spearman results for group ", group_i, " (", days, "-day)")
  next
}

# summarise spearman across start dates
spearman_summary <- aggregate(
  spearman ~ individual,
  data = spearman_long,
  FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd   = if (sum(!is.na(x)) >= 2) sd(x, na.rm = TRUE) else NA)
  }
)

# unpack the matrix column
spearman_summary$mean_spearman <- spearman_summary$spearman[, "mean"]
spearman_summary$sd_spearman   <- spearman_summary$spearman[, "sd"]
spearman_summary$spearman <- NULL

no_data_inds <- spearman_sd$individual[is.na(spearman_sd$sd_spearman)]

if (length(no_data_inds) > 0) {
  message(
    "Individuals with insufficient data for SD (",
    group_i, "): ",
    paste(no_data_inds, collapse = ", ")
  )
}

plot_df <- spearman_sd[!is.na(spearman_sd$sd_spearman), ]

if (nrow(plot_df) == 0) next 

#colnames(spearman_sd)[2] <- "sd_spearman"



png(height = 450, width = 1000, units = 'px', filename = paste0(plot_dir, 'spearman_scores for ', group_i, " - ", days,  ' day subsampling.png'))

par(mfrow = c(1, 2))

visualize_lion_network(total_matrix, group$id_info)

par(mar = c(15,5,4,5))
#bar_mids <- barplot(spearman_scores,
                    # names.arg = paste("Focal", seq_along(focal_matrices[[group_i]])),
                    # col = "slateblue",
                    # ylim = c(0, 1),
                    # main = paste0("Spearman rank correlation for " ,group_i, " group ", days, "-day interval"),
                    # ylab = "Rank", xaxt = "n")

# Add correct axis labels at the bar midpoints
#axis(1, at = bar_mids, labels = group$id_info$name, las = 2, cex.axis = 1)


bar_mids <- barplot(
  spearman_summary$mean_spearman,
  names.arg = spearman_summary$individual,
  ylim = c(-1, 1),
  col = "steelblue",
  ylab = "Mean Spearman rank (± SD)",
  main = paste0(
    "Mean Spearman rank vs full network\n",
    group_i, " (", days, "-day subsampling)"
  ),
  las = 2
)

# add SD error bars
with(spearman_summary, {
  ok <- !is.na(sd_spearman)
  
  arrows(
    x0 = bar_mids[ok],
    y0 = mean_spearman[ok] - sd_spearman[ok],
    x1 = bar_mids[ok],
    y1 = mean_spearman[ok] + sd_spearman[ok],
    angle = 90,
    code = 3,
    length = 0.05,
    lwd = 2
  )
})

abline(h = 0, lty = 2, col = "grey40")

dev.off()
}
}


}
#------------------------------------------------





#------------------------------------------------

#quantify the variance of subgroup size and number of subgroups in start date
#boxplot(sub_df$sub_count ~ as.numeric(sub_df$month))

#ggplot(sub_df, aes(month, sub_count, fill = as.factor(focal_ind)))+
#  geom_boxplot()






































