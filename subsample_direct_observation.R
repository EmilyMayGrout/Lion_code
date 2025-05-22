#This script will subsample the 15-min lion GPS data to imitate the direct behavioural observations field biologists collect to determine how good this method is at describing their subgrouping 
#Every 2 weeks, people go out between 10am and midday to record the location of a radio collared individual and the individuals within observation distance (R) as being in the subgroup (to get average subgroup size over time)


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
library(fields)
library(asnipe)

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


#Remove time Simba died

# Find the index of the last non-NA value
last_valid_index <- max(which(!is.na(xs[3,])))

# Get the corresponding timestamp
ts[last_valid_index]

#insert NA to after Simba died
xs[3, last_valid_index:ncol(xs)] <- NA
ys[3, last_valid_index:ncol(xs)] <- NA



#----SUBSAMPLING-------------------


# Make empty dataframe for the subsampled data
sub_df <- data.frame(focal_ind = character(),
                     datetime = as.POSIXct(character()),
                     sub_com = I(list()),   # Use list column for storing vectors
                     stringsAsFactors = FALSE)

# Observation times (as_hms is from hms or lubridate)
obs_time <- as_hms(c("10:00:00", "10:30:00", "11:00:00", "11:30:00", "12:00:00"))

# Dates
dates <- unique(as_date(ts))
first_date <- dates[1:29]

# Loop each 2 week observation (for different start dates), then loop through each ind as the focal to get the number of individuals in their subgroup for that time 
for (j in obs_time) {
  j_time <- as_hms(j)
  
  for (i in seq_along(first_date)) {
    two_week_subsampling <- dates[seq(i, length(dates), 14)]
    
    for (p in seq_along(two_week_subsampling)) {
      date_time_jip <- as.POSIXct(paste(two_week_subsampling[p], j_time), 
                                  format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
      
      ts_jip <- which(ts == date_time_jip)
      if (length(ts_jip) == 0) next  # Skip if datetime doesn't exist in ts
      
      for (ind in 1:n_inds) {
        subgroup_jip <- subgroup_data$ind_subgroup_membership[ind, ts_jip]
        
        if (is.na(subgroup_jip)) next
        
        grp_members_jip <- which(subgroup_data$ind_subgroup_membership[, ts_jip] == subgroup_jip)
        grp_members_jip <- grp_members_jip[grp_members_jip != ind]
        
        # Append row to sub_df
        sub_df <- rbind(sub_df, data.frame(
          focal_ind = ind,
          datetime = date_time_jip,
          sub_com = I(list(grp_members_jip)), 
          date_sample = i, # Store list
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

sub_df$sub_count <- sapply(sub_df$sub_com, length)+1 #to get total subgroup size, need to dd one to add the focal individual
sub_df$date <- as.Date(sub_df$datetime)
sub_df$time <- as_hms(sub_df$datetime)
sub_df$month <- format(sub_df$datetime, "%m")

#look at the distribution of subgroup sizes for one individual
#just look at one of the dates in the dateshift 
sub_df_filt <- sub_df[sub_df$date_sample == 4,]
#sub_df_filt <- sub_df[sub_df_filt$time == as_hms("10:00:00"),]


#extract the mean and SE of subgroup size for each individual for each month for plotting
summary_df <- sub_df_filt %>%
  mutate(month = lubridate::month(datetime, label = TRUE)) %>%
  group_by(focal_ind, month, time) %>%
  summarise(
    mean_sub = mean(sub_count),
    se_sub = sd(sub_count) / sqrt(n()),
    .groups = "drop"
  )


gg <- ggplot(summary_df, aes(x = month, y = mean_sub, color = factor(focal_ind), group = focal_ind)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_sub - se_sub, ymax = mean_sub + se_sub),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  theme_classic() +
  facet_wrap(~time)+
  labs(
    title = "Monthly Mean Subgroup Sizes by Individual",
    x = "Month",
    y = "Mean Subgroup Size",
    color = "Focal Individual"
  ) +
  theme(axis.text.x = element_text(hjust = 1))

ggsave(filename = paste0(plot_dir, '2_week_sampling_subgroup_size.png'), plot = gg, width = 10, height = 7, dpi = 300)




#run anova between the different start dates to see how this affects the results

# Prepare storage dataframe
anova_results <- data.frame(
  date_sample = character(),
  variable = character(),
  df = numeric(),
  sum_sq = numeric(),
  mean_sq = numeric(),
  f_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each unique date_sample
for(i in unique(sub_df$date_sample)) {
  
  sub_df_filt <- sub_df[sub_df$date_sample == i, ]
  
  # Check if enough data to fit model
  if(nrow(sub_df_filt) > 2 && length(unique(sub_df_filt$focal_ind)) > 1) {
    
    model <- aov(sub_count ~ focal_ind + time, data = sub_df_filt)
    model_summary <- summary(model)[[1]]
    
    # Loop over rows in summary (one per factor)
    for(j in 1:nrow(model_summary)) {
      anova_results <- rbind(anova_results, data.frame(
        date_sample = i,
        variable = rownames(model_summary)[j],
        df = model_summary$Df[j],
        sum_sq = model_summary$`Sum Sq`[j],
        mean_sq = model_summary$`Mean Sq`[j],
        f_value = model_summary$`F value`[j],
        p_value = model_summary$`Pr(>F)`[j]
      ))
    }
  }
}


anova_results_filt <- anova_results[anova_results$variable != "Residuals  ",]

#Show how the statistical significance of each variable changes over time

#Plot of p-values over date_sample for each variable:
gg <- ggplot(anova_results_filt, 
       aes(x = as.Date(date_sample), y = p_value, color = variable)) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "slateblue3")+
  geom_point() +
  labs(title = "ANOVA p-values over time",
       y = "p-value", x = "Date") +
  theme_classic()

ggsave(filename = paste0(plot_dir, 'predictors_significance_for_subgroup_size.png'), plot = gg, width = 10, height = 6, dpi = 300)

#this plot shows the significance of sub count for the focal ind depending on the date the samples are taken
#the time of sampling only occassionally has a significant affect on the subgroup size
#the focal individual significantly affects the sub-count for most start dates, but this varies depending on when the start date is - sometimes it doesn't matter on when the start date is.....
#



anova_model <- aov(sub_count ~ focal_ind + time + date_sample, data = sub_df)
summary(anova_model)

#----------MAKE PROP TIME TOGETHER PER FOCAL------------------
#
#first need to go through the sub_df dataframe for each focal, and recreate the matrix used for prop time together

#filter to just the 10am data collection for now
sub_df$time <- as_hms(sub_df$datetime)
sub_df_10am <- sub_df[sub_df$time == as_hms("10:00:00"),]
sub_df_10am_datesamp1 <- sub_df_10am[sub_df_10am$date_sample == 1,]


# Create a list to store one matrix per focal individual
subgroup_matrices <- list()

# Loop over each focal individual
for(i in unique(sub_df_10am_datesamp1$focal_ind)) {
  
  sub_df_focal_i <- sub_df_10am_datesamp1[sub_df_10am_datesamp1$focal_ind == i, ]
  n_timepoints <- nrow(sub_df_focal_i)
  
    # Initialize matrix: rows = individuals, columns = timepoints
  matrix_i <- matrix(0, nrow = n_inds, ncol = n_timepoints,
                     dimnames = list(1:n_inds, paste0("t", 1:n_timepoints)))
  
  for(j in 1:n_timepoints) {
    
    # Always mark the focal individual present (value = 1)
    matrix_i[i, j] <- 1
    
    # Get individuals in subgroup at time j
    sub_comp <- unlist(sub_df_focal_i$sub_com[j])
    
    if (length(sub_comp) > 0 && !all(is.na(sub_comp))) {
      matrix_i[sub_comp, j] <- 1
    }
    # Others remain NA
  }
  
  
  #because Simba died, need to put NAs into the matrices so these times are not caluclated as proportion of time apart
  if(i != 3){
    
    #find columns after Simba's death and add NAs to Simba's row
    matrix_i[3, which(sub_df_focal_i$datetime > ts[last_valid_index]) ] <- NA
    
  }
  
  # Store in list, named by focal individual
  subgroup_matrices[[as.character(i)]] <- matrix_i
}


#each matrix has the individuals involved for that focal, ignoring the subgroup ids for inds not in the same subgroup


# Create a list to store focal matrices
focal_matrices <- list()

# Matrix to store results
ff_net_matrix <- matrix(NA, nrow = n_inds, ncol = n_inds)

# Timepoints where focal (ID 1) is observed

for(p in 1:n_inds){
  
focal_i_matrix <- subgroup_matrices[[p]]
focal_obs <- !is.na(focal_i_matrix[p, ])

for (i in 1:n_inds) {
  for (j in 1:n_inds) {
    
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
new_order <- c(1:3,5,4,6)
ffnet_reorder <- ff_net_matrix[new_order, new_order]

# Store in list, named by focal individual
focal_matrices[[as.character(p)]] <- ffnet_reorder

name <- lion_ids$name[p]

png(height = 400, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_network_focal_', name, '.png'))

par(mfrow = c(1, 2))
visualize_network_matrix_trago(ffnet_reorder, lion_ids[new_order,])


# Remove the focal individual from comparison
values <- ff_net_matrix[p, -p]
names <- lion_ids$name[-p]

bar_mids <- barplot(values, col = "slateblue",
                    main = paste("Proportion of Time in Subgroup with", name),
                    ylim = c(0, 1), xaxt = "n", cex.axis = 1.5)

# Add correct axis labels at the bar midpoints
axis(1, at = bar_mids, labels = names, las = 2, cex.axis = 1.5)

dev.off()
}





#compare which focals have closest result to the total dataset

R = 100
subgroup_data <- get_subgroup_data(xs, ys, R)

ff_net_all <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net_all[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net_all) <- NA
new_order <- c(1:3,5,4,6)
ffnet_reorder_all <- ff_net_all[new_order, new_order]

visualize_network_matrix_trago(ffnet_reorder_all, lion_ids[new_order,])


#--------compare rank of sub sampled matrices to full matrix-----------------

real_matrix <- ffnet_reorder_all

# Function to extract vectorized upper triangle of a matrix (excluding diagonal)
get_upper_triangle <- function(mat) {
  mat[lower.tri(mat)] <- NA  # keep upper triangle only
  return(as.vector(mat[upper.tri(mat)]))  # extract as vector
}

# Vector from real data
real_vector <- get_upper_triangle(real_matrix)

# Initialize result vector
spearman_scores <- numeric(length(focal_matrices))

for (i in seq_along(focal_matrices)) {
  focal_vec <- get_upper_triangle(focal_matrices[[i]])
  
  # Ensure comparison only for dyads with valid values in both
  valid_idx <- !is.na(real_vector) & !is.na(focal_vec)
  
  # Spearman correlation
  if (sum(valid_idx) >= 10) {
    spearman_scores[i] <- cor(real_vector[valid_idx], focal_vec[valid_idx], method = "spearman")
  } else {
    spearman_scores[i] <- NA  # not enough valid points
  }
}

png(height = 400, width = 500, units = 'px', filename = paste0(plot_dir,'spearman_scores.png'))

par(mar = c(8,5,4,5))
bar_mids <- barplot(spearman_scores,
        names.arg = paste("Focal", seq_along(focal_matrices)),
        col = "slateblue",
        ylim = c(0, 1),
        main = "Spearman rank correlation between focal's and full group matrix",
        ylab = "Rank", xaxt = "n")

# Add correct axis labels at the bar midpoints
axis(1, at = bar_mids, labels = lion_ids$name[new_order], las = 2, cex.axis = 1)

dev.off()







