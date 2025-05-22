#helpful for comparing differences between dataframes
#anti_join(data.frame(t = splits_df2$t), data.frame(t = splits_df$t), by = "t")



detect_merges_across_sample_rates <- function(subsampled_results, include_singletons = TRUE) {
  merge_list_samplingrates <- list()
  element_names <- names(subsampled_results)
  
  for (name in element_names) {
    samplerate <- subsampled_results[[name]]
    merge <- c()
    n_times_cut <- ncol(samplerate$matrix)
    
    for (t in 2:n_times_cut) {
      merge_group <- samplerate$matrix[, t]
      subgroups_previously <- samplerate$matrix[, t - 1]
      
      if (sum(!is.na(merge_group)) == 0 | sum(!is.na(subgroups_previously)) == 0) next
      
      merge_group[is.na(subgroups_previously)] <- NA
      subgroups_previously[is.na(merge_group)] <- NA
      
      n_merge_group <- length(unique(merge_group[!is.na(merge_group)]))
      n_subgroups_previously <- length(unique(subgroups_previously[!is.na(subgroups_previously)]))
      
      singletons_now <- sum(table(merge_group) == 1)
      singletons_later <- sum(table(subgroups_previously) == 1)
      
      if (!include_singletons) {
        if (n_merge_group == 1 &&
            n_subgroups_previously > 1 &&
            singletons_later == 0) {
          merge <- c(merge, t)
        }
        if (n_merge_group > 1 &&
            (singletons_now + 1) == n_merge_group &&
            n_subgroups_previously > n_merge_group &&
            singletons_now == singletons_later) {
          merge <- c(merge, t)
        }
      } else {
        if (n_merge_group < n_subgroups_previously) {
          merge <- c(merge, t)
        }
      }
    }
    
    if (length(merge) == 0) {
      merge_list_samplingrates[[name]] <- data.frame(
        t = integer(),
        merge_group = I(list()),
        sub1 = I(list()),
        sub2 = I(list()),
        sub3 = I(list()),
        sub4 = I(list()),
        sub5 = I(list()),
        subsample_indices = I(list()),
        datetime = I(list()),
        n_merge = integer(),
        n_sub1 = integer(),
        n_sub2 = integer(),
        n_sub3 = integer()
      )
      next
    }
    
    merge_df <- data.frame(t = merge,
                           merge_group = NA, sub1 = NA, sub2 = NA,
                           sub3 = NA, sub4 = NA, sub5 = NA,
                           subsample_indices = NA)
    
    for (i in 1:nrow(merge_df)) {
      t <- merge_df$t[i]
      merge_group <- samplerate$matrix[, t]
      subgroups_previously <- samplerate$matrix[, t - 1]
      
      merge_group[is.na(subgroups_previously)] <- NA
      subgroups_previously[is.na(merge_group)] <- NA
      
      changed_indices <- which(merge_group != subgroups_previously)
      group_id_changers <- merge_group[changed_indices]
      prev_id_changers <- subgroups_previously[changed_indices]
      
      merge_df$merge_group[i] <- list(which(merge_group == group_id_changers[1]))
      group_ids_before <- unique(prev_id_changers)
      
      for (j in 1:length(group_ids_before)) {
        group_id <- group_ids_before[j]
        inds_in_group <- which(subgroups_previously == group_id)
        orig_inds_in_group <- intersect(inds_in_group, changed_indices)
        
        merge_df[[paste0("sub", j)]][i] <- list(orig_inds_in_group)
      }
      
      merge_df$subsample_indices[i] <- list(samplerate$subsample_indices[t])
      merge_df$datetime[i] <- list(samplerate$subsampled_ts[t + 1])
    }
    
    merge_df$n_merge <- sapply(merge_df$merge_group, function(x) sum(!is.na(x)))
    merge_df$n_sub1 <- sapply(merge_df$sub1, function(x) sum(!is.na(x)))
    merge_df$n_sub2 <- sapply(merge_df$sub2, function(x) sum(!is.na(x)))
    merge_df$n_sub3 <- sapply(merge_df$sub3, function(x) sum(!is.na(x)))
    
    merge_list_samplingrates[[name]] <- merge_df
  }
  
  return(merge_list_samplingrates)
}


#------------------------------------------------------------------------
#extracting merges from a single matrix
detect_merges <- function(ind_subgroup_membership) {
  
  n_times <- ncol(subgroup_data$ind_subgroup_membership)
  merges <- c()
  
  for (t in 1:(n_times - 1)) {
    
    now <- subgroup_data$ind_subgroup_membership[, t]
    before <- subgroup_data$ind_subgroup_membership[, t - 1]
    
    # Keep only individuals with data in both time steps
    valid_inds <- !is.na(now) & !is.na(before)
    if (sum(valid_inds) == 0) next
    
    now_valid <- now[valid_inds]
    before_valid <- before[valid_inds]
    
    # Find the modal group at time t-1
    modal_group <- as.numeric(names(which.max(table(before_valid))))
    modal_inds <- which(valid_inds)[before_valid == modal_group]
    
    # Count subgroups at t and t+1
    n_subgroups_now <- length(unique(now_valid))
    n_subgroups_before <- length(unique(before_valid))
    
    # Check what groups these individuals came from at t-1
    before_groups <- before[modal_inds]
    unique_before_groups <- unique(before_groups[!is.na(before_groups)])
    
    if (length(unique_before_groups) > 1 | n_subgroups_before > n_subgroups_now) {
      merges <- c(merges, t)
    }
  }
  
  
  # Create a data frame of merges, with contributing groups and merged group
  merges_df <- data.frame(t = merges, merged_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)
  
  for(i in 1:nrow(merges_df)){
    t <- merges_df$t[i]
    subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
    subgroups_before <- subgroup_data$ind_subgroup_membership[,t-1]
    
    # Ensure NAs are synced
    subgroups_now[is.na(subgroups_before)] <- NA
    subgroups_before[is.na(subgroups_now)] <- NA
    
    
    # Identify which groups at t-1 merged at time t
    merge_groups <- sapply(unique(subgroups_now), function(h) {
      inds <- which(subgroups_now == h)
      length(unique(subgroups_before[inds])) > 1
    })
    
    
    # Get the group IDs in subgroups_now that are results of merges
    merging_groups <- unique(subgroups_now)[merge_groups]
    
    # Get all individuals involved in any merging group
    merged_inds <- which(subgroups_now %in% merging_groups)
    
    # Find their prior groupings
    prev_group_ids <- unique(subgroups_before[merged_inds])
    
    # Store full merged membership
    merges_df$merged_group[i] <- list(merged_inds)
    
    # Store each contributing subgroup
    for(j in 1:length(prev_group_ids)){
      
      prev_group_ids_j <- which(subgroups_before == prev_group_ids[j])
      
      if(j==1){ 
        merges_df$sub1[i] <- list(prev_group_ids_j) }
      if(j==2){ 
        merges_df$sub2[i] <- list(prev_group_ids_j) }
      if(j==3){ 
        merges_df$sub3[i] <- list(prev_group_ids_j) }
      if(j==4){ 
        merges_df$sub4[i] <- list(prev_group_ids_j) }
      if(j==5){ 
        merges_df$sub5[i] <- list(prev_group_ids_j) 
      }
    }
  }
  # Add sizes of each group
  merges_df$n_merged <- sapply(merges_df$merged_group, function(x){sum(!is.na(x))})
  merges_df$n_sub1  <- sapply(merges_df$sub1, function(x){sum(!is.na(x))})
  merges_df$n_sub2  <- sapply(merges_df$sub2, function(x){sum(!is.na(x))})
  merges_df$n_sub3  <- sapply(merges_df$sub3, function(x){sum(!is.na(x))})
  merges_df$n_sub4  <- sapply(merges_df$sub4, function(x){sum(!is.na(x))})
  merges_df$n_sub5  <- sapply(merges_df$sub5, function(x){sum(!is.na(x))})
  
  return(merges_df)
  
}

detect_splits_across_sample_rates <- function(subsampled_results, include_singletons = TRUE) {
  split_list_samplingrates <- list()
  element_names <- names(subsampled_results)
  
  for (name in element_names) {
    samplerate <- subsampled_results[[name]]
    splits <- c()
    n_times_cut <- ncol(samplerate$matrix)
    
    for (t in 1:(n_times_cut - 1)) {
      subgroups_now <- samplerate$matrix[, t]
      subgroups_later <- samplerate$matrix[, t + 1]
      
      if (sum(!is.na(subgroups_now)) == 0 | sum(!is.na(subgroups_later)) == 0) next
      
      subgroups_now[is.na(subgroups_later)] <- NA
      subgroups_later[is.na(subgroups_now)] <- NA
      
      n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
      n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
      
      singletons_now <- sum(table(subgroups_now) == 1)
      singletons_later <- sum(table(subgroups_later) == 1)
      
      # Check what groups these individuals belong to at t+1
      later_groups <- later[modal_inds]
      unique_later_groups <- unique(later_groups[!is.na(later_groups)])
      
      if (!include_singletons) {
        if ((n_subgroups_now == 1) &&
            n_subgroups_later > n_subgroups_now &&
            singletons_later == 0) {
          splits <- c(splits, t)
        }
        if (n_subgroups_now > 1 &&
            ((singletons_now + 1) == n_subgroups_now) &&
            n_subgroups_later > n_subgroups_now &&
            singletons_now == singletons_later) {
          splits <- c(splits, t)
        }
      } else {
        if (n_subgroups_later > n_subgroups_now) {
          splits <- c(splits, t)
       
       }
      }
     }
    
    if (length(splits) == 0) {
      split_list_samplingrates[[name]] <- data.frame(
        t = integer(),
        split_group = I(list()),
        sub1 = I(list()),
        sub2 = I(list()),
        sub3 = I(list()),
        sub4 = I(list()),
        sub5 = I(list()),
        subsample_indices = I(list()),
        datetime = I(list()),
        n_split = integer(),
        n_sub1 = integer(),
        n_sub2 = integer(),
        n_sub3 = integer()
      )
      next
    }
    
    split_df <- data.frame(t = splits,
                           split_group = NA, sub1 = NA, sub2 = NA,
                           sub3 = NA, sub4 = NA, sub5 = NA,
                           subsample_indices = NA)
    
    for (i in 1:nrow(split_df)) {
      t <- split_df$t[i]
      subgroups_now <- samplerate$matrix[, t]
      subgroups_later <- samplerate$matrix[, t + 1]
      
      subgroups_now[is.na(subgroups_later)] <- NA
      subgroups_later[is.na(subgroups_now)] <- NA
      
      changed_indices <- which(subgroups_now != subgroups_later)
      group_id_changers <- subgroups_now[changed_indices]
      new_id_changers <- subgroups_later[changed_indices]
      
      split_df$split_group[i] <- list(which(subgroups_now == group_id_changers[1]))
      group_ids_after <- unique(new_id_changers)
      
      for (j in 1:length(group_ids_after)) {
        group_id <- group_ids_after[j]
        inds_in_group <- which(subgroups_later == group_id)
        orig_inds_in_group <- intersect(inds_in_group, changed_indices)
        
        split_df[[paste0("sub", j)]][i] <- list(orig_inds_in_group)
      }
      
      split_df$subsample_indices[i] <- list(samplerate$subsample_indices[t])
      split_df$datetime[i] <- list(samplerate$subsampled_ts[t + 1])
    }
    
    split_df$n_split <- sapply(split_df$split_group, function(x) sum(!is.na(x)))
    split_df$n_sub1 <- sapply(split_df$sub1, function(x) sum(!is.na(x)))
    split_df$n_sub2 <- sapply(split_df$sub2, function(x) sum(!is.na(x)))
    split_df$n_sub3 <- sapply(split_df$sub3, function(x) sum(!is.na(x)))
    
    split_list_samplingrates[[name]] <- split_df
  }
  
  return(split_list_samplingrates)
}


#---------------------------------------------------------------------------------
#function for one matrix:
detect_splits <- function(ind_subgroup_membership) {
  
  n_times <- ncol(ind_subgroup_membership)
  splits <- c()
  
  for (t in 1:(n_times - 1)) {
    now <- subgroup_data$ind_subgroup_membership[, t]
    later <- subgroup_data$ind_subgroup_membership[, t + 1]
    
    # Keep only individuals with data in both time steps
    valid_inds <- !is.na(now) & !is.na(later)
    if (sum(valid_inds) == 0) next
    
    now_valid <- now[valid_inds]
    later_valid <- later[valid_inds]
    
    # Find the modal group in 'now'
    modal_group <- as.numeric(names(which.max(table(now_valid)))) #finds the subgroup ID at time t that has the most individuals
    modal_inds <- which(valid_inds)[now_valid == modal_group]
    
    #if its not the modal group that split
    #get number of subgroups now and in next time step (later)
    n_subgroups_now <- length(unique(now_valid))
    n_subgroups_later <- length(unique(later_valid))
    
    # Check what groups these individuals belong to at t+1
    later_groups <- later[modal_inds]
    unique_later_groups <- unique(later_groups[!is.na(later_groups)])
    
    if (length(unique_later_groups) > 1 | n_subgroups_later > n_subgroups_now) {
      splits <- c(splits, t)
    }
  }
  
  
  #make a data frame of splits, with original group and subgroups
  splits_df <- data.frame(t = splits, orig_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)
  for(i in 1:nrow(splits_df)){
    t <- splits_df$t[i]
    subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
    subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
    
    #transfer NAs from now to later and later to now so erroneous splits are not detected
    subgroups_now[which(is.na(subgroups_later))] <- NA
    subgroups_later[which(is.na(subgroups_now))] <- NA
    
    split_groups <- sapply(unique(subgroups_now), function(g) {
      inds <- which(subgroups_now == g)
      length(unique(subgroups_later[inds])) > 1
    })
    
    # Get the group IDs in `subgroups_now` that split
    splitting_group <- unique(subgroups_now)[split_groups]
    
    # Find individuals in the original group
    orig_change_inds <- which(subgroups_now %in% splitting_group)
    
    #find the groups where the original members went
    group_ids_later <- unique(subgroups_later[orig_change_inds])
    
    #store original group membership in data frame
    splits_df$orig_group[i] <- list(orig_change_inds)
    
    
    for(j in 1:length(group_ids_later)){
      
      group_id <- group_ids_later[j]
      inds_in_group <- which(subgroups_later == group_id)
      
      orig_inds_in_group <- intersect(inds_in_group, orig_change_inds) #only count the original group members 
      
      #put lists into a data frame
      if(j==1){
        splits_df$sub1[i] <- list(orig_inds_in_group)
      }
      if(j==2){
        splits_df$sub2[i] <- list(orig_inds_in_group)
      }
      if(j==3){
        splits_df$sub3[i] <- list(orig_inds_in_group)
      }
      if(j==4){
        splits_df$sub4[i] <- list(orig_inds_in_group)
      }
      if(j==5){
        splits_df$sub5[i] <- list(orig_inds_in_group)
      }
    }
  }
  
  #number in each subgroup
  splits_df$n_orig <- sapply(splits_df$orig_group, function(x){return(sum(!is.na(x)))})
  splits_df$n_sub1 <- sapply(splits_df$sub1, function(x){return(sum(!is.na(x)))})
  splits_df$n_sub2 <- sapply(splits_df$sub2, function(x){return(sum(!is.na(x)))})
  splits_df$n_sub3 <- sapply(splits_df$sub3, function(x){return(sum(!is.na(x)))})
  splits_df$n_sub4 <- sapply(splits_df$sub4, function(x){return(sum(!is.na(x)))})
  splits_df$n_sub5 <- sapply(splits_df$sub5, function(x){return(sum(!is.na(x)))})
  
  return(splits_df)
}


#----------------------------------------------------------------







#split_list_samplingrates2 <- detect_splits(subsampled_results, include_singletons = T)


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

# Function to extract the number of 15-min steps skipped
get_steps_skipped <- function(name) {
  if (grepl("days", name)) {
    # Extract number of days
    days <- as.numeric(gsub(".*_(\\d+)days.*", "\\1", name))
    return(days * 24 * 4)  # 4 samples per hour * 24 hours
  } else {
    # Extract HH_MM_SS format
    time_str <- sub("subsampled_matrix_", "", name)
    parts <- unlist(strsplit(time_str, "_"))
    hours <- as.numeric(parts[1])
    minutes <- as.numeric(parts[2])
    seconds <- as.numeric(parts[3])
    total_minutes <- hours * 60 + minutes + seconds / 60
    return(total_minutes / 15)
  }
}


subsample_multiple_intervals <- function(ind_subgroup_membership, ts, intervals) {
  # Create a list to store the results
  subsampled_results <- list()
  
  # Loop through each interval and subsample
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
  
  return(subsampled_results)
}


subsample_multiple_day_intervals <- function(ind_subgroup_membership, ts, day_intervals) {
  # Create a list to store the results
  subsampled_results <- list()
  
  # Loop through each interval and subsample
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
  
  return(subsampled_results)
}



compute_merge_split_time_differences <- function(split_list_samplingrates, merge_list_samplingrates) {
  
  time_diff_list <- list()
  element_names <- names(split_list_samplingrates)
  
  for (n in element_names) {
    
    splits_df <- split_list_samplingrates[[n]]
    merge_df <- merge_list_samplingrates[[n]]
    
    if (length(merge_df) == 0 || length(splits_df) == 0) next
    
    # Label event types
    merge_df$event <- "merge"
    splits_df$event <- "split"
    
    # Extract relevant columns
    merge_df_subset <- merge_df[, c("subsample_indices", "event")]
    splits_df_subset <- splits_df[, c("subsample_indices", "event")]
    
    # Combine and sort
    time_diff <- rbind(splits_df_subset, merge_df_subset) %>%
      mutate(subsample_indices = as.numeric(subsample_indices)) %>%
      arrange(subsample_indices)
    
    # Initialize keep flags
    time_diff$keep <- 0
    time_diff$keep2 <- 0
    
    # Mark rows to keep based on event sequence
    # duration is from the last merge to the first split
    # and from the first split to the next merge
    
    for (i in 1:(nrow(time_diff) - 1)) {
      if (time_diff$event[i] == "merge" && time_diff$event[i + 1] == "split") {
        time_diff$keep[i] <- 1
        time_diff$keep2[i + 1] <- 1
      }
    }
    
    # Filter rows
    time_diff$keep3 <- time_diff$keep + time_diff$keep2
    time_diff <- time_diff[time_diff$keep3 == 1, c("subsample_indices", "event")]
    
    # Split by event
    splits_time_diff <- time_diff[time_diff$event == "split", ]
    merge_time_diff <- time_diff[time_diff$event == "merge", ]
    
    time_diff <- data.frame(cbind(splits_time_diff, merge_time_diff))
    time_diff <- time_diff[, -c(2, 4)]
    colnames(time_diff) <- c("splits_t", "merge_t")
    
    # Calculate differences
    time_diff$diff <- as.numeric(time_diff$splits_t) - as.numeric(time_diff$merge_t)
    time_diff$diff_time_hour <- round((time_diff$diff * 10) / 60, 1)
    
    # Store result
    time_diff_list[[n]] <- time_diff
  }
  
  return(time_diff_list)
}



# Function to calculate dyadic distances
calculate_dyadic_distances <- function(data) {
  distances <- data %>%
    inner_join(data, by = "Time", suffix = c(".1", ".2")) %>%
    filter(Individual.1 < Individual.2) %>%
    mutate(Distance = sqrt((X_UTM.1 - X_UTM.2)^2 + (Y_UTM.1 - Y_UTM.2)^2)) %>%
    select(Time, Individual.1, Individual.2, Distance)
  
  return(distances)
}


# Function to generate pairwise data frame for a specific individual
generate_distance_df <- function(individual_id) {
  individual_data <- dyadic_distances %>%
    filter(Individual.1 == individual_id | Individual.2 == individual_id) %>%
    mutate(Other_Individual = if_else(Individual.1 == individual_id, Individual.2, Individual.1)) %>%
    select(Time, Other_Individual, Distance) %>%
    mutate(Focal_Individual = individual_id)
  
  return(individual_data)
}


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
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%m-%d") # Adjust the datetime scale
}





compute_dyadic_metrics <- function(data) {
  inds <- sort(unique(c(data$Focal_Individual, data$Other_Individual)))
  
  dyad_array <- array(NA_real_, 
                      dim = c(length(inds), length(inds), 3),
                      dimnames = list(inds, inds, c("mean", "median", "sd"))
  )
  
  for (i in inds) {
    for (j in inds) {
      pair_data <- data$Distance[data$Focal_Individual == i & data$Other_Individual == j]
      
      if (length(pair_data) > 0) {
        dyad_array[i, j, "mean"]   <- mean(pair_data, na.rm = TRUE)
        dyad_array[i, j, "median"] <- median(pair_data, na.rm = TRUE)
        dyad_array[i, j, "sd"]     <- sd(pair_data, na.rm = TRUE)
      }
    }
  }
  
  return(dyad_array)
}
