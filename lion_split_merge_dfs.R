#getting merge dataframe to compare duration of splits

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/lion/data/processed/"
code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
coati_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/lion/results/level0/'
gps_file <- "lion_xy_15min_level0.RData" #level0 is when Venus is not removed
id_file <- 'lion_ids.RData'

R <- 100

library(fields)
library(viridis)
library(hms)

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

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, R)



#-----SPLITS---------


#run through time by time and check if each time is a split
splits <- c()
for(t in 1:(n_times-1)){
  
  #get subgroup membership now and later
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(subgroups_now))==0 | sum(!is.na(subgroups_later))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
  n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(subgroups_now)==1)
  singletons_later <- sum(table(subgroups_later)==1)
  
  #determine if this time step is a split
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if((n_subgroups_now == 1 )
     & n_subgroups_later > n_subgroups_now
     & singletons_later==0
  ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & ((singletons_now + 1) == n_subgroups_now )
     & n_subgroups_later > n_subgroups_now
     & singletons_now == singletons_later
  ){
    splits <- c(splits, t)
  }
  
  
}

#make a data frame of splits, with original group and subgroups
splits_df <- data.frame(t = splits, orig_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)
for(i in 1:nrow(splits_df)){
  t <- splits_df$t[i]
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #transfer NAs from now to later and later to now so errenous splits are not detected
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(subgroups_now)
  orig_subgroup_members <- which(subgroups_now == orig_subgroup_number)
  
  #store original group membership in data frame
  splits_df$orig_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_later[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_later==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
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

#save dataframe as its needed for the figure1_fissionfusion_plot.R
save(splits_df, file = paste0(data_dir, "lion_splits_df.Rdata"))

#--------------------------------------------------------


#--------MERGE------------------


#run through time by time and check if each time is a merge
merge <- c()

for(t in 1:(n_times-1)){
  
  #get subgroup membership now and previous
  merge_group <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(merge_group))==0 | sum(!is.na(subgroups_previously))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  merge_group[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(merge_group))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_merge_group <- length(unique(merge_group[!is.na(merge_group)]))
  n_subgroups_previously <- length(unique(subgroups_previously[!is.na(subgroups_previously)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(merge_group)==1)
  singletons_later <- sum(table(subgroups_previously)==1)
  
  #determine if this time step is a merge
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if(n_merge_group==1 
     & n_subgroups_previously >1
     & singletons_later==0
  ){
    merge <- c(merge, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_merge_group > 1 
     & (singletons_now + 1) == n_merge_group
     & n_subgroups_previously > n_merge_group
     & singletons_now == singletons_later
  ){
    merge <- c(merge, t)
  }
  
  
}

#this seems to work


#make a data frame of merges, with merged group and subgroups
merge_df <- data.frame(t = merge, merge_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)


i=5
for(i in 1:nrow(merge_df)){
  t <- merge_df$t[i]
  merge_group <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #transfer NAs from now to later and later to now
  merge_group[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(merge_group))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(merge_group)
  orig_subgroup_members <- which(merge_group == orig_subgroup_number)
  
  #store original group membership in data frame
  merge_df$merge_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_previously[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_previously==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
    #really hacky shit to get R to put lists into a data frame :(
    if(j==1){
      merge_df$sub1[i] <- list(orig_inds_in_group)
    }
    if(j==2){
      merge_df$sub2[i] <- list(orig_inds_in_group)
    }
    if(j==3){
      merge_df$sub3[i] <- list(orig_inds_in_group)
    }
    if(j==4){
      merge_df$sub4[i] <- list(orig_inds_in_group)
    }
    if(j==5){
      merge_df$sub5[i] <- list(orig_inds_in_group)
    }
  }
}

#number in each subgroup
merge_df$n_merge <- sapply(merge_df$merge_group, function(x){return(sum(!is.na(x)))})
merge_df$n_sub1 <- sapply(merge_df$sub1, function(x){return(sum(!is.na(x)))})
merge_df$n_sub2 <- sapply(merge_df$sub2, function(x){return(sum(!is.na(x)))})
merge_df$n_sub3 <- sapply(merge_df$sub3, function(x){return(sum(!is.na(x)))})

save(merge_df, file = paste0(data_dir, "lion_merge_df.Rdata"))  





