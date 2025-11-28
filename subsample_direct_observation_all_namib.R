#this script is using all the namibia groups to run the subsampling analysis to simulate traditional field methods
#compare accuracy of results when visiting at 6am, 10am, and 2pm

library(hms)
library(ggplot2)
library(dplyr)
library(fields)
library(asnipe)

code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
plot_dir <-  'C:/Users/egrout/Dropbox/lion/results/level1/'

setwd(code_dir)
source('coati_function_library_V1.R')
source('lion_functions.R')

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
load('allnamib_xy_2hour_level1.RData') #level 1 is filtering to gps at the correct 2h/4h interval, getting the closest fix within 20 mins of the desired hours
load('allnamib_ids_withdod.RData')

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

make_plot <- F
if(make_plot == T){
png(height = 1100, width = 1200, units = 'px', filename = paste0(plot_dir,'allnamib_network_60m_2yr_withdod.png'))

visualize_lion_network(ffnet_reorder, all_ids[new_order,])

dev.off()
}


#next is split data for each group and coalition
subgroup_data_perpride <- list()

for(i in unique(all_ids$pride_id)){
  
  if(is.na(i)) next #skip if pride id is NA
  
  #index xs and ys for each prides individuals
  xs_pride_i <- xs[which(all_ids$pride_id == i),, drop = FALSE]
  ys_pride_i <- ys[which(all_ids$pride_id == i),, drop = FALSE]
  
  subgroup_data_i <- get_subgroup_data(xs_pride_i, ys_pride_i, R)
  subgroup_data_perpride[[i]] <- subgroup_data_i
  
  
}


for(group in names(subgroup_data_perpride)){

  n_inds <- nrow(subgroup_data_perpride[[group]]$ind_subgroup_membership)
  
  if(length(n_inds) < 2) next
  
  ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

  #going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      
      #getting subgroup id for individual i and j
      sub_ids_i <- (subgroup_data_perpride[[group]]$ind_subgroup_membership)[i,]
      sub_ids_j <- (subgroup_data_perpride[[group]]$ind_subgroup_membership)[j,]
      
      #computing edge weight (fraction of time in same subgroup)
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1:n_inds)
ffnet_reorder <- ff_net[new_order, new_order]

make_plot <- T
if(make_plot == T){
  #png(height = 1100, width = 1200, units = 'px', filename = paste0(plot_dir,'allnamib_network_60m_2yr_withdod.png'))
  
  visualize_lion_network(ffnet_reorder, all_ids[new_order,])
  
  #dev.off()
}
}






