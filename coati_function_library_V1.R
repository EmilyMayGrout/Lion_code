#library of general functions

#LIBRARIES
library(dbscan)
#library(rgdal)
library(sp)
library(lubridate)
library(stringr)

#GRAPHICS FOR WINDOWS AND MAC
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

#--------------------------------------------------------------------------

#This function computes information about subgroup membership over time
#Inputs:
# R [numeric]: radius for DBSCAN
# xs, ys [n_inds x n_times matrix]: x and y locations of individuals 
#Outputs:
# subgroup_data: object that contains:
#   $ind_subgroup_membership [n_inds x n_times matrix]: matrix giving the subgroup id for each individual at each time
#   $n_subgroups [n_times vector]: vector of the number of subgroups over time
#   $subgroup_counts [max_n_subgroups x n_times matrix]: matrix giving the number of individuals in each subgroup (with NAs when there are fewer subgroups than the max)
#   $R [numeric]: radius used in DBSCAN
get_subgroup_data <- function(xs, ys, R){
  
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  sub_groups <- matrix(NA, nrow = n_inds, ncol = n_times)
  
  for (t_idx in 1:n_times){
    
    #getting the x and y coordinates for each time point
    x <- xs[ , t_idx]
    y <- ys[, t_idx]
    
    non_nas <- which(!is.na(x))
    
    x_non_na <- x[non_nas]
    y_non_na <- y[non_nas]
    
    n_tracked <- length(non_nas)
    
    if(n_tracked >= 1){
      #scan for each time column to get the cluster number for each individual
      groups <- dbscan(x= cbind(x_non_na, y_non_na), eps = R, minPts = 1)
      #saving the cluster for each time into the matrix
      sub_groups[non_nas, t_idx] <- groups$cluster
    }
    
  }
  
  #get max number of groups per time
  ngroups <- suppressWarnings(apply(sub_groups, 2, max, na.rm = T)) #suppress warnings because if all NAs it gives a warning, this is fine because we deal with it in the next step by removing the infinite maxes
  #remove infinite values
  ngroups[is.infinite(ngroups)] <- NA

  #now want to look at size of the sub-groups as the histogram shows 2 subgroups being most frequent 
  #but this might be because of one individual not being with the group
  
  group_max <- max(ngroups, na.rm = T)
  
  subgroup_counts <- matrix(NA, nrow = group_max, ncol = n_times)
  
  #for loop for getting the number of individuals in each sub group every 10 mins (when the radius for defining within group is 50m, there are maximum 5 groups in full dataset)
  for(i in 1:n_times){
    
    #going through each time stamp
    current_group <- sub_groups[, i]
    #number of individuals in each subgroup 
    counts <- table(current_group)
    #to remove NA's when there is no data from all individuals
    if(length(counts) >= 1){
      #putting the data into the matrix
      subgroup_counts[1:length(counts),i] <- sort(counts)
      
    }
  }
  
  #store all data in a list
  subgroup_data <- list()
  subgroup_data$ind_subgroup_membership <- sub_groups
  subgroup_data$n_subgroups <- ngroups
  subgroup_data$subgroup_counts <- subgroup_counts
  subgroup_data$R <- R
  
  return(subgroup_data)
  
}

#-----------------------------------------------------------------

#jet.colors <- colorRampPalette(c("black" ,"red4","lightpink3", "lightpink1", "#FFCA00"))

#Make a network visualization with colors
#INPUT:
#net[adjacency matrix of the network]
#coati_ids[dataframe of the coati ids]
#OUTPUT:
#plot
visualize_network_matrix_galaxy <- function(net, coati_ids){
  
  zmin <- min(net, na.rm=T)
  zmax <- max(net, na.rm=T)
  par(mgp=c(3, 1, 0), mar=c(11,11,3,3)) #bottom, left, top, and right ##CHANGE THESE VALUES IF DOING FORLOOP OF THE MATRIX PLOTS AND THEY DON'T FIT mar=c(11,8,3,5)) and legend.mar = 12
  image.plot(net, col = viridis(256), zlim=c(zmin,zmax), xaxt= 'n', yaxt = 'n', legend.cex = 7, legend.width = 1.3,legend.mar = 6, axis.args=list(cex.axis=2))
  axis(1, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2, cex.axis=1.8)
  axis(2, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2,  cex.axis=1.8)
  
  #points(rep(-.07, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  #points(seq(0,1,length.out=nrow(net)),rep(-.07,n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  
  #for the all group matrices - genetics
  points(rep(-.07, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  points(seq(0,1,length.out=nrow(net)),rep(-.07,n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  
}
#for all groups genetics matrix, changed rep to -.018 from -.07

#changed the size of the labels with cex.axis = 1.5, default is 1
#also changed the size of the legend axis values with: axis.args=list(cex.axis=2), remove if want default
#legend width was 1.3, changed to 5 for the multiple matrix plots
#cex in points change points size
#-0.8 to move the coloured ID points out of the matrix plot
#legend.cex = 5 changed to 7

visualize_network_matrix_presedente <- function(net, coati_ids){
  
  zmin <- min(net, na.rm=T)
  zmax <- max(net, na.rm=T)
  par(mgp=c(3, 1, 0), mar=c(11,11,3,3)) #bottom, left, top, and right ##CHANGE THESE VALUES IF DOING FORLOOP OF THE MATRIX PLOTS AND THEY DON'T FIT mar=c(11,8,3,5)) and change legend mar from 6 to 9
  #par(mgp=c(3, 1, 0), mar=c(11,8,3,8))
  image.plot(net, col = viridis(256), zlim=c(zmin,zmax), xaxt= 'n', yaxt = 'n', legend.cex = 7, legend.width = 1.3,legend.mar = 6, axis.args=list(cex.axis=2))
  axis(1, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2, cex.axis=1.8)
  axis(2, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2,  cex.axis=1.8)
  
  points(rep(-.05, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  points(seq(0,1,length.out=nrow(net)),rep(-.05,n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
}

visualize_network_matrix_trago <- function(net, coati_ids){
  
  zmin <- min(net, na.rm=T)
  zmax <- max(net, na.rm=T)
  par(mgp=c(3, 1, 0), mar=c(9,9,3,1)) #bottom, left, top, and right
  #change legend.mar to 8 for plotting multiple together, not legend axis
  image.plot(net, col = viridis(256), zlim=c(zmin,zmax), xaxt= 'n', yaxt = 'n', legend.cex = 7, 
             legend.width = 1.3,legend.mar = 6, legend.line = 1, axis.args=list(cex.axis=2))
  axis(1, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2, cex.axis=1.8)
  axis(2, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2,  cex.axis=1.8)
  
  points(rep(-.1, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2.5)
  points(seq(0,1,length.out=nrow(net)),rep(-.1,n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2.5)
}



#---------------------------------------------------------------------------

#Function for getting distance between individuals over time
#INPUT:
#xs, ys, radius
#OUTPUT
#object that contains:
#   $distance over time: distance between each individual at each time point
#   $proximity_network: probability that i and j are within a distance of R
#   $r_within [numeric]: radius used for proximity network

get_proximity_data <- function(xs, ys, r_within){
  
  #get n_inds and n_times
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #get distance between individuals at these subset times
  
  #make array to store data into
  dist_over_time <- array(NA, dim = c(n_inds, n_inds, n_times))
  
  #for loop through each individual with each individual for each time point (where all individuals are together and have gps point)
  for(t in 1:ncol(xs)){
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        
        #get x and y coordinates to calculate distance using pythagorus
        xcoord_i <- xs[i,t]
        xcoord_j <- xs[j,t]
        dx <- (xcoord_i - xcoord_j)
        ycoord_i <- ys[i,t]
        ycoord_j <- ys[j,t]
        dy <- (ycoord_i - ycoord_j)
        dist <- sqrt((dx)^2 + (dy)^2)
        
        #put the statements into the array at the correct position in the correct dimension
        dist_over_time[i, j, t] <- dist
      }
    }
  }
  
  #threshold the distances to determine if individuals are in proximity at each time
  net_over_time <- dist_over_time < r_within
  
  #construct proximity network
  proximity_net <- apply(net_over_time, MARGIN = c(1,2), FUN = mean, na.rm=T)
  diag(proximity_net) <- NA
  
  #store output
  proximity_data <- list()
  proximity_data$proximity_net <- proximity_net
  proximity_data$dist_over_time <- dist_over_time
  proximity_data$r_within <- r_within
  
  #return output
  return(proximity_data)
  
}

#-----------------------------------------------------------------

#Randomise the splits and construct new randomized data frame
#Keep original group the same, keep subgroup sizes the same, but randomize who goes to which subgroup
#INPUTS:
# splits_df: data frame with information about splits from real data
#OUTPUTS:
# splits_df_random: randomised version of splits data frame
randomise_splits <- function(splits_df){
  
  #start with the splits_df data frame from the data
  splits_df_rand <- splits_df
  
  #loop over each row and randomize who goes where
  for (i in 1:nrow(splits_df_rand)){
    
    #get the original group and sizes of subgroups for one row
    row <- splits_df[i,]
    orig_group <- row$orig_group[[1]]
    n_sub1 <- row$n_sub1
    n_sub2 <- row$n_sub2
    n_sub3 <- row$n_sub3
    
    #shuffle the original group to a random order
    orig_group_shuffled <- sample(orig_group)
    
    #allocate them to groups
    sub1 <- orig_group_shuffled[1:n_sub1]
    sub2 <- orig_group_shuffled[(length(sub1)+1):(length(sub1)+n_sub2)]
    
    splits_df_rand$sub1[i] <- list(sub1)
    splits_df_rand$sub2[i] <- list(sub2)
    
    if(n_sub3 > 0){
      sub3 <- orig_group_shuffled[(length(sub1)+length(sub2)+1):(n_sub1+n_sub2+n_sub3)]
      splits_df_rand$sub3[i] <- list(sub3)
    }
    
  }
  
  return(splits_df_rand)
  
}


#-----------------------------------------------------------------

#Helper function for get_consistency
#Returns q if q <= 0.5 and 1-q if q > 0.5
dist_to_0_or_1 <- function(q){
  
  if(is.na(q)){
    return(NA)
  }
  if(q > 0.5) {return(1-q)
  }else{return(q)}
  
}

#-----------------------------------------------------------------

#Compute our metric of consistency, C
#Let q_dyad = p_dyad_together
#m_dyad = 1-q_dyad if q_dyad > 0.5 or = q_dyad if q_dyad < 0.5 (so it's the distance to 0 or 1, whichever is closer)
#C = 1 - 2*mean(m_dyad)
#At the end, C is a measure of consistency where C = 0 if all q_dyad are at 0.5 and C = 1 if all q_dyad are either 0 or 1
#Why do we multiply by 2 and subtract from 1? The multiplying by 2 puts the number between 0 and 1 instead of between 0 and 0.5. The 1 minus makes it be consistency instead of inconsistency
#INPUTS:
# p_dyad_together: [matrix] of probabilities of splitting together given you were both in the original group
#OUTPUTS:
# C: [numeric] metric of consistency
get_consistency <- function(p_dyad_together){
  
  m_dyad <- sapply(p_dyad_together, FUN = dist_to_0_or_1)
  C <- 1 - 2*mean(m_dyad, na.rm=T)
  return(C)
}

#-----------------------------------------------------------------

#function that takes in a splits dataframe and outputs a p_dyad_together matrix (probability of splitting together for each dyad)
#INPUTS:
# splits_df_local: [data frame] containing information on group splits
# n_inds_local: [numeric] number of individuals in the full group
#OUTPUTS:
# p_dyad_together_local: [matrix n_inds_local x n_inds_local]: probability of being together in a split, given you were both in the original group 
get_p_dyad_together <- function(splits_df_local, n_inds_local){
  
  #COMPUTE METRIC OF P(STAY TOGETHER | originally together) for each dyad
  p_dyad_together_local <- array(NA, dim = c(n_inds_local,n_inds_local))
  #loop over dyads
  for(i in 1:(n_inds_local-1)){
    for(j in (i+1):n_inds_local){
      
      #find rows where they were both in the original group
      originally_together_rows <- which(unlist(lapply(splits_df_local$orig_group, FUN = function(x){return(i %in% x & j %in% x)})))
      
      #how many times do they end up in the same group
      still_together <- 0
      for(r in originally_together_rows){
        if(i %in% splits_df_local$sub1[r][[1]] & j%in% splits_df_local$sub1[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub2[r][[1]] & j %in% splits_df_local$sub2[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub3[r][[1]] & j %in% splits_df_local$sub3[r][[1]]){
          still_together <- still_together + 1
        }
        
      }
      
      p_dyad_together_local[i,j] <- still_together / length(originally_together_rows)
      
    }
  }
  return(p_dyad_together_local)
  
}




