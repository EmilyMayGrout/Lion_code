
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(sf)
library(tidyverse)
library(fields)
library(viridis)
library(tidygraph)
library(reshape2)
library(RColorBrewer)

code_dir <- 'C:/Users/egrout/Dropbox/lion/Lion_code/'
plot_dir <-  'C:/Users/egrout/Dropbox/lion/results/level0/'

setwd(code_dir)
source('coati_function_library_V1.R')
source('lion_functions.R')

#load data downloaded from movebank 
all_nab <- read.csv("../data/raw/all_namibia.csv")

#Directories
indir <-"../data/raw/gps/all_nab/"
outdir <- "../data/processed/"
metadatadir <- "../data/raw/metadata/"

setwd(indir)

#looking at which inds have data for each month for each year

#put time in posixct format
all_nab$timestamp <- as.POSIXct(all_nab$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")
all_nab$month_year <- format(all_nab$timestamp, "%Y-%m")

presence_matrix <- all_nab %>%
  group_by(individual.local.identifier, month_year) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = month_year, values_from = present, values_fill = 0)

# Convert to long format for ggplot
presence_long <- presence_matrix %>%
  pivot_longer(-individual.local.identifier, names_to = "month_year", values_to = "present")

# Plot heatmap
gg <- ggplot(presence_long, aes(x = month_year, y = individual.local.identifier, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  theme_classic() +
  labs(x = "Month-Year", y = "Individual", fill = "Data Present") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg

#ggsave(filename = "C:/Users/egrout/Dropbox/lion/results/level0/data_presence_all_namib.png", plot = gg, width = 15, height = 8, dpi = 300)

count_matrix <- all_nab %>%
  group_by(individual.local.identifier, month_year) %>%
  summarise(count = n(), .groups = "drop")

count_matrix$month_year <- factor(count_matrix$month_year, levels = sort(unique(count_matrix$month_year)))

pp <- ggplot(count_matrix, aes(x = month_year, y = individual.local.identifier, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey90") +
  theme_bw() +
  labs(
    title = "Data Point Counts per Individual per Month-Year",
    x = "Month-Year",
    y = "Individual",
    fill = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pp
#ggsave(filename = "C:/Users/egrout/Dropbox/lion/results/level0/data_presence_all_namib3.png", plot = pp, width = 15, height = 10, dpi = 300)


#-----------------------------------------------

#look at the proportion of time within 500m at 10am (UTC time), midday local time

load('C:/Users/egrout/Dropbox/lion/data/processed/allnamib_xy_2hour_level0.RData')
load('C:/Users/egrout/Dropbox/lion/data/processed/allnamib_latlon_2hour_level0.RData') 


lion_ids <- data.frame(gsub('.txt','',all_files))
colnames(lion_ids) <- "name"

all_ids <- read.csv("C:/Users/egrout/Dropbox/lion/data/raw/metadata/namib_lion_ids.csv", header = F)

#missing metadata for "NPL-28" "NPL-27" "NPL-40"

colnames(all_ids) <- c("name", "sex", "group_type", "pride_id", "yob")

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#subset to times when the data is every 15 mins
non_na_counts <- colSums(!is.na(xs))

run_script <- F

R = 500

if(run_script == T){
subgroup_data <- get_subgroup_data(xs, ys, R)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

split_by_year <- F
year <- unique(format(ts, "%Y"))
ff_net_year <- array(NA, dim = c(n_inds, n_inds, length(year)))
dimnames(ff_net_year)[[3]] <- year

for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    if(split_by_year){
      for(y in 1:length(seq_along(year))){
        
        y <- year[y]
        
        # Logical index for matching year-month
        ts_filter <- format(ts, "%Y") == y
        
        sub_ids_i <- subgroup_data$ind_subgroup_membership[i, ts_filter]
        sub_ids_j <- subgroup_data$ind_subgroup_membership[j, ts_filter]
        
        ff_net_year[i, j, y] <- mean(sub_ids_i == sub_ids_j, na.rm = TRUE)
        
      }
    }
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA

rownames(ff_net) <- lion_ids$name
colnames(ff_net) <- lion_ids$name

#save(ff_net, ff_net_year, file = "C:/Users/egrout/Dropbox/lion/data/processed/ff_net_500m.RData")
}

load("C:/Users/egrout/Dropbox/lion/data/processed/ff_net_500m.RData")


#reordering the lion_ids by pride ID
lion_ids1 <- lion_ids[order(lion_ids$pride_id), ]
#changing ff_net order to the pride ID order in lion_ids1
ff_net1 <- ff_net[lion_ids1$name, lion_ids1$name]

#plot just the prides, just the coalitions, or all inds
#plot partial matrix of females up one side and males along the other -- so focus on pride-coalition associations

groups <- unique(na.omit(lion_ids1$pride_id))
colors <- setNames(rainbow(length(groups)), groups)
individual_colors <- setNames(colors[as.character(lion_ids1$pride_id)], lion_ids1$name)

generate_longData <- function(ff_net1, lion_ids1, group_type = "all", individual_colors) {
  lion_ids2 <- lion_ids1[!is.na(lion_ids1$group_type), ]
  
  if (group_type == "pride") {
    inds <- lion_ids2$name[lion_ids2$group_type == "pride"]
  } else if (group_type == "coalition") {
    inds <- lion_ids2$name[lion_ids2$group_type == "coalition"]
  } else {
    inds <- lion_ids1$name
  }
  
  # Subset color mapping
  colours_sub <- individual_colors[inds]
  
  # Subset matrix and mask lower triangle
  sub_net <- ff_net1[inds, inds]
  sub_net[lower.tri(sub_net)] <- NA
  
  longData <- reshape2::melt(sub_net, varnames = c("Var1", "Var2"), value.name = "value")
  longData$Var1 <- factor(longData$Var1, levels = inds)
  longData$Var2 <- factor(longData$Var2, levels = inds)
  
  return(list(data = longData, inds = inds, colours = colours_sub))
}


#generate subsets
group_type <- "all"
matrix_data <- generate_longData(ff_net1, lion_ids1, group_type = group_type, individual_colors)

n <- length(unique(matrix_data$data$Var1))
a <- ggplot(matrix_data$data, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "slateblue4", high = "yellow", na.value = "white") +
  
  # Add vertical and horizontal grid lines between cells
  geom_vline(xintercept = seq(0.5, n + 0.5, by = 1), color = "snow2", linewidth = 0.2) +
  geom_hline(yintercept = seq(0.5, n + 0.5, by = 1), color = "snow2", linewidth = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.4, color = matrix_data$colours),
    axis.text.y = element_text(size = 9, color = matrix_data$colours),
    panel.grid = element_blank(),  # turn off default grid
    panel.border = element_blank() # keep it clean
  ) +
  labs(x = "", y = "", fill = "Proportion", title = paste0("Proportion of time ", group_type, " individuals are within ", R, "m of one another"))

a

if(group_type == "all"){
  # Create pride-level color mapping (one color per pride)
  pride_levels <- unique(na.omit(lion_ids1$pride_id))
  pride_colors <- setNames(rainbow(length(pride_levels)), pride_levels)
  
  # Dummy data for the legend (not actually plotted)
  pride_legend_df <- data.frame(pride_id = pride_levels)
  
  gg <- a +
    # Add dummy points to create legend entries
    geom_point(data = pride_legend_df,
               aes(x = NaN, y = NaN, color = pride_id),
               inherit.aes = FALSE, show.legend = TRUE) +
    # Manual color scale matching pride_id to colors
    scale_color_manual(name = "Pride", values = pride_colors) +
    guides(fill = guide_colorbar(title.position = "top"),
           color = guide_legend(title.position = "top")) +
    theme(
      legend.position = "right",
      legend.box = "horizontal",
      legend.title = element_text(hjust = 0.2)
    )
  
  gg
}


#rerunning plots for prides only, coalitions only, and all to put next to each other
#all <- p+c+gg
#all

#ggsave(filename = paste0("C:/Users/egrout/Dropbox/lion/results/level0/proptogether_", R, "m_namib_correct.png"), plot = all, width = 20, height = 6, dpi = 300)



#plot females on one axes and males on the other

# Subset names based on sex
lion_ids2 <- lion_ids1[!is.na(lion_ids1$sex), ]
males <- lion_ids2$name[lion_ids2$sex == "M"]
females <- lion_ids2$name[lion_ids2$sex == "F"]

# Subset the matrix for males (rows) vs females (columns)
mf_matrix <- ff_net1[males, females]

# Melt the matrix
longData_mf <- reshape2::melt(mf_matrix, varnames = c("Var1", "Var2"))

# Set the factor levels to preserve order
longData_mf$Var1 <- factor(longData_mf$Var1, levels = males)
longData_mf$Var2 <- factor(longData_mf$Var2, levels = females)

# Plot dimensions
n_row <- length(males)
n_col <- length(females)

#use individual_colors for axis text color
p_mf <- ggplot(longData_mf, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "slateblue4", high = "yellow", na.value = "white") +
  
  # Add custom grid lines
  geom_vline(xintercept = seq(0.5, n_col + 0.5, by = 1), color = "snow2", linewidth = 0.2) +
  geom_hline(yintercept = seq(0.5, n_row + 0.5, by = 1), color = "snow2", linewidth = 0.2) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 9, angle = 90, hjust = 1, vjust = 0.4,
      color = individual_colors[as.character(females)]  # female colors
    ),
    axis.text.y = element_text(
      size = 9,
      color = individual_colors[as.character(males)]  # male colors
    ),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  labs(
    x = "Females", y = "Males",
    fill = "Proportion",
    title = paste0("Proportion of time Males vs Females are within ", R, "m")
  )

p_mf

#ggsave(filename = paste0("C:/Users/egrout/Dropbox/lion/results/level0/males_vs_females_", R, "m_namib.png"), plot = p_mf, width = 7, height = 6, dpi = 300)


#plotting the raw movement data

max_xs <- max(xs, na.rm = T)
max_ys <- max(ys, na.rm = T)
min_xs <- min(xs, na.rm = T)
min_ys <- min(ys, na.rm = T)

#using original ids order
groups <- unique(na.omit(lion_ids$pride_id))
colors <- setNames(rainbow(length(groups)), groups)
individual_colors <- setNames(colors[as.character(lion_ids$pride_id)], lion_ids$name)

colors <- individual_colors[lion_ids$name]

# Set up an empty plot with the right limits
plot(NA, xlim = c(min_xs, max_xs), ylim = c(min_ys, max_ys), xlab = "X", ylab = "Y", type = "n")

# Add lines for each individual (each row in xs and ys)
for (i in 1:nrow(xs)) {
  lines(xs[i, ], ys[i, ], col = colors[i])  # You can replace `col = i` with a custom color vector if needed
}

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------


#clean the environment to look at each prides subgrouping dynamics in detail

rm(list=c("a", "combined_df", "count_matrix", "gg", "latlon_i", "lats", "lion_ids1", "lion_ids2", "longData_mf", "lons", "matrix_data", "p_mf", "pp", "presence_long", "presence_matrix", "pride_legend_df", "tagdata", "utm_i"))


#reordering the lion_ids by pride ID
lion_ids1 <- lion_ids[order(lion_ids$pride_id), ]
#changing ff_net order to the pride ID order in lion_ids1
ff_net1 <- ff_net[lion_ids1$name, lion_ids1$name]

ff_net1[ff_net1==0] <- NaN

#png(height = 1500, width = 1550, units = 'px', filename = paste0(plot_dir,'allnamib_500m_network.png'))
visualize_lion_network(ff_net1, lion_ids1)
#dev.off()













