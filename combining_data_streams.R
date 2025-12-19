#joining the high res Omukutu pride data to the movebank data from all namib after both were preprocessed to create the xs and ys matrices
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")

load('allnamib_xy_2hour_level2.RData') #level 1 is filtering to gps at the correct 2h/4h interval, getting the closest fix within 20 mins of the desired hours, level 2 is finding closest GPS within 20 mins of target hour at night and the closest point within +130 mins and -20 mins of day target hours (when they shouldn't move much)
load('allnamib_ids_withdod_V2.RData')

xs_all <- xs
ys_all <- ys
ts_all <- ts

#load Omukutu pride data
load("lion_xy_15min_level0_2yr.RData")
load("lion_ids.RData")

xs_omukutu <- xs
ys_omukutu <- ys
omukutu_ids <- lion_ids
ts_omukutu <- ts

rm("xs", "ys", "ts","lion_ids")

#next is replacing the omukutu inds in the xs_all and ys_all with the data from the xs_omukutu and ys_omukutu

#indices in ts_omukutu that match ts_all
ts_idx <- which(ts_omukutu %in% ts_all)

#find corresponding column indices in ts_all
ts_all_idx <- match(ts_omukutu[ts_idx], ts_all)

#get row for Omukutu IDs in all_ids 
omukutu_rows <- which(all_ids$name %in% omukutu_ids$tag_id)

#match Omukutu IDs to all_ids
match_idx <- match(all_ids$name, omukutu_ids$tag_id)

#keep only rows that actually match for replacing
valid <- which(!is.na(match_idx))

#replace omukutu to its high res data using matched rows
xs_all[valid, ts_all_idx] <- xs_omukutu[match_idx[valid], ts_idx]
ys_all[valid, ts_all_idx] <- ys_omukutu[match_idx[valid], ts_idx]

#add NAs for 8, 12 and 16 hours in omukutu
ts_hours <- hour(ts_all)
cols_to_NA <- which(ts_hours %in% c(8, 12, 16))
xs_all[omukutu_rows, cols_to_NA] <- NA
ys_all[omukutu_rows, cols_to_NA] <- NA


#some collars record on uneven numbers so not sure what to do about these for downstream analyses, for now will only use the data collected at 0,2,4,6,10,14,18,22

#saving the combined dataset for full group analyses
xs <- xs_all
ys <- ys_all
ts <- ts_all

setwd("C:/Users/egrout/Dropbox/lion/data/processed/")
save(list=c('xs','ys','ts'), file = 'allnamib_omukutu_xy_2hour_level3.RData')



#looking at the target hours different individuals have data for per month
n_inds_all <- nrow(xs_all)
hours <- 0:23

# initialize matrix
hour_counts_all <- matrix(0, nrow = n_inds_all, ncol = length(hours))
rownames(hour_counts_all) <- all_ids$name
colnames(hour_counts_all) <- hours

# fill counts per hour
for(i in 1:n_inds_all){
  has_data <- !is.na(xs_all[i, ]) & !is.na(ys_all[i, ])
  hrs <- hour(ts_all[has_data])
  hour_counts_all[i, ] <- table(factor(hrs, levels = hours))
}

# convert to long format for ggplot
hour_counts_long_all <- as.data.frame(hour_counts_all) %>%
  tibble::rownames_to_column(var = "individual") %>%
  pivot_longer(cols = -individual, names_to = "hour", values_to = "count")

ggplot(hour_counts_long_all, aes(x = hour, y = individual, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "grey90", high = "steelblue") +
  theme_minimal() +
  labs(x = "Hour of day", y = "Individual", fill = "Data points") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Identify the row for XPL-138
row_idx <- which(all_ids$name == "XPL-138")

# Extract timestamps and presence for this individual
ts_ind <- ts_all
data_present <- !is.na(xs_all[row_idx, ]) | !is.na(ys_all[row_idx, ])

# Keep only times with data
ts_present <- ts_ind[data_present]

# Extract year-month and hour
month_year <- format(ts_present, "%Y-%m")
hours <- hour(ts_present)

# Create a table of hours per month
hours_per_month <- tapply(hours, month_year, function(x) sort(unique(x)))
hours_per_month

#the dates where the usual 2h/4h sampling for XPL-138 and XPL-137 are not matching, so may be unable to use this group for long-term trends 




