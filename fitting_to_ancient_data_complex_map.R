# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)
library(rnaturalearth)
library(sp)
library(ggplot2)
library(data.table)
library(sf)
library(geosphere)
library(scatterpie)
library(ggmap)
library(Rmpfr)
library(robustbase)
library(MASS)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in simulation files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
wd = ("/path/to/simulation/data")

setwd(wd)
map_size_km = 3700
map = ne_countries(scale = "medium", returnclass = "sf")

# Map extents
y_min = 0       # Minimum y-coordinate
y_max = map_size_km    # Maximum y-coordinate
lat_min = 35   # Minimum latitude
lat_max = 72    # Maximum latitude

x_min = 0       # Minimum y-coordinate
x_max = map_size_km    # Maximum y-coordinate
lon_min = -10   # Minimum latitude
lon_max = 30    # Maximum latitude

# -------------------------------
# If landscape is complex map
# -------------------------------

# Specify and select files
non_square_file_names = list.files(".", pattern="sim_pop*", full.names = TRUE)

if (length(non_square_file_names) != 0)
{
  # Sort Files
  non_square_file_names = mixedsort(non_square_file_names)
  # Read in File
  non_square_files = lapply(non_square_file_names, read.csv)
}

# -----------------------------------------
# If there's an ancestry sample file
# -----------------------------------------

# Specify and select files
ancestry_sample_names = list.files(".", pattern="sim_ancestry_sample*", full.names = TRUE)

if (length(ancestry_sample_names) != 0)
{
  # Sort Files
  ancestry_sample_names = mixedsort(ancestry_sample_names)
  # Read in File
  ancestry_sample_files = lapply(ancestry_sample_names, function(x) {
    head(read.csv(x, header = T), 100000)
  })
}


# -----------------------------------------
# Read in parameter file
# -----------------------------------------

# Read in Param File
param_file_name = "param_inputs_clean.txt"
param_file = lapply(param_file_name, read.csv)

# ----------------------------------------------------------------------------------------------------------------
# Adding parameter combinations to the ancestry sample table
# ----------------------------------------------------------------------------------------------------------------

if (length(ancestry_sample_names) != 0)
{
  ancestry_sample_data = lapply(1:length(ancestry_sample_files), function(x) 
  {
    params = strsplit(as.character(param_file[[1]][x,]), split = "   ")[[1]]
    params = matrix(params, nrow = 1)
    colnames(params) = gsub(sapply(strsplit(params, split = " = "), function(x) x[1]), pattern = " ", replacement = "_")
    data = ancestry_sample_files[[x]]
    
    cbind(params, data)
  }
  )
  
  # Create data table with simulation output data
  ancestry_sample_data = rbindlist(ancestry_sample_data)
  
  # Add columns that remove extra strings in parameter set
  ancestry_sample_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_sample_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_sample_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_sample_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Movement_X_n := as.numeric(gsub(Movement_X, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Movement_Y_n := as.numeric(gsub(Movement_Y, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_sample_data[, Map_Style_n := as.numeric(gsub(Map_Style, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Water_Crossings_n := as.numeric(gsub(Water_Crossings, pattern = ".*\\s(.+)", replacement = "\\1"))]
  
  ancestry_sample_data$Farming_Ancestry = ancestry_sample_data$Farming_Ancestry * (494/496)
  
  # Set > 1 Values to 1
  ancestry_sample_data$Farming_Ancestry[ancestry_sample_data$Farming_Ancestry > 1] = 1.0
  
  ancestry_sample_data_farmers_only = ancestry_sample_data[ancestry_sample_data$Individual_Z != 0]
  
  ancestry_sample_data_farmers_only = data.table(ancestry_sample_data_farmers_only)
}

# ----------------------------------------------------------------------------------------------------------------
# Adding distances from lower corner and calculating ancestry by distance
# ----------------------------------------------------------------------------------------------------------------

#  (x1, y1) is the bottom right hand corner of the map (i.e., farmer origin)
x1 = map_size_km
y1 = 0

# Number of bins
num_parts = 100

# Calculate the distance of each individual from the origin
# sqrt((x2 – x1)^2 + (y2 – y1)^2)
ancestry_sample_data_farmers_only = ancestry_sample_data_farmers_only %>% mutate(Distance_from_Origin = sqrt((Individual_X-x1)^2 + (Individual_Y-y1)^2))

#Max Distance From Origin Sampled 
max_dist = max(ancestry_sample_data_farmers_only$Distance_from_Origin)

ancestry_sample_data_farmers_only$approximated_lat = lat_min + (ancestry_sample_data_farmers_only$Individual_Y - y_min) * (lat_max - lat_min) / (y_max - y_min)
ancestry_sample_data_farmers_only$approximated_long = lon_min + (ancestry_sample_data_farmers_only$Individual_X - x_min) * (lon_max - lon_min) / (x_max - x_min)

ancestry_sample_data_farmers_only_southern_inds = ancestry_sample_data_farmers_only[ancestry_sample_data_farmers_only$approximated_lat <= 45]
ancestry_sample_data_farmers_only_northern_inds = ancestry_sample_data_farmers_only[ancestry_sample_data_farmers_only$approximated_lat > 45]

#possible_cut_bins = c(0,262,524,786,1048,1310,1572,1834,2096,2358,2620,2882,3144,3406,3668,3930,4192,4454,4716,4978,5240)

possible_cut_bins = c(0, 52.4, 104.8, 157.2, 209.6, 262.0, 314.4, 366.8, 419.2, 471.6, 524.0, 576.4, 628.8, 681.2, 733.6, 786.0, 838.4, 890.8, 943.2, 995.6, 1048.0, 1100.4, 1152.8, 1205.2, 1257.6, 1310.0, 1362.4, 1414.8, 1467.2, 1519.6, 1572.0, 1624.4, 
                      1676.8, 1729.2, 1781.6, 1834.0, 1886.4, 1938.8, 1991.2, 2043.6, 2096.0, 2148.4, 2200.8, 2253.2, 2305.6, 2358.0, 2410.4, 2462.8, 2515.2, 2567.6, 2620.0, 2672.4, 2724.8, 2777.2, 2829.6, 2882.0, 2934.4, 2986.8, 3039.2, 3091.6, 3144.0, 
                      3196.4, 3248.8, 3301.2, 3353.6, 3406.0, 3458.4, 3510.8, 3563.2, 3615.6, 3668.0, 3720.4, 3772.8, 3825.2, 3877.6, 3930.0, 3982.4, 4034.8, 4087.2, 4139.6, 4192.0, 4244.4, 4296.8, 4349.2, 4401.6, 4454.0, 4506.4, 4558.8, 4611.2, 4663.6, 
                      4716.0, 4768.4, 4820.8, 4873.2, 4925.6, 4978.0, 5030.4, 5082.8, 5135.2, 5187.6, 5240.0)


# Bin based on distance from origin
ancestry_sample_data_farmers_only = ancestry_sample_data_farmers_only %>% 
  mutate(Distance_Bin = cut(Distance_from_Origin, breaks=possible_cut_bins, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                                                                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                                                                     37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                                                                                     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                                                                                     71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
                                                                                     88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)))
#labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))
ancestry_sample_data_farmers_only_southern_inds = ancestry_sample_data_farmers_only_southern_inds %>% 
  mutate(Distance_Bin = cut(Distance_from_Origin, breaks=possible_cut_bins, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                                                                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                                                                     37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                                                                                     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                                                                                     71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
                                                                                     88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)))
#labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))
ancestry_sample_data_farmers_only_northern_inds = ancestry_sample_data_farmers_only_northern_inds %>% 
  mutate(Distance_Bin = cut(Distance_from_Origin, breaks=possible_cut_bins, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                                                                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                                                                     37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                                                                                     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                                                                                     71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
                                                                                     88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)))
#labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))

# Calculate where the first bin midpoint point is
bin_size = (5240 / num_parts)

# Create column with mid point of each bin
ancestry_sample_data_farmers_only = ancestry_sample_data_farmers_only %>% mutate(Mid_Point_km = (as.numeric(Distance_Bin)*bin_size - bin_size/2))
ancestry_sample_data_farmers_only_southern_inds = ancestry_sample_data_farmers_only_southern_inds %>% mutate(Mid_Point_km = (as.numeric(Distance_Bin)*bin_size - bin_size/2))
ancestry_sample_data_farmers_only_northern_inds = ancestry_sample_data_farmers_only_northern_inds %>% mutate(Mid_Point_km = (as.numeric(Distance_Bin)*bin_size - bin_size/2))

# Group by param combination and calculate mean ancestry in each bin for each param combo
ancestry_sample_data_farmers_only = ancestry_sample_data_farmers_only %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Mean_Ancestry = mean(Farming_Ancestry))
ancestry_sample_data_farmers_only_southern_inds = ancestry_sample_data_farmers_only_southern_inds %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Mean_Ancestry = mean(Farming_Ancestry))
ancestry_sample_data_farmers_only_northern_inds = ancestry_sample_data_farmers_only_northern_inds %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Mean_Ancestry = mean(Farming_Ancestry))

# Get variance of sampling in each bin
ancestry_sample_data_farmers_only = ancestry_sample_data_farmers_only %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Ancestry_Variance = var(Farming_Ancestry))
ancestry_sample_data_farmers_only_southern_inds = ancestry_sample_data_farmers_only_southern_inds %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Ancestry_Variance = var(Farming_Ancestry))
ancestry_sample_data_farmers_only_northern_inds = ancestry_sample_data_farmers_only_northern_inds %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Ancestry_Variance = var(Farming_Ancestry))

setwd("/path/to/ancient/dna/data")

all_ancient_inds = data.table(read.table("LaPolice__SI_Table__qpAdm_estimates.txt", header = T, sep="\t"))

#' Keep only plausible model individuals with s.e < 0.022
filtered_inds = all_ancient_inds[se.p3 < 0.022 & p.value >= 0.01 &
                             (weight.p1 >= 0 & weight.p1 <=1) &
                             (weight.p2 >= 0 & weight.p2 <=1) &
                             (weight.p3 >= 0 & weight.p3 <=1)]

# Remove individuals with > 5% Steppe ancestry
filtered_inds <- filtered_inds[weight.p1 <= 0.05]

headers = data.table(read.table("qpadm_headers_complex.csv", header = T, sep=","))

# Steppe Ancestry Histogram (no ancestry > 5%)
#hist(filtered_inds$weight.p1)
# WHG Ancestry Historgram
#hist(filtered_inds$weight.p2)
# EEF Ancestry Historgram
#hist(filtered_inds$weight.p3)


# Get the distance of each sample frm the origin (Ankara)
library(geosphere)
# Convert 'lat' and 'lon' to numeric
filtered_inds[, `:=`(
  lat = as.numeric(as.character(lat)),
  lon = as.numeric(as.character(lon))
)]

# Define Ankara's coordinates
Ankara_coords <- c(lon = 35, lat = 40)
# Calculate the distance to Ankara in kilometers
# Assuming 'd' is your data.table and has columns 'lat' and 'lon'
filtered_inds[, distFromAnkara_km := distHaversine(cbind(lon, lat), Ankara_coords) / 1000]

filtered_inds$time_period[filtered_inds$Age <= 6500 & filtered_inds$Age >= 5500] = "Middle"
filtered_inds$time_period[filtered_inds$Age < 5500] = "Late"
filtered_inds$time_period[filtered_inds$Age > 6500] = "Early"

filtered_inds$location_bin[filtered_inds$lat <= 45] = "<= 45 Degrees"
filtered_inds$location_bin[filtered_inds$lat > 45] = "> 45 Degrees"

#filtered_inds = filtered_inds %>% 
  #mutate(Distance_Bin = cut(distFromAnkara_km, breaks=possible_cut_bins, labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))

filtered_inds = filtered_inds %>% 
  mutate(Distance_Bin = cut(distFromAnkara_km, breaks=possible_cut_bins, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                                                                                  20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                                                                  37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                                                                                  54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                                                                                  71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
                                                                                  88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)))

# This will give you a table of farming ancestry estimates and the 
# distance of each sample from the farming origin
filtered_inds[,.(distFromAnkara_km, weight.p3)] 

filtered_inds_southern_inds = filtered_inds[filtered_inds$location_bin == "<= 45 Degrees"]
filtered_inds_northern_inds = filtered_inds[filtered_inds$location_bin == "> 45 Degrees"]

ancestry_sample_data_farmers_only = data.table(ancestry_sample_data_farmers_only)
ancestry_sample_data_farmers_only_southern_inds = data.table(ancestry_sample_data_farmers_only_southern_inds)
ancestry_sample_data_farmers_only_northern_inds = data.table(ancestry_sample_data_farmers_only_northern_inds)

filtered_inds$approximated_y = y_min + (filtered_inds$lat - lat_min) * (y_max - y_min) / (lat_max - lat_min)
filtered_inds$approximated_x = x_min + (filtered_inds$long - lon_min) * (x_max - x_min) / (lon_max - lon_min)

filtered_inds_southern_inds$approximated_y = y_min + (filtered_inds_southern_inds$lat - lat_min) * (y_max - y_min) / (lat_max - lat_min)
filtered_inds_southern_inds$approximated_x = x_min + (filtered_inds_southern_inds$long - lon_min) * (x_max - x_min) / (lon_max - lon_min)

filtered_inds_northern_inds$approximated_y = y_min + (filtered_inds_northern_inds$lat - lat_min) * (y_max - y_min) / (lat_max - lat_min)
filtered_inds_northern_inds$approximated_x = x_min + (filtered_inds_northern_inds$long - lon_min) * (x_max - x_min) / (lon_max - lon_min)

# Read in Param File
#param_file = data.table(read.table("param_inputs_clean_filtered_triple_maps.csv", header = T, sep=","))
param_file = data.table(read.table("param_inputs_clean_filtered_complex_map_with_single_best_fit_run.csv", header = T, sep=","))

# ----------------------------------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------------------------------

interpolate = function(ancient_x_dist, sim_x_dist, sim_farming_ancestry){
  
  interpolation = approx(sim_x_dist, sim_farming_ancestry, xout=ancient_x_dist, method="linear", rule = 2)
  
  # Pull out the distance from the approximation list
  approximated_ancestry = interpolation[[2]]
  
  View(interpolation)
  
  # Return value
  return(approximated_ancestry)
  
}

main_ancestry_dt = headers

for (param in param_file$Learning_Prob_n)
  #for (param in param_file$Assortative_Mating_n)
{
  individual_param_run = ancestry_sample_data_farmers_only[Learning_Prob_n == param]
  
  sim_farming_ancestry = (individual_param_run$Bin_Mean_Ancestry)
  sim_x_dist = individual_param_run$Mid_Point_km
  ancient_x_dist = filtered_inds$distFromAnkara_km
  
  interpolated = interpolate(ancient_x_dist, sim_x_dist, sim_farming_ancestry)
  
  filtered_inds[, approximated_ancestry := as.numeric(0)]
  filtered_inds[, approximated_ancestry := interpolated]
  
  filtered_inds[, distance_from_sim := as.numeric(0)]
  filtered_inds[, distance_from_sim := (weight.p3 - approximated_ancestry)]
  
  filtered_inds[, distance_from_sim_squared := as.numeric(0)]
  filtered_inds[, distance_from_sim_squared := (distance_from_sim^2)]
  
  filtered_inds[, mean_distance_from_sim_squared := as.numeric(0)]
  filtered_inds[, mean_distance_from_sim_squared := mean(distance_from_sim_squared)]
  
  filtered_inds[, log_likelihood := as.numeric(0)]
  filtered_inds[, log_likelihood := (dnorm(weight.p3, approximated_ancestry, se.p3, log = T))]
  
  filtered_inds[, sum_log_likelihood := as.numeric(0)]
  filtered_inds[, sum_log_likelihood := sum(log_likelihood)]
  
  filtered_inds[, exponential_sum_log_likelihood := as.numeric(0)]
  filtered_inds[, exponential_sum_log_likelihood := (exp(sum_log_likelihood))]
  
  filtered_inds[, Learning_Prob_n := as.numeric(0)]
  filtered_inds[, Learning_Prob_n := param]
  
  #filtered_inds[, Assortative_Mating_n := as.numeric(0)]
  #filtered_inds[, Assortative_Mating_n := param]
  
  main_ancestry_dt = rbind(main_ancestry_dt, filtered_inds)
  
}

ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))
#ggplot(main_ancestry_dt) + geom_point(aes(Assortative_Mating_n, sum_log_likelihood))

ggplot(main_ancestry_dt) + geom_point(aes(distFromAnkara_km, weight.p3))

main_ancestry_dt_best_fit_quad = main_ancestry_dt[Learning_Prob_n <= 0.0025]
#main_ancestry_dt_best_fit_quad = main_ancestry_dt[Assortative_Mating_n >= 0.9]

quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)
#quadraticModel = lm(sum_log_likelihood ~ Assortative_Mating_n + I(Assortative_Mating_n^2), data=main_ancestry_dt_best_fit_quad)

l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])
#l_doubleprime = -(quadraticModel[["coefficients"]][["I(Assortative_Mating_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

LearningValues = seq(0, 0.006, 0.0001)
#AMValues = seq(0.9, 1, 0.001)

sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))
#sum_log_likelihoodPredict = predict(quadraticModel,list(Assortative_Mating_n=AMValues, Assortative_Mating_n2=AMValues^2))

# This just gets an approximate best fitting value- which of your tested params fits the best:
approximate_prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]

# To get the actual best fitting prediction put this in the console:
#           quadraticModel["coefficients"]

# Or do this to pull them out as variables

coefficients = quadraticModel$coefficients
coefficients = as.vector(coefficients)
Learning_Prob_n_coefficient = coefficients[2]
ILearning_Prob_n2_coefficient = coefficients[3]

# Then the prediction will equal this equation. Get the Learning_Probs from the coefficients
# prediction = ((-Learning_Prob_n) / (2*[I(Learning_Prob_n^2)]))

# Or do it out with code, not by hand

prediction = ((-as.numeric(Learning_Prob_n_coefficient)) / (2*as.numeric(ILearning_Prob_n2_coefficient)))

confidence_interval_floor = prediction-confidence_interval
confidence_interval_ceiling = prediction+confidence_interval

prediction_rounded <- round(prediction, 6)
confidence_interval_floor_rounded <- round(confidence_interval_floor, 6)
confidence_interval_ceiling_rounded <- round(confidence_interval_ceiling, 6)

formatted_prediction_with_conf = paste(prediction_rounded, sprintf("[%s, %s]", confidence_interval_floor_rounded, confidence_interval_ceiling_rounded))

max_like = ((quadraticModel$coefficients[1]) - ((quadraticModel$coefficients[2])^2/(4*(quadraticModel$coefficients[3]))))

setwd("/path/to/plotting/output/directory")

#if (!file.exists("refit")) {
  #dir.create("refit")
#}

#setwd(wd_refit)

ggplot(main_ancestry_dt) + geom_point(aes(distFromAnkara_km, weight.p3))

#create scatterplot of original data values
png(file="sum_log_likelihoodPredict.png", width = 6, height = 6, units = "in", res = 900)
plot(main_ancestry_dt_best_fit_quad$Learning_Prob_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16, xlab="Learning Rate", ylab="Sum Log Likelihood")
#add predicted lines based on quadratic regression model
lines(LearningValues, sum_log_likelihoodPredict, col='blue')
dev.off()

# Fit regression to the ancient DNA data:

#Standard lm                            
#empirical_slope_regression = lm(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)

#Robust regression package 1  
#empirical_slope_robust_regression = rlm(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)

#Robust regression package 2
empirical_slope_robust_regression2_northern_inds = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds_northern_inds)
empirical_slope_robust_regression2_southern_inds = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds_southern_inds)
empirical_slope_robust_regression2 = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)
empirical_slope_robust_regression2_early_inds = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds[time_period == "Early"])
empirical_slope_robust_regression2_middle_inds = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds[time_period == "Middle"])
empirical_slope_robust_regression2_late_inds = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds[time_period == "Late"])

lmrob_coefficients = empirical_slope_robust_regression2$coefficients
lmrob_coefficients = as.vector(lmrob_coefficients)
lmrob_intercept = lmrob_coefficients[1]
lmrob_slope = lmrob_coefficients[2]

lmrob_coefficients_northern_inds = empirical_slope_robust_regression2_northern_inds$coefficients
lmrob_coefficients_northern_inds = as.vector(lmrob_coefficients_northern_inds)
lmrob_intercept_northern_inds = lmrob_coefficients_northern_inds[1]
lmrob_slope_northern_inds = lmrob_coefficients_northern_inds[2]

lmrob_coefficients_southern_inds = empirical_slope_robust_regression2_southern_inds$coefficients
lmrob_coefficients_southern_inds = as.vector(lmrob_coefficients_southern_inds)
lmrob_intercept_southern_inds = lmrob_coefficients_southern_inds[1]
lmrob_slope_southern_inds = lmrob_coefficients_southern_inds[2]

lmrob_coefficients_early_inds = empirical_slope_robust_regression2_early_inds$coefficients
lmrob_coefficients_early_inds = as.vector(lmrob_coefficients_early_inds)
lmrob_intercept_early_inds = lmrob_coefficients_early_inds[1]
lmrob_slope_early_inds = lmrob_coefficients_early_inds[2]

lmrob_coefficients_middle_inds = empirical_slope_robust_regression2_middle_inds$coefficients
lmrob_coefficients_middle_inds = as.vector(lmrob_coefficients_middle_inds)
lmrob_intercept_middle_inds = lmrob_coefficients_middle_inds[1]
lmrob_slope_middle_inds = lmrob_coefficients_middle_inds[2]

lmrob_coefficients_late_inds = empirical_slope_robust_regression2_late_inds$coefficients
lmrob_coefficients_late_inds = as.vector(lmrob_coefficients_late_inds)
lmrob_intercept_late_inds = lmrob_coefficients_late_inds[1]
lmrob_slope_late_inds = lmrob_coefficients_late_inds[2]

formatted_prediction_with_conf
