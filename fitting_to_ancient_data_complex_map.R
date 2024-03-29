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
  
  # Set < 0 Values to 0
  ancestry_sample_data$Farming_Ancestry[ancestry_sample_data$Farming_Ancestry < 0.0] = 0.0
  
  # Set > 1 Values to 1
  ancestry_sample_data$Farming_Ancestry[ancestry_sample_data$Farming_Ancestry > 1] = 1.0
}

# ----------------------------------------------------------------------------------------------------------------
# Adding distances from lower corner and calculating ancestry by distance
# ----------------------------------------------------------------------------------------------------------------

#  (x1, y1) is the bottom right hand corner of the map (i.e., farmer origin)
x1 = map_size_km
y1 = 0

# Number of bins
num_parts = 20

# Calculate the distance of each individual from the origin
# sqrt((x2 – x1)^2 + (y2 – y1)^2)
ancestry_sample_data = ancestry_sample_data %>% mutate(Distance_from_Origin = sqrt((Individual_X-x1)^2 + (Individual_Y-y1)^2))

#Max Distance From Origin Sampled 
max_dist = max(ancestry_sample_data$Distance_from_Origin)

possible_cut_bins = c(0,262,524,786,1048,1310,1572,1834,2096,2358,2620,2882,3144,3406,3668,3930,4192,4454,4716,4978,5240)

# Bin based on distance from origin
ancestry_sample_data = ancestry_sample_data %>% 
  mutate(Distance_Bin = cut(Distance_from_Origin, breaks=possible_cut_bins, labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))

# Calculate where the first bin midpoint point is
bin_size = (max_dist / num_parts)

# Create column with mid point of each bin
ancestry_sample_data = ancestry_sample_data %>% mutate(Mid_Point_km = (as.numeric(Distance_Bin)*bin_size - bin_size/2))

# Group by param combination and calculate mean ancestry in each bin for each param combo
ancestry_sample_data = ancestry_sample_data %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Mean_Ancestry = mean(Farming_Ancestry))

# Get variance of sampling in each bin
ancestry_sample_data = ancestry_sample_data %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Ancestry_Variance = var(Farming_Ancestry))

# Convert to DT
ancestry_sample_data_dt = setDT(ancestry_sample_data)


# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

map = ne_countries(scale = "medium", returnclass = "sf")

setwd("/path/to/empirical/data")

ancient_data = data.table(read.table("qpAdm_Fixed__EuroHG_AnatoliaN.txt", header = T))

headers = data.table(read.table("qpadm_headers.csv", header = T, sep=","))

ancient_data_filtered = ancient_data[GroupMeanLat < 55 & GroupMeanBP > 4500 & GroupMeanBP < 8000 & GroupNinds > 2]

# Remove greenland individual

ancient_data_filtered = ancient_data_filtered[GroupMeanLong>-20]

# Comp dist from Ankara (29, 41)

ancient_data_filtered[, distFromAnkara := distm(c(GroupMeanLong, GroupMeanLat), rev(c(39.93, 32.85)), fun = distHaversine), target]

ancient_data_filtered[, distFromAnkara_km := distFromAnkara*0.001]

ancient_data_filtered[, time_slice := cut(GroupMeanBP/1000, breaks = seq(12,1,-0.5))]

ancient_data_filtered_df = ancient_data_filtered %>% 
  mutate(Distance_Bin = cut(distFromAnkara_km, breaks=possible_cut_bins, labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)))

# Convert back to DT
ancient_data_filtered = setDT(ancient_data_filtered_df)

# Read in Param File
param_file = data.table(read.table("param_inputs_clean_filtered.txt", header = T, sep=","))

# ----------------------------------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------------------------------

interpolate = function(ancient_x_dist, sim_x_dist, sim_farming_ancestry){
  
  interpolation = approx(sim_x_dist, sim_farming_ancestry, xout=ancient_x_dist, method="linear")
  
  # Pull out the distance from the approximation list
  approximated_ancestry = interpolation[[2]]
  
  View(interpolation)
  
  # Return value
  return(approximated_ancestry)
  
}

main_ancestry_dt = headers

for (learning_param in param_file$Learning_Prob_n)
{
  # Convert back to DT
  ancient_data_filtered = setDT(ancient_data_filtered)
  
  individual_param_run = ancestry_sample_data_dt[Learning_Prob_n == learning_param]
  
  sim_farming_ancestry = (individual_param_run$Bin_Mean_Ancestry)
  sim_x_dist = individual_param_run$Mid_Point_km
  ancient_x_dist = ancient_data_filtered$distFromAnkara_km
  
  interpolated = interpolate(ancient_x_dist, sim_x_dist, sim_farming_ancestry)
  
  ancient_data_filtered[, approximated_ancestry := as.numeric(0)]
  ancient_data_filtered[, approximated_ancestry := interpolated]
  
  ancient_data_filtered[, distance_from_sim := as.numeric(0)]
  ancient_data_filtered[, distance_from_sim := (scaled_weight_Anatolia_N - approximated_ancestry)]
  
  ancient_data_filtered[, distance_from_sim_squared := as.numeric(0)]
  ancient_data_filtered[, distance_from_sim_squared := (distance_from_sim^2)]
  
  ancient_data_filtered[, mean_distance_from_sim_squared := as.numeric(0)]
  ancient_data_filtered[, mean_distance_from_sim_squared := mean(distance_from_sim_squared)]
  
  # Convert back to DT
  ancient_data_filtered = setDT(ancient_data_filtered)
  
  ancient_data_filtered = merge(ancient_data_filtered, unique(ancestry_sample_data_dt[Learning_Prob_n == learning_param, .(Distance_Bin, Bin_Ancestry_Variance)], by="Distance_Bin"), by = "Distance_Bin")
  
  ancient_data_filtered = subset(ancient_data_filtered, select = -c(Bin_Ancestry_Variance.x))
  colnames(ancient_data_filtered)[colnames(ancient_data_filtered) == "Bin_Ancestry_Variance.y"] <- "Bin_Ancestry_Variance"
  
  ancient_data_filtered[, total_sd := as.numeric(0)]
  #ancient_data_filtered[, total_sd := sqrt((ancient_data_filtered$se.p1^2) + ancient_data_filtered$Bin_Ancestry_Variance)]
  ancient_data_filtered[, total_sd := ancient_data_filtered$se.p1]
                        
  ancient_data_filtered[, log_likelihood := as.numeric(0)]
  ancient_data_filtered[, log_likelihood := (dnorm(scaled_weight_Anatolia_N, approximated_ancestry, total_sd, log = T))]
  
  ancient_data_filtered[, sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, sum_log_likelihood := sum(log_likelihood)]
  
  ancient_data_filtered[, exponential_sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, exponential_sum_log_likelihood := (exp(sum_log_likelihood))]
  
  ancient_data_filtered[, Learning_Prob_n := as.numeric(0)]
  ancient_data_filtered[, Learning_Prob_n := learning_param]
  
  main_ancestry_dt = rbind(main_ancestry_dt, ancient_data_filtered)
  
}

ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))

main_ancestry_dt_best_fit_quad = main_ancestry_dt[Learning_Prob_n < 0.0025]

quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)

l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

LearningValues = seq(0, 0.006, 0.0001)

sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))

# This just gets an approximate best fitting value- which of your tested params fits the best:
#           approximate_prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]

# To get the actual best fitting prediction put this in the console:
#           quadraticModel["coefficients"]

# Then the prediction will equal this equation. Get the Learning_Probs from the coefficients
# prediction = ((-Learning_Prob_n) / (2*[I(Learning_Prob_n^2)]))

setwd(wd)

#create scatterplot of original data values
pdf(file="sum_log_likelihoodPredict.pdf")
plot(main_ancestry_dt_best_fit_quad$Learning_Prob_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16, xlab="Learning Rate", ylab="Sum Log Likelihood")
#add predicted lines based on quadratic regression model
lines(LearningValues, sum_log_likelihoodPredict, col='blue')
dev.off()

# Fit regression to the ancient DNA data:

#Standard lm                            
empirical_slope_regression = lm(scaled_weight_Anatolia_N ~ I(distFromAnkara/1000), data=ancient_data_filtered)
#Robust regression package 1  
empirical_slope_robust_regression = rlm(scaled_weight_Anatolia_N ~ I(distFromAnkara/1000), data=ancient_data_filtered)
#Robust regression package 2  
empirical_slope_robust_regression2 = lmrob(scaled_weight_Anatolia_N ~ I(distFromAnkara/1000), data=ancient_data_filtered)

