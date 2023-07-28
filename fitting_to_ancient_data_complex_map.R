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

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in simulation files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/complex_map")
map_size_km = 3700

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
    head(read.csv(x, header = T), 10000)
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

# Bin based on distance from origin
ancestry_sample_data = ancestry_sample_data %>% mutate(Distance_Bin = ntile(Distance_from_Origin, n=num_parts))

# Calculate where the first partition midpoint point is
partition_size = (max_dist / num_parts)

# Create column with mid point of each bin
ancestry_sample_data = ancestry_sample_data %>% mutate(Mid_Point_km = (Distance_Bin*partition_size - partition_size/2))

# Group by param combination and calculate mean ancestry in each bin for each param combo
ancestry_sample_data = ancestry_sample_data %>% group_by(Learning_Prob_n, Distance_Bin) %>%
  mutate(Bin_Mean_Ancestry = mean(Farming_Ancestry))

ancestry_sample_data_dt = setDT(ancestry_sample_data)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

map = ne_countries(scale = "medium", returnclass = "sf")

setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/empirical_data")

ancient_data = data.table(read.table("qpAdm_Fixed__EuroHG_AnatoliaN.txt", header = T))

headers = data.table(read.table("qpadm_headers.csv", header = T, sep=","))

ancient_data_filtered = ancient_data[GroupMeanLat < 55 & GroupMeanBP > 4500 & GroupMeanBP < 8000 & GroupNinds > 2]

# Remove greenland individual

ancient_data_filtered = ancient_data_filtered[GroupMeanLong>-20]

# Comp dist from Ankara (29, 41)

ancient_data_filtered[, distFromAnkara := distm(c(GroupMeanLong, GroupMeanLat), rev(c(39.93, 32.85)), fun = distHaversine), target]

ancient_data_filtered[, distFromAnkara_km := distFromAnkara*0.001]

ancient_data_filtered[, time_slice := cut(GroupMeanBP/1000, breaks = seq(12,1,-0.5))]

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

for (param in param_file$Learning_Prob_n)
{
  individual_param_run = ancestry_sample_data_dt[Learning_Prob_n == param]
  
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
  
  #total_varience = sqrt((ancient_data_filtered$se.p1)^2 * nrow(ancient_data_filtered) + (var(ancestry_sample_data_dt$Bin_Mean_Ancestry)))
    
  ancient_data_filtered[, log_likelihood := as.numeric(0)]
  ancient_data_filtered[, log_likelihood := (dnorm(scaled_weight_Anatolia_N, approximated_ancestry, se.p1, log = T))]
  
  ancient_data_filtered[, sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, sum_log_likelihood := sum(log_likelihood)]
  
  ancient_data_filtered[, exponential_sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, exponential_sum_log_likelihood := (exp(sum_log_likelihood))]
  
  ancient_data_filtered[, Learning_Prob_n := as.numeric(0)]
  ancient_data_filtered[, Learning_Prob_n := param]
  
  main_ancestry_dt = rbind(main_ancestry_dt, ancient_data_filtered)
  
}

ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))

main_ancestry_dt_best_fit_quad = main_ancestry_dt[sum_log_likelihood > -20000]

quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)

l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

LearningValues = seq(0, 0.006, 0.0001)

sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))

prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]

setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/complex_map")

#create scatterplot of original data values
plot(main_ancestry_dt_best_fit_quad$Learning_Prob_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16)
#add predicted lines based on quadratic regression model
lines(LearningValues, sum_log_likelihoodPredict, col='blue')

# ----------------------------------------------------------------------------------------------------------------
# Plot data
# ----------------------------------------------------------------------------------------------------------------

remaining_ancestry_distance = ggplot(ancestry_sample_data) + 
  geom_line(aes(Mid_Point_km, Bin_Mean_Ancestry, col = factor(Learning_Prob_n), group=(Learning_Prob_n))) + 
  theme_bw() + labs(x = "Distance From Origin (km)") + labs(y = "Remaining Farming Ancestry") + ggtitle("Remaining Farmer Ancestry", 
                                                                                                        subtitle = "(When Zero Hunter Gatherers Remain)")
ggsave("remaining_ancestry_distance.png", plot = remaining_ancestry_distance, units = "in", width = 14, height = 8, device="png", dpi=700)

plot_ancestry_map = ggplot(ancestry_sample_data) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
  facet_grid(Learning_Prob_n ~ .) + ggtitle("Farming Ancestry by X & Y Coordinate (Visualization of Ancestry Distribution on Landscape)") + theme_bw() 
ggsave("plot_ancestry_map.png", plot = plot_ancestry_map, units = "in", width = 20, height = 14, device="png", dpi=700)

learning_empirical_overlay = ggplot() +
  geom_line(data = ancestry_sample_data, aes(Mid_Point_km, Bin_Mean_Ancestry, col = factor(Learning_Prob_n)), linewidth = 1.25) +
  geom_point(data = main_ancestry_dt, aes(distFromAnkara/1000, scaled_weight_Anatolia_N)) + 
  geom_smooth(n=1000, span=10) + coord_cartesian(ylim=c(0,1)) + theme_bw() + labs(col = "Learning Probability") +
  labs(y = "Remaining Farming Ancestry") + labs(x = "Distance From Origin (km)")

ggsave("learning_empirical_overlay.png", plot = learning_empirical_overlay, units = "in", width = 14, height = 8, device="png", dpi=700)
