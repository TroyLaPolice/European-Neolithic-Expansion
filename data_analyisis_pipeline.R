# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/different_mvmt_ranges_5-10km_different_maps/square_runs")
map_size_km = 3700

# -------------------------------
# If landscape is simple square
# -------------------------------

# Specify and select files
square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)

if (length(square_file_names) != 0)
{
  # Sort Files
  square_file_names = mixedsort(square_file_names)
  # Read in File
  square_files = lapply(square_file_names, read.csv)
}

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
# If there's an ancestry distribution file
# -----------------------------------------

# Specify and select files
ancestry_dist_files_names = list.files(".", pattern="sim_ancestry_distribution*", full.names = TRUE)

if (length(ancestry_dist_files_names) != 0)
{
  # Sort Files
  ancestry_dist_files_names = mixedsort(ancestry_dist_files_names)
  # Read in File
  ancestry_dist_files = lapply(ancestry_dist_files_names, read.csv)
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
  ancestry_sample_files = lapply(ancestry_sample_names, read.csv)
}

# -----------------------------------------
# Read in parameter file
# -----------------------------------------

# Read in Param File
param_file_name = "param_inputs_clean.txt"
param_file = lapply(param_file_name, read.csv)

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to "long" format and combining different parameter combinations for the large data table
# ----------------------------------------------------------------------------------------------------------------

if (length(square_file_names) != 0){
  files = square_files}
if (length(non_square_file_names) != 0){
  files = non_square_files}

all_sim_data = lapply(1:length(files), function(x) 
  {
    params = strsplit(as.character(param_file[[1]][x,]), split = "   ")[[1]]
    params = matrix(params, nrow = 1)
    colnames(params) = gsub(sapply(strsplit(params, split = " = "), function(x) x[1]), pattern = " ", replacement = "_")
    data = files[[x]]
    
    if (length(square_file_names) != 0)
    {
    partitionedVariablesID = grepl("\\d", names(data))
    partitionedVariables = names(data)[partitionedVariablesID]
    partitionedVariablesUnique = unique(gsub(names(data)[partitionedVariablesID], pattern = "\\d+", replacement = ""))
    partitionNumber = as.numeric(gsub(names(data)[partitionedVariablesID], pattern = "[a-zA-Z_]", replacement = ""))
    numberOfPartitions = max(partitionNumber)
    
    separate_partitions_ID = rep(1:length(partitionedVariablesUnique), each = numberOfPartitions)
    separate_partitions_list = lapply(1:length(partitionedVariablesUnique), function(x) partitionedVariables[separate_partitions_ID == x])
    
    data_melted = data.table::melt(data = data.table(data), measure = separate_partitions_list, value.name = partitionedVariablesUnique, variable.name = "Partition")
    }
    
    if (length(non_square_file_names) != 0){
      data_melted = data}
    
    cbind(params, data_melted)

  }
  )
  
  # Create data table with simulation output data
  all_sim_data = rbindlist(all_sim_data)

  # Add columns that remove extra strings in parameter set
  all_sim_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
  all_sim_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
  all_sim_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
  all_sim_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
  
  colnames(all_sim_data)[colnames(all_sim_data) == "Farmer_Ancestry_Farmers"] = "Farmer_Ancestry_All_Farmers"

# ----------------------------------------------------------------------------------------------------------------
# Adding parameter combinations to the ancestry distribution table
# ----------------------------------------------------------------------------------------------------------------

if (length(ancestry_dist_files_names) != 0)
{
  ancestry_dist_data = lapply(1:length(ancestry_dist_files), function(x) 
  {
    params = strsplit(as.character(param_file[[1]][x,]), split = "   ")[[1]]
    params = matrix(params, nrow = 1)
    colnames(params) = gsub(sapply(strsplit(params, split = " = "), function(x) x[1]), pattern = " ", replacement = "_")
    data = ancestry_dist_files[[x]]
    
    cbind(params, data)
  }
  )
  
  # Create data table with simulation output data
  ancestry_dist_data = rbindlist(ancestry_dist_data)
  
  # Add columns that remove extra strings in parameter set
  ancestry_dist_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_dist_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_dist_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
  ancestry_dist_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
}

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
}

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to include km values for partitions
# ----------------------------------------------------------------------------------------------------------------

if (length(square_file_names) != 0)
{
  # Calculate the number of partitions used
  num_parts = max(as.integer(all_sim_data$Partition))

  # Calculate where the first partition midpoint point is
  partition_size = (map_size_km / num_parts)

  all_sim_data[, Mid_Point_km := as.numeric(Partition)*partition_size - partition_size/2]
}
  
# ----------------------------------------------------------------------------------------------------------------
# Function for approximating distance where 50% of the population is farmer
# ----------------------------------------------------------------------------------------------------------------

interpolate = function(km, value){
  
  if (all(value < 0.5) | all(value > 0.5)) return(NA)
  
  # Approximate a distance on the x-axis in km when the ratio of farmers to HGs is 50/50 
  # or the ancestry is at least 50% farmer (depending on what is passed to the function)
  interpolation = approx(value, km, xout=0.5, method="linear")
  
  # Pull out the distance from the approximation list
  approximated_x = interpolation[[2]]
  
  # Return value
  return(approximated_x)
  
}

# ----------------------------------------------------------------------------------------------------------------
# Run on all data and create new columns for distance and calculate speed
# ----------------------------------------------------------------------------------------------------------------

if (length(square_file_names) != 0)
{
  # Add column to data for the distance traveled on the x-axis in km when the ratio of farmers to HGs is 50/50
  all_sim_data[, km_50perc := as.numeric(0)]
  all_sim_data[, km_50perc := interpolate(Mid_Point_km, RatioFarmerToHG_Partition), 
               .(Year, Learning_Prob, Assortative_Mating, Movement)]
  
  # Run a linear model to calculate the speed of wave and add as a column
  all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Learning_Prob, Assortative_Mating, Movement)]
}
  
# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

# ----------------------------------------------------------------------------------------------------------------
# Plot data
# ----------------------------------------------------------------------------------------------------------------

remaining_ancestry = ggplot(all_sim_data[TotalHGs == 0]) + 
  geom_point(aes(as.numeric(Assortative_Mating_n), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + geom_line(aes(as.numeric(Assortative_Mating_n), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), 
                             group = factor(Learning_Prob_n))) + labs(x = "Percentage of Assortative Mating (0 = None, 1 = Full)") + 
  labs(y = "Remaining Farming Ancestry") + ggtitle("Remaining Farmer Ancestry", 
                                                   subtitle = "(When Zero Hunter Gatherers Remain)")
ggsave("remaining_ancestry.png", plot = remaining_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)

overall_ancestry_farmers = ggplot(all_sim_data) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating_n) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry") + ggtitle("Farming Ancestry x Time (Just Farming Population)", 
                                         subtitle = "Faceted by Assortative Mating Preference")
ggsave("overall_ancestry_farmers.png", plot = overall_ancestry_farmers, units = "in", width = 10, height = 11, device="png", dpi=700)

overall_ancestry_all = ggplot(all_sim_data) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating_n) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry") + ggtitle("Farming Ancestry x Time (Whole Population)", 
                                         subtitle = "Faceted by Assortative Mating Preference")
ggsave("overall_ancestry_all.png", plot = overall_ancestry_all, units = "in", width = 10, height = 11, device="png", dpi=700)

# -------------------------------
# If landscape is simple square
# -------------------------------

if (length(square_file_names) != 0)
{
  speed = ggplot(all_sim_data[Learning_Prob_n != 1 & Learning_Prob_n != 0.01 & Learning_Prob_n != 0.1]) + 
    geom_point(aes(as.numeric(Assortative_Mating_n), speedOfWave, col=factor(Learning_Prob_n), pch=factor(Movement_n))) + theme_bw() + 
    labs(y = "Speed of Wave (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 1 = Full)") +
    ggtitle("Speed of Wave (km per year)")
  ggsave("speed.png", plot = speed, units = "in", width = 8, height = 6, device="png", dpi=700)
  
  ancestry = 
    ggplot(all_sim_data[Year %% 200 == 0]) + 
    geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
    facet_grid(as.numeric(Assortative_Mating_n) ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
    labs(y = "Percent Farming Ancestry in Farmers") + ggtitle("Farming Ancestry x Distance from Origin Point (Just Farming Population)", 
                                                              subtitle = "Faceted by Learning Probability ~ Assortative Mating Preference")
  ggsave("farming_ancestry_in_farmers.png", plot = ancestry, units = "in", width = 14, height = 7, device="png", dpi=700)
  
  ancestry_whole_pop = 
    ggplot(all_sim_data[Year %% 200 == 0]) + 
    geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = Year, group = factor(Year))) + 
    facet_grid(as.numeric(Assortative_Mating_n) ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
    labs(y = "Percent Farming Ancestry in Population") + ggtitle("Farming Ancestry x Distance from Origin Point (Whole Population)", 
                                                                 subtitle = "Faceted by Learning Probability ~ Assortative Mating Preference")
  ggsave("farming_ancestry_in_whole_pop.png", plot = ancestry_whole_pop, units = "in", width = 14, height = 7, device="png", dpi=700)
  
  wave = ggplot(all_sim_data[Year %% 200 == 0]) + 
    geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col = Year, group = factor(Year))) + 
    facet_grid(Assortative_Mating_n ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
    labs(y = "Percent Farmers") + ggtitle("Percent Farmers (Behaviorally) x Distance from Origin Point", 
                                          subtitle = "Faceted by Learning Probability ~ Assortative Mating Preference")
  ggsave("wave.png", plot = wave, units = "in", width = 14, height = 8, device="png", dpi=700)
  
  wave_ancestry_overlay = ggplot(all_sim_data[Year %% 500 == 0]) + 
    geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col = "RatioFarmerToHG", group = factor(Year))) + 
    facet_grid(Assortative_Mating_n ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
    geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = "Farming Ancestry", group = factor(Year))) + 
    labs(y = "Percent Farmers") + ggtitle("Percent Farmers (Behaviorally) and Farming Ancestry x Distance from Origin Point", 
                                          subtitle = "Faceted by Learning Probability ~ Assortative Mating Preference")
  ggsave("wave_ancestry_overlay.png", plot = wave_ancestry_overlay, units = "in", width = 14, height = 8, device="png", dpi=700)
}

# -----------------------------------------
# If there's an ancestry sample file
# -----------------------------------------

if (length(ancestry_sample_names) != 0)
{
  assort_ancestry_by_x = ggplot(ancestry_sample_data[Learning_Prob_n == 0]) + geom_point(aes(Individual_X, Farming_Ancestry, col=Individual_Z)) + 
    facet_grid(cut(Year, seq(0,6000,1000))~Assortative_Mating) + ggtitle("Farming Ancestry by X Coordinate", 
                                                                        subtitle = "Faceted by Assortative Mating ~ 1,000 Year Time Bins (Learning Prob == 0)") + theme_bw() 
  ggsave("assort_ancestry_by_x.png", plot = assort_ancestry_by_x, units = "in", width = 20, height = 14, device="png", dpi=700)
  
  learning_ancestry_by_x = ggplot(ancestry_sample_data[Assortative_Mating_n == 1]) + geom_point(aes(Individual_X, Farming_Ancestry, col=Individual_Z)) + 
    facet_grid(cut(Year, seq(0,6000,1000))~Learning_Prob_n) + ggtitle("Farming Ancestry by X Coordinate", 
                                                                      subtitle = "Faceted by Learning Probability ~ 1,000 Year Time Bins (Assortative Mating == Full)") + theme_bw() 
  ggsave("learning_ancestry_by_x.png", plot = learning_ancestry_by_x, units = "in", width = 20, height = 14, device="png", dpi=700)
  
  assort_ancestry_2d = ggplot(ancestry_sample_data[Learning_Prob_n == 0]) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
    facet_grid(cut(Year, seq(0,6000,1000))~Assortative_Mating_n) + ggtitle("Farming Ancestry by X & Y Coordinate (Visualization of Ancestry Distribution on Landscape)", 
                                                                      subtitle = "Faceted by Assortative Mating ~ 1,000 Year Time Bins (Learning Prob == 0)") + theme_bw() 
  ggsave("assort_ancestry_2d.png", plot = assort_ancestry_2d, units = "in", width = 20, height = 14, device="png", dpi=700)
  
  learning_ancestry_2d = ggplot(ancestry_sample_data[Assortative_Mating_n == 1]) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
    facet_grid(cut(Year, seq(0,6000,1000))~Learning_Prob_n) + ggtitle("Farming Ancestry by X & Y Coordinate (Visualization of Ancestry Distribution on Landscape)", 
                                                                            subtitle = "Faceted by Learning Probability ~ 1,000 Year Time Bins (Assortative Mating == Full)") + theme_bw() 
  ggsave("learning_ancestry_2d.png", plot = learning_ancestry_2d, units = "in", width = 20, height = 14, device="png", dpi=700)
  
  singlular_ancestry_2d = ggplot(ancestry_sample_data[Learning_Prob_n == 0.0005]) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
    facet_grid(cut(Year, seq(0,6000,2000)) ~ .) + ggtitle("Ancestry Distribution on Landscape", subtitle = "2,000 Year Time Bins") + theme_bw() 
  ggsave("singlular_ancestry_2d.png", plot = singlular_ancestry_2d, units = "in", width = 6.85, height = 14, device="png", dpi=700)
}
