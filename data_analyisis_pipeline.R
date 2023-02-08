# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/collab_runs/ancestry_dist_more_inds_allkm")
map_size_km = 3700

# Specify and select files
square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)
square_file_names = mixedsort(square_file_names)

ancestry_dist_files_names = list.files(".", pattern="sim_ancestry_distribution*", full.names = TRUE)
ancestry_dist_files_names = mixedsort(ancestry_dist_files_names)

ancestry_sample_names = list.files(".", pattern="sim_ancestry_sample*", full.names = TRUE)
ancestry_sample_names = mixedsort(ancestry_sample_names)

param_file_name = "param_inputs_clean.txt"

# Read in Files
square_files = lapply(square_file_names, read.csv)
ancestry_dist_files = lapply(ancestry_dist_files_names, read.csv)
ancestry_sample_files = lapply(ancestry_sample_names, read.csv)
param_file = lapply(param_file_name, read.csv)

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to "long" format and combining different parameter combinations for the large data table
# ----------------------------------------------------------------------------------------------------------------

all_sim_data = lapply(1:length(square_files), function(x) 
{
  params = strsplit(as.character(param_file[[1]][x,]), split = "   ")[[1]]
  params = matrix(params, nrow = 1)
  colnames(params) = gsub(sapply(strsplit(params, split = " = "), function(x) x[1]), pattern = " ", replacement = "_")
  data = square_files[[x]]
  
  partitionedVariablesID = grepl("\\d", names(data))
  partitionedVariables = names(data)[partitionedVariablesID]
  partitionedVariablesUnique = unique(gsub(names(data)[partitionedVariablesID], pattern = "\\d+", replacement = ""))
  partitionNumber = as.numeric(gsub(names(data)[partitionedVariablesID], pattern = "[a-zA-Z_]", replacement = ""))
  numberOfPartitions = max(partitionNumber)
  
  separate_partitions_ID = rep(1:length(partitionedVariablesUnique), each = numberOfPartitions)
  separate_partitions_list = lapply(1:length(partitionedVariablesUnique), function(x) partitionedVariables[separate_partitions_ID == x])
  
  data_melted = data.table::melt(data = data.table(data), measure = separate_partitions_list, value.name = partitionedVariablesUnique, variable.name = "Partition")
  
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

# ----------------------------------------------------------------------------------------------------------------
# Adding parameter combinations to the ancestry distribution table
# ----------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------
# Adding parameter combinations to the ancestry sample table
# ----------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to include km values for partitions
# ----------------------------------------------------------------------------------------------------------------

# Calculate the number of partitions used
num_parts = max(as.integer(all_sim_data$Partition))

# Calculate where the first partition midpoint point is
partition_size = (map_size_km / num_parts)

all_sim_data[, Mid_Point_km := as.numeric(Partition)*partition_size - partition_size/2]

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

# Add column to data for the distance traveled on the x-axis in km when the ratio of farmers to HGs is 50/50
all_sim_data[, km_50perc := as.numeric(0)]
all_sim_data[, km_50perc := interpolate(Mid_Point_km, RatioFarmerToHG_Partition), 
             .(Year, Learning_Prob, Assortative_Mating, Movement)]

# Run a linear model to calculate the speed of wave and add as a column
all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Learning_Prob, Assortative_Mating, Movement)]

# ----------------------------------------------------------------------------------------------------------------
# Run on all data and create new columns for distance and calculate ancestry cline
# ----------------------------------------------------------------------------------------------------------------

# Add column to data for the distance traveled on the x-axis in km when the farming ancestry is at least 50 percent
all_sim_data[, Ancestry_Cline := as.numeric(0)]
all_sim_data[, Ancestry_Cline := interpolate(Mid_Point_km, Farmer_Ancestry_Partition_All), 
             .(Year, Learning_Prob, Assortative_Mating, Movement)]

# Run a linear model to calculate the speed of wave and add as a column
all_sim_data[, speedOfAncestry := lm(Ancestry_Cline ~ Year)$coeff[2], .(Learning_Prob, Assortative_Mating, Movement)]

# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

# ----------------------------------------------------------------------------------------------------------------
# Plot to make visual checks
# ----------------------------------------------------------------------------------------------------------------

speed = ggplot(all_sim_data[Learning_Prob_n != 1 & Learning_Prob_n != 0.01 & Learning_Prob_n != 0.1]) + 
  geom_point(aes(as.numeric(Assortative_Mating_n), speedOfWave, col=factor(Learning_Prob_n), pch=factor(Movement_n))) + theme_bw() + 
  labs(y = "Speed of Wave (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")
ggsave("speed.png", plot = speed, units = "in", width = 8, height = 6, device="png", dpi=700)

ancestry_speed = ggplot(all_sim_data) + 
  geom_point(aes(as.numeric(Assortative_Mating_n), speedOfAncestry, col=factor(Learning_Prob_n), pch=factor(Movement_n))) + theme_bw() + 
  labs(y = "Speed of Ancestry Expansion (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")
ggsave("ancestry_speed.png", plot = ancestry_speed, units = "in", width = 8, height = 6, device="png", dpi=700)

ancestry = 
  ggplot(all_sim_data[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(as.numeric(Assortative_Mating_n) ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry in Farmers")
ggsave("ancestry.png", plot = ancestry, units = "in", width = 11, height = 7, device="png", dpi=700)

ancestry_whole_pop = 
  ggplot(all_sim_data[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = Year, group = factor(Year))) + 
  facet_grid(as.numeric(Assortative_Mating_n) ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry in Population")
ggsave("ancestry_whole_pop.png", plot = ancestry_whole_pop, units = "in", width = 11, height = 7, device="png", dpi=700)

wave = ggplot(all_sim_data[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col = Year, group = factor(Year))) + 
  facet_grid(Assortative_Mating_n ~ Learning_Prob_n) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farmers")
ggsave("wave.png", plot = wave, units = "in", width = 10, height = 8, device="png", dpi=700)

remaining_ancestry = ggplot(all_sim_data[TotalHGs == 0]) + 
  geom_point(aes(as.numeric(Assortative_Mating_n), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + geom_line(aes(as.numeric(Assortative_Mating_n), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)") + 
  labs(y = "Remaining Farming Ancestry")
ggsave("remaining_ancestry.png", plot = remaining_ancestry, units = "in", width = 5, height = 4, device="png", dpi=700)

overall_ancestry_farmers = ggplot(all_sim_data) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating_n) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry")
ggsave("overall_ancestry_farmers.png", plot = overall_ancestry_farmers, units = "in", width = 9, height = 6, device="png", dpi=700)

overall_ancestry_all = ggplot(all_sim_data) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All, col = factor(Learning_Prob_n), group = factor(Learning_Prob_n))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating_n) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry")
ggsave("overall_ancestry_all.png", plot = overall_ancestry_all, units = "in", width = 9, height = 6, device="png", dpi=700)

overall_ancestry_all2 = ggplot(all_sim_data[Learning_Prob_n == 0.0]) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All, col = factor(as.numeric(Assortative_Mating_n)), group = factor(as.numeric(Assortative_Mating_n)))) + 
  theme_bw() + labs(x = "Year") + 
  labs(y = "Farming Ancestry")
ggsave("overall_ancestry_all2.png", plot = overall_ancestry_all2, units = "in", width = 9, height = 6, device="png", dpi=700)

assort_ancestry_by_x = ggplot(ancestry_sample_data[Learning_Prob_n == 0]) + geom_point(aes(Individual_X, Farming_Ancestry, col=Individual_Z)) + 
  facet_grid(cut(Year, seq(0,6000,1000))~Assortative_Mating)
ggsave("assort_ancestry_by_x.png", plot = assort_ancestry_by_x, units = "in", width = 21, height = 11, device="png", dpi=700)

learning_ancestry_by_x = ggplot(ancestry_sample_data[Assortative_Mating_n == 1]) + geom_point(aes(Individual_X, Farming_Ancestry, col=Individual_Z)) + 
  facet_grid(cut(Year, seq(0,6000,1000))~Learning_Prob_n)
ggsave("learning_ancestry_by_x.png", plot = learning_ancestry_by_x, units = "in", width = 21, height = 11, device="png", dpi=700)

assort_ancestry_2d = ggplot(ancestry_sample_data[Assortative_Mating_n == 1]) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
  facet_grid(cut(Year, seq(0,6000,1000))~Learning_Prob_n)
ggsave("assort_ancestry_2d.png", plot = assort_ancestry_2d, units = "in", width = 21, height = 11, device="png", dpi=700)

learning_ancestry_2d = ggplot(ancestry_sample_data[Learning_Prob_n == 0]) + geom_point(aes(Individual_X, Individual_Y, col = Farming_Ancestry, pch=factor(Individual_Z))) + 
  facet_grid(cut(Year, seq(0,6000,1000))~Assortative_Mating_n)
ggsave("learning_ancestry_2d.png", plot = learning_ancestry_2d, units = "in", width = 21, height = 11, device="png", dpi=700)
