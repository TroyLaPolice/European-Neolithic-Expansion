# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/9-12")
map_size_km = 3700

square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)
square_file_names = mixedsort(square_file_names)
param_file_name = "param_inputs_clean.txt"

# Read in Files
square_files = lapply(square_file_names, read.csv)
param_file = lapply(param_file_name, read.csv)

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to "long" format and combining different parameter combinations
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

all_sim_data = rbindlist(all_sim_data)

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

interpolate = function(km, ratio){
  
  if (all(ratio < 0.5) | all(ratio > 0.5)) return(NA)
  
  # Approximate a distance on the x-axis in km when the ratio of farmers to HGs is 50/50
  interpolation = approx(ratio, km, xout=0.5, method="linear")
  
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
all_sim_data[, km_50perc := interpolate(km = Mid_Point_km, ratio = RatioFarmerToHG_Partition), 
                       .(Year, Scale_Factor, Downscale, Movement_SD, Learning_Prob, Assortative_Mating)]

# Run a linear model to calculate the speed of wave and add as a column
all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Downscale, Scale_Factor, Movement_SD, Learning_Prob, Assortative_Mating)]

# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

# Simplify to pull out speeds for each param
sim_data_simple = all_sim_data[, !c("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmerToHG", "Farmer_Ancestry_All", "Farmer_Ancestry_All_Farmers", "Farmer_Ancestry_All_HGs", "Num_Repro_Age_Inds", "NewBirths", "ReproFreq", "Partition", "All_in_Partition", "Farmers_in_Partition", "HGs_in_Partition", "RatioFarmerToHG_Partition", "Farmer_Ancestry_Partition_Farmers", "Farmer_Ancestry_Partition_HGs", "Farmer_Ancestry_Partition_All", "Mid_Point_km", "km_50perc")]
sim_data_simple = unique(sim_data_simple, by = "speedOfWave")

# Write to output file
fwrite(sim_data_simple, file = "sim_data_simple.csv", append = FALSE, quote = "auto", sep = ",")

# ----------------------------------------------------------------------------------------------------------------
# Plot to make visual checks
# ----------------------------------------------------------------------------------------------------------------

plot1 = ggplot(all_sim_data) + 
  geom_point(aes(Year, km_50perc, col= Learning_Prob)) + facet_grid(Scale_Factor ~ Downscale)

plot2 = ggplot(all_sim_data[Learning_Prob == "Learning_Prob = 0.01" & Movement_SD == "Movement_SD = 10" & Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col= factor(Year))) + facet_grid(Scale_Factor ~ Downscale)

plot3 = ggplot(all_sim_data) + 
  geom_point(aes(Movement_SD, speedOfWave, col= Movement_SD)) + facet_grid(Scale_Factor ~ Downscale)

plot5 = ggplot(all_sim_data) + 
  geom_point(aes(Learning_Prob, speedOfWave, col= Learning_Prob)) + facet_grid(Scale_Factor ~ Downscale)

no_assort_ancestry = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = None" & Year %% 50 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "No Assortative Mating") + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry")

some_assort_ancestry = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = 0.8" & Year %% 50 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Some Assortative Mating") + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry")

full_assort_ancestry = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = Full Assortative" & Year %% 50 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Full Assortative Mating") + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry")

ggsave("plot1.png", plot = plot1, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot2.png", plot = plot2, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot3.png", plot = plot3, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot4.png", plot = no_assort_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot5.png", plot = plot5, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot6.png", plot = some_assort_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("plot7.png", plot = full_assort_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)
