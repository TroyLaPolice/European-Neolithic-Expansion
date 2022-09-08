# Load Libraries
library(data.table)
library(ggplot2)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/8-4")
map_size_km = 3700

square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)
square_file_names = mixedsort(square_file_names)
param_file_name = "initial_run_params_clean.txt"

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
                       .(Year, Scale_Factor, Movement_SD, Learning_Prob, Assortative_Mating)]

# Run a linear model to calculate the speed of wave and add as a column
all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Scale_Factor, Movement_SD, Learning_Prob, Assortative_Mating)]

# ----------------------------------------------------------------------------------------------------------------
# Plot to make visual checks
# ----------------------------------------------------------------------------------------------------------------

ggplot(all_sim_data[Learning_Prob == "Learning Prob = 0.0" & Assortative_Mating == "Assortative Mating = None"]) + 
  geom_point(aes(Year, km_50perc, col= Movement_SD)) + facet_grid(Scale_Factor ~ .)

ggplot(all_sim_data[Learning_Prob == "Learning Prob = 0.0" & Assortative_Mating == "Assortative Mating = None" & Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col= factor(Year))) + facet_grid(Scale_Factor ~ Movement_SD)

ggplot(all_sim_data) + 
  geom_point(aes(Movement_SD, speedOfWave, col= Learning_Prob)) + facet_grid(Scale_Factor ~ Assortative_Mating)


