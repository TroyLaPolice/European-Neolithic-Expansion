# Load Libraries
library(tidyr)
library(gtools)
library(data.table)
library(dplyr)

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/scal_fact_test")
num_parts = 20
map_size_km = 3700

square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)
square_file_names = mixedsort(square_file_names)
param_file_name = "initial_run_params_clean.txt"

# Read in Files
square_files = lapply(square_file_names, read.csv)
param_file = lapply(param_file_name, read.csv)


# Transforming data to "long" format and combining different parameter combinations
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

# Transforming data to include km values for partitions

# Add exta column with partition numbers to be replaced with km values
all_sim_data_km_values = cbind(all_sim_data, all_sim_data$Partition)
colnames(all_sim_data_km_values)[21] = "Mid_Point_km"

# Calculate where the first partition midpoint point is
first_partition_midpoint = ((map_size_km / num_parts) / 2)

#Set this to be the previous midpoint before calculating further in the loop
prev_midpt = first_partition_midpoint

# Initialize list of midpoints starting with first partition
midpoints = list(first_partition_midpoint)

# Loop through partitions appending midpoints
for (partition in 1:(num_parts-1)) {
  partition_midpoint = first_partition_midpoint * 2 + prev_midpt
  midpoints = append(midpoints, partition_midpoint)
  prev_midpt = partition_midpoint
}

# Convert entries in DF to character from factor for replacement
all_sim_data_km_values$Mid_Point_km = as.character(all_sim_data_km_values$Mid_Point_km)

# Replace partition numbers with corresponding km vlaues for midpoint
for (partition in 1:(num_parts)) {
  all_sim_data_km_values$Mid_Point_km[all_sim_data_km_values$Mid_Point_km == partition] = midpoints[partition]
}

# Turn back to factor
all_sim_data_km_values$Mid_Point_km = as.factor(all_sim_data_km_values$Mid_Point_km)
