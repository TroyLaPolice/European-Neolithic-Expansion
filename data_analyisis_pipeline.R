# Load Libraries
library(data.table)

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

# Add extra column with partition numbers to be replaced with km values
all_sim_data_km_values = cbind(all_sim_data, all_sim_data$Partition)
names(all_sim_data_km_values)[names(all_sim_data_km_values) == "V2"] = "Mid_Point_km"

# Calculate the number of partitions used
num_parts = max(as.integer(all_sim_data$Partition))

# Calculate where the first partition midpoint point is
first_partition_midpoint = ((map_size_km / num_parts) / 2)

#Set this to be the previous midpoint before calculating further in the loop
prev_midpt = first_partition_midpoint

# Initialize list of midpoints starting with first partition
midpoints = first_partition_midpoint

# Loop through partitions appending midpoints
for (partition in 1:(num_parts-1)) {
  partition_midpoint = first_partition_midpoint * 2 + prev_midpt
  midpoints = append(midpoints, partition_midpoint)
  prev_midpt = partition_midpoint
}

# Convert entries in DF to character from factor for replacement
all_sim_data_km_values$Mid_Point_km = as.character(all_sim_data_km_values$Mid_Point_km)

# Replace partition numbers with corresponding km values for midpoint
for (partition in 1:(num_parts)) {
  all_sim_data_km_values$Mid_Point_km[all_sim_data_km_values$Mid_Point_km == partition] = midpoints[partition]
}

# Turn back to factor
all_sim_data_km_values$Mid_Point_km = as.factor(all_sim_data_km_values$Mid_Point_km)

# ----------------------------------------------------------------------------------------------------------------

test_year_df = all_sim_data_km_values[Year == 600 & 
                                        Scale_Factor == "Scale Factor = 0.5" & 
                                        Movement_SD == "Movement SD = 10km" & 
                                        Learning_Prob == "Learning Prob = 0.0" & 
                                        Assortative_Mating == "Assortative Mating = None"]

test_year_df_simple = data.table(as.character(test_year_df$Mid_Point_km), as.character(test_year_df$RatioFarmerToHG_Partition))
colnames(test_year_df_simple) = c("km", "ratio")

# ----------------------------------------------------------------------------------------------------------------

interpolate = function(data){
  
  # Approximate a distance from the x-axis in km when the ratio of farmers to HGs is 50/50
  interpolation = approx(test_year_df_simple$ratio, test_year_df_simple$km, xout=0.5, method="linear")

  # Pull out the distance from the approximation list
  approximated_x = interpolation[[2]]

  plot(test_year_df_simple, pch = 19)
  points(approximated_x, 0.5, col='red', pch=19)
  
  # Return value
  return(approximated_x)
  
}

interpolate(test_year_df_simple)
