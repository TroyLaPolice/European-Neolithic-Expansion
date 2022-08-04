# Load Libraries
library(ggplot2)
library(tidyr)
library(gtools)
library(data.table)

setwd("C:/Users/troyl/OneDrive/School/Documents/Grad School/Huber Lab/cluster_plots/8-1")

square_file_names = list.files(".", pattern="sim_sq*", full.names = TRUE)
square_file_names = mixedsort(square_file_names)
param_file_name = "initial_run_params_clean.txt"

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

# Plotting (examples)
no_assort_ratio = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = None"]) + geom_point(aes(Year, RatioFarmerToHG)) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "No Assortative Mating") + labs(y = "Population Ratio of Farmers to HGs")

full_assort_ratio = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = Full Assortative"]) + geom_point(aes(Year, RatioFarmerToHG)) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Full Assortative Mating") + labs(y = "Population Ratio of Farmers to HGs")

full_assort_parts = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = Full Assortative" & Year %% 50 == 0]) + geom_line(aes(as.numeric(Partition), Farmers_in_Partition, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Full Assortative Mating") + labs(x = "Partition") + labs(y = "Number of Farmers")

no_assort_parts = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = None" & Year %% 50 == 0]) + geom_line(aes(as.numeric(Partition), Farmers_in_Partition, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "No Assortative Mating") + labs(x = "Partition") + labs(y = "Number of Farmers")

full_assort_parts_ancestry = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = Full Assortative" & Year %% 50 == 0]) + geom_line(aes(as.numeric(Partition), Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Full Assortative Mating") + labs(x = "Partition") + labs(y = "Percent Farming Ancestry")

no_assort_parts_ancestry = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = None" & Year %% 50 == 0]) + geom_line(aes(as.numeric(Partition), Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "No Assortative Mating") + labs(x = "Partition") + labs(y = "Percent Farming Ancestry") + ylim(0,1.0)

no_assort_pop = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = None"]) + geom_point(colour="blue", aes(Year, TotalFarmers, col = "Farmers")) + geom_point(colour="red", aes(Year, TotalHGs, col = "Hunter Gatherers")) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "No Assortative Mating") + labs(y = "Population Size")

full_assort_pop = ggplot(all_sim_data[Assortative_Mating == "Assortative Mating = Full Assortative"]) + geom_point(colour="blue", aes(Year, TotalFarmers, col = "Farmers")) + geom_point(colour="red", aes(Year, TotalHGs, col = "Hunter Gatherers")) + 
  facet_grid(Movement_SD ~ Learning_Prob) + theme_bw() + labs(title = "Full Assortative Mating") + labs(y = "Population Size")

# Save plots
ggsave("full_assort_parts_out.png", plot = full_assort_parts, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("no_assort_parts_out.png", plot = no_assort_parts, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("no_assort_ratio_out.png", plot = no_assort_ratio, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("full_assort_ratio_out.png", plot = full_assort_ratio, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("full_assort_parts_ancestry_out.png", plot = full_assort_parts_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("no_assort_parts_ancestry_out.png", plot = no_assort_parts_ancestry, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("no_assort_pop_out.png", plot = no_assort_pop, units = "in", width = 10, height = 8, device="png", dpi=700)
ggsave("full_assort_pop_out.png", plot = full_assort_pop, units = "in", width = 10, height = 8, device="png", dpi=700)
