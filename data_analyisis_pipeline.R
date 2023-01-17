# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/cluster_runs/collab_runs/equation_model_param_match")
map_size_km = 3700

# Specify and select files
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

# Create data table with simulation output data
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
             .(Year, Replicate, Competition, Learning_Prob, Assortative_Mating)]

# Run a linear model to calculate the speed of wave and add as a column
all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Replicate, Competition, Learning_Prob, Assortative_Mating)]

# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

# ----------------------------------------------------------------------------------------------------------------
# Run on all data and create new columns for distance and calculate ancestry cline
# ----------------------------------------------------------------------------------------------------------------

# Add column to data for the distance traveled on the x-axis in km when the farming ancestry is at least 50 percent
all_sim_data[, Ancestry_Cline := as.numeric(0)]
all_sim_data[, Ancestry_Cline := interpolate(Mid_Point_km, Farmer_Ancestry_Partition_All), 
             .(Year, Replicate, Competition, Learning_Prob, Assortative_Mating)]

# Run a linear model to calculate the speed of wave and add as a column
#all_sim_data[, speedOfAncestry := lm(Ancestry_Cline ~ Year)$coeff[2], .(Replicate, Competition, Learning_Prob, Assortative_Mating)]

# ----------------------------------------------------------------------------------------------------------------
# Simplify and write outputs
# ----------------------------------------------------------------------------------------------------------------

# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

# Simplify to pull out speeds for each param
sim_data_simple = all_sim_data[, !c("Year", "PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmerToHG", "Farmer_Ancestry_All", "Farmer_Ancestry_All_HGs", "Num_Repro_Age_Inds", "NewBirths", "ReproFreq", "Partition", "All_in_Partition", "Farmers_in_Partition", "HGs_in_Partition", "Farmer_Ancestry_Partition_Farmers", "Farmer_Ancestry_Partition_HGs", "Farmer_Ancestry_Partition_All", "Mid_Point_km", "km_50perc", "Ancestry_Cline")]
sim_data_simple = unique(sim_data_simple, by = "speedOfWave")

# Write to output file
fwrite(sim_data_simple, file = "sim_data_simple.csv", append = FALSE, quote = "auto", sep = ",")

# Simplify to pull out speeds for each param
sim_data_simple = all_sim_data[, !c("PopulationSize", "TotalFarmers", "RatioFarmerToHG", "Farmer_Ancestry_All_HGs", "Num_Repro_Age_Inds", "NewBirths", "ReproFreq", "Partition", "All_in_Partition", "Farmers_in_Partition", "HGs_in_Partition", "Farmer_Ancestry_Partition_HGs", "km_50perc", "Ancestry_Cline")]

sim_data_simple_as_numeric = sim_data_simple %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 1.0", 100))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.99999", 99.999))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.999999", 99.9999))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.0", 0.0))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.8", 80))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.9", 90))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.95", 95))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.99", 99))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.999", 99.9))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Assortative_Mating = replace(Assortative_Mating, Assortative_Mating == "Assortative_Mating = 0.9999", 99.99))

sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0", 0.0))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.01", 0.01))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.1", 0.1))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 1", 1.0))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.000001", 0.000001))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.00001", 0.00001))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0001", 0.0001))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.001", 0.001))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.00075", 0.00075))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0009", 0.0009))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.00075", 0.00075))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0001", 0.0001))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.00005", 0.00005))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0015", 0.0015))
sim_data_simple_as_numeric = sim_data_simple_as_numeric %>% 
  mutate(Learning_Prob = replace(Learning_Prob, Learning_Prob == "Learning_Prob = 0.0025", 0.0025))

# Write to output file
fwrite(sim_data_simple, file = "sim_data_simple.csv", append = FALSE, quote = "auto", sep = ",")
# Write to output file
fwrite(sim_data_simple_as_numeric, file = "sim_data_simple_as_numeric.csv", append = FALSE, quote = "auto", sep = ",")

# ----------------------------------------------------------------------------------------------------------------
# Plot to make visual checks
# ----------------------------------------------------------------------------------------------------------------

ancestry = 
  ggplot(sim_data_simple_as_numeric[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_Farmers, col = Year, group = factor(Year))) + 
  facet_grid(as.numeric(Assortative_Mating) ~ Learning_Prob) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry in Farmers")

ggsave("ancestry.png", plot = ancestry, units = "in", width = 11, height = 7, device="png", dpi=700)

ancestry_whole_pop = 
  ggplot(sim_data_simple_as_numeric[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = Year, group = factor(Year))) + 
  facet_grid(as.numeric(Assortative_Mating) ~ Learning_Prob) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farming Ancestry in Population")

ggsave("ancestry_whole_pop.png", plot = ancestry_whole_pop, units = "in", width = 11, height = 7, device="png", dpi=700)

ancestry_speed_all_comp = ggplot(sim_data_simple_as_numeric[Competition != "Competition = Group"]) + 
  geom_point(aes(as.numeric(Assortative_Mating), speedOfAncestry, col=factor(Learning_Prob))) + theme_bw() + 
  labs(y = "Speed of Ancestry Expansion (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")

#ggsave("ancestry_speed.png", plot = ancestry_speed_all_comp, units = "in", width = 8, height = 6, device="png", dpi=700)

# ancestry_speed_group_comp = ggplot(sim_data_simple_as_numeric[Competition != "Competition = All"]) + 
#   geom_point(aes(as.numeric(Assortative_Mating), speedOfAncestry, col=factor(Learning_Prob))) + theme_bw() + 
#   labs(y = "Speed of Ancestry Expansion (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")

# ggsave("ancestry_speed_group_comp.png", plot = ancestry_speed_group_comp, units = "in", width = 8, height = 6, device="png", dpi=700)

speed = ggplot(sim_data_simple_as_numeric) + 
  geom_point(aes(as.numeric(Assortative_Mating), speedOfWave, col=factor(Learning_Prob))) + theme_bw() + 
  labs(y = "Speed of Wave (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")

ggsave("speed.png", plot = speed, units = "in", width = 8, height = 6, device="png", dpi=700)

# speed_comp = ggplot(sim_data_simple_as_numeric) + 
#   geom_point(aes(as.numeric(Assortative_Mating), speedOfWave, col=Competition)) + theme_bw() + 
#   labs(y = "Speed of Wave (km per year)") + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)")
# 
# ggsave("speed_comp.png", plot = speed_comp, units = "in", width = 8, height = 6, device="png", dpi=700)

# speed_reps_group = ggplot(sim_data_simple_as_numeric[Competition != "Competition = All" & Assortative_Mating == 99.999 & Learning_Prob != 0.001]) + 
#   geom_point(aes(as.numeric(Learning_Prob), speedOfWave, col=Replicate)) + theme_bw() + 
#   labs(y = "Speed of Wave (km per year)") + labs(x = "Learning Probability")
# 
# ggsave("speed_reps_group.png", plot = speed_reps_group, units = "in", width = 8, height = 6, device="png", dpi=700)
# 
# speed_reps_all = ggplot(sim_data_simple_as_numeric[Competition != "Competition = Group" & Assortative_Mating == 99.999 & Learning_Prob != 0.001]) + 
#   geom_point(aes(as.numeric(Learning_Prob), speedOfWave, col=Replicate)) + theme_bw() + 
#   labs(y = "Speed of Wave (km per year)") + labs(x = "Learning Probability")
# 
# ggsave("speed_reps_all.png", plot = speed_reps_all, units = "in", width = 8, height = 6, device="png", dpi=700)

wave = ggplot(sim_data_simple_as_numeric[Year %% 200 == 0]) + 
  geom_line(aes(Mid_Point_km, RatioFarmerToHG_Partition, col = Year, group = factor(Year))) + 
  facet_grid(Assortative_Mating ~ Learning_Prob) + theme_bw() + labs(x = "Distance From Origin (km)") + 
  labs(y = "Percent Farmers")

ggsave("wave.png", plot = wave, units = "in", width = 10, height = 8, device="png", dpi=700)

remaining_ancestry = ggplot(sim_data_simple_as_numeric[TotalHGs == 0]) + 
  geom_point(aes(as.numeric(Assortative_Mating), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob), group = factor(Learning_Prob))) + 
  theme_bw() + geom_line(aes(as.numeric(Assortative_Mating), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob), group = factor(Learning_Prob))) + labs(x = "Percentage of Assortative Mating (0 = None, 100 = Full)") + 
  labs(y = "Remaining Farming Ancestry")

ggsave("remaining_ancestry.png", plot = remaining_ancestry, units = "in", width = 5, height = 4, device="png", dpi=700)

overall_ancestry_farmers = ggplot(sim_data_simple_as_numeric) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All_Farmers, col = factor(Learning_Prob), group = factor(Learning_Prob))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry")

ggsave("overall_ancestry_farmers.png", plot = overall_ancestry_farmers, units = "in", width = 9, height = 6, device="png", dpi=700)

overall_ancestry_all = ggplot(sim_data_simple_as_numeric) + 
  geom_line(aes(as.numeric(Year), Farmer_Ancestry_All, col = factor(Learning_Prob), group = factor(Learning_Prob))) + 
  theme_bw() + facet_grid(as.numeric(Assortative_Mating) ~ .) + labs(x = "Year") + 
  labs(y = "Farming Ancestry")

ggsave("overall_ancestry_all.png", plot = overall_ancestry_all, units = "in", width = 9, height = 6, device="png", dpi=700)


