# Load Libraries
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

# Set input params
wd = ("/path/to/simulation/data")

setwd(wd)

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
    binedVariablesID = grepl("\\d", names(data))
    binedVariables = names(data)[binedVariablesID]
    binedVariablesUnique = unique(gsub(names(data)[binedVariablesID], pattern = "\\d+", replacement = ""))
    binNumber = as.numeric(gsub(names(data)[binedVariablesID], pattern = "[a-zA-Z_]", replacement = ""))
    numberOfBins = max(binNumber)
    
    separate_bins_ID = rep(1:length(binedVariablesUnique), each = numberOfBins)
    separate_bins_list = lapply(1:length(binedVariablesUnique), function(x) binedVariables[separate_bins_ID == x])
    
    data_melted = data.table::melt(data = data.table(data), measure = separate_bins_list, value.name = binedVariablesUnique, variable.name = "Bin")
  }
  
  if (length(non_square_file_names) != 0){
    data_melted = data}
  
  cbind(params, data_melted)
  
}
)

# Create data table with simulation output data
all_sim_data = rbindlist(all_sim_data)

# Set > 1 Values to 1
all_sim_data$Farmer_Ancestry_All[all_sim_data$Farmer_Ancestry_All > 1] = 1
all_sim_data$Farmer_Ancestry_All_Farmers[all_sim_data$Farmer_Ancestry_All_Farmers > 1] = 1
all_sim_data$Farmer_Ancestry_All_HGs[all_sim_data$Farmer_Ancestry_All_HGs > 1] = 1
all_sim_data$Farmer_Ancestry_Bin_Farmers[all_sim_data$Farmer_Ancestry_Bin_Farmers > 1] = 1
all_sim_data$Farmer_Ancestry_Bin_HGs[all_sim_data$Farmer_Ancestry_Bin_HGs > 1] = 1
all_sim_data$Farmer_Ancestry_Bin_All[all_sim_data$Farmer_Ancestry_Bin_All > 1] = 1

# Add columns that remove extra strings in parameter set
#all_sim_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
all_sim_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
all_sim_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
all_sim_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
#all_sim_data[, Movement_X_n := as.numeric(gsub(Movement_X, pattern = ".*\\s(.+)", replacement = "\\1"))]
#all_sim_data[, Movement_Y_n := as.numeric(gsub(Movement_Y, pattern = ".*\\s(.+)", replacement = "\\1"))]
all_sim_data[, Map_Style_n := as.numeric(gsub(Map_Style, pattern = ".*\\s(.+)", replacement = "\\1"))]
#all_sim_data[, Water_Crossings_n := as.numeric(gsub(Water_Crossings, pattern = ".*\\s(.+)", replacement = "\\1"))]
#all_sim_data[, Replicate_n := as.numeric(gsub(Replicate, pattern = ".*\\s(.+)", replacement = "\\1"))]

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
    
    #cbind(params, data)
  }
  )
  
  # Create data table with simulation output data
  ancestry_dist_data = rbindlist(ancestry_dist_data)
  
  # Add columns that remove extra strings in parameter set
  #ancestry_dist_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Movement_X_n := as.numeric(gsub(Movement_X, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Movement_Y_n := as.numeric(gsub(Movement_Y, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Map_Style_n := as.numeric(gsub(Map_Style, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_dist_data[, Water_Crossings_n := as.numeric(gsub(Water_Crossings, pattern = ".*\\s(.+)", replacement = "\\1"))]
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
    
    #cbind(params, data)
  }
  )
  
  # Create data table with simulation output data
  ancestry_sample_data = rbindlist(ancestry_sample_data)
  
  # Add columns that remove extra strings in parameter set
  #ancestry_sample_data[, Downscale_n := as.numeric(gsub(Downscale, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Learning_Prob_n := as.numeric(gsub(Learning_Prob, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Assortative_Mating_n := as.numeric(gsub(Assortative_Mating, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Movement_n := as.numeric(gsub(Movement, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Movement_X_n := as.numeric(gsub(Movement_X, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Movement_Y_n := as.numeric(gsub(Movement_Y, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Map_Style_n := as.numeric(gsub(Map_Style, pattern = ".*\\s(.+)", replacement = "\\1"))]
  #ancestry_sample_data[, Water_Crossings_n := as.numeric(gsub(Water_Crossings, pattern = ".*\\s(.+)", replacement = "\\1"))]
}

# ----------------------------------------------------------------------------------------------------------------
# Transforming data to include km values for bins
# ----------------------------------------------------------------------------------------------------------------

if (length(square_file_names) != 0)
{
  # Calculate the number of bins used
  num_parts = max(as.integer(all_sim_data$Bin))
  
  # Calculate where the first bin midpoint point is
  bin_size = (map_size_km / num_parts)
  
  all_sim_data[, Mid_Point_km := as.numeric(Bin)*bin_size - bin_size/2]
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
  all_sim_data[, km_50perc := interpolate(Mid_Point_km, RatioFarmerToHG_Bin), 
               .(Year, Learning_Prob, Assortative_Mating, Movement)]
  
  # Run a linear model to calculate the speed of wave and add as a column
  all_sim_data[, speedOfWave := lm(km_50perc ~ Year)$coeff[2], .(Learning_Prob, Assortative_Mating, Movement)]
}

# ----------------------------------------------------------------------------------------------------------------
# Run on all data and create new columns for percentage of cultural effect
# ----------------------------------------------------------------------------------------------------------------
  
if (length(square_file_names) != 0)
{
  # Pull Out Speed of Fully Demic Model
  speed_fully_demic_cols = all_sim_data[Learning_Prob_n == 0.0 & Assortative_Mating_n == 1.0]
  speed_fully_demic = speed_fully_demic_cols$speedOfWave[1]
  
  # Add column for Cultural Effect
  all_sim_data[, cultural_effect := as.numeric(0)]
  
  # Calculate Cultural Effect
  all_sim_data[, cultural_effect := ((speedOfWave-speed_fully_demic)/speedOfWave)]

  # Set < 0 Values to 0
  all_sim_data$cultural_effect[all_sim_data$cultural_effect < 0.0] = 0.0
}

# Write to output file
fwrite(all_sim_data, file = "all_sim_data.csv", append = FALSE, quote = "auto", sep = ",")

first_complete_year = as.data.frame(all_sim_data)

first_complete_year = subset(first_complete_year, TotalHGs==0)
first_complete_year =  subset(first_complete_year, 
              subset = !duplicated(first_complete_year[c("TotalHGs","Learning_Prob_n", "Movement_n", "Assortative_Mating_n", "Mid_Point_km")]))
              
first_complete_year = setDT(first_complete_year)

setwd("/path/to/plotting/output/directory")

# ----------------------------------------------------------------------------------------------------------------
# Plot data
# ----------------------------------------------------------------------------------------------------------------

custom_palette <- c(
  "#3cb371",  # Medium Sea Green,
  "#0077ff",  # Electric Blue
  "#e41a1c",  # Vibrant Red
  "#ff9980",  # Light Coral Orange
  "#d95f02",  # Burnt Orange
  "#f9a8d4",  # Light Pastel Pink
  "#e7298a",  # Bright Magenta
  "#927fbf",  # Cool Lavender
  "#008080",  # Deep Teal
  "#99d6d6",  # Light Teal
  "#4b297b"   # Darker Eggplant
)
