# Load Libraries
library(ggplot2)
library(tidyr)

# ********************************************************************************************************************************************
# Set WD, input and output file name additions, as well as trimming range for non-informative data
# ********************************************************************************************************************************************

setwd("/home/tml5905/Documents/HuberLab/HunterGatherFarmerInteractions/longrun_test/7-8")

input_file_extention = "_test2"  # Should be the same as the ending of your SLiM output files
                                #   if endings are default this will be an empty string.

output_file_extention = "_test2" # Addition to attach to the output files
                                #   Best practice is to keep these the same for simplicity.

trimming_range = -1

# Create file names
square_input_name = paste("sim_square_wave_stats_per_year", input_file_extention, ".txt", sep = "")
general_input_name = paste("sim_pop_stats_per_year", input_file_extention, ".txt", sep = "")
general_input_with_deaths_name = paste("sim_pop_stats_per_year_with_death", input_file_extention, ".csv", sep = "")

# Check if files exist and read in those that do
if ((file.exists(square_input_name)) && (file.exists(general_input_with_deaths_name))){
  square_input_file = read.csv(square_input_name, sep = " ", header = TRUE)
  general_input_file = read.csv(general_input_with_deaths_name, sep = ",", header = TRUE)
} else if ((file.exists(square_input_name)) && (file.exists(general_input_name))){
  square_input_file = read.csv(square_input_name, sep = " ", header = TRUE)
  general_input_file = read.csv(general_input_name, sep = " ", header = TRUE)
} else if (file.exists(square_input_name)){
  square_input_file = read.csv(square_input_name, sep = " ", header = TRUE)
} else if (file.exists(general_input_with_deaths_name)){
  general_input_file = read.csv(general_input_with_deaths_name, sep = ",", header = TRUE)
} else if (file.exists(general_input_name)){
  general_input_file = read.csv(general_input_name, sep = " ", header = TRUE)
} else{
  print("Requested files do not exist")
}

# ************************************************************************************************
# Melt the wide file dfs into long format
# ************************************************************************************************

# If run on a square
if (file.exists(square_input_name)){
square_input_file_tidyr = tidyr::pivot_longer(square_input_file, cols=c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                                                                        "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                                                                        "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                                                                        "Farmers_in_Partition10", "RatioFarmerToHG_Partition1", "HGs_in_Partition1", 
                                                                        "HGs_in_Partition2", "HGs_in_Partition3", "HGs_in_Partition4", 
                                                                        "HGs_in_Partition5", "HGs_in_Partition6", "HGs_in_Partition7", 
                                                                        "HGs_in_Partition8", "HGs_in_Partition9", "HGs_in_Partition10",
                                                                        "RatioFarmerToHG_Partition2", "RatioFarmerToHG_Partition3", "RatioFarmerToHG_Partition4", 
                                                                        "RatioFarmerToHG_Partition5",  "RatioFarmerToHG_Partition6", "RatioFarmerToHG_Partition7", 
                                                                        "RatioFarmerToHG_Partition8", "RatioFarmerToHG_Partition9", "RatioFarmerToHG_Partition10"),
                                      names_to = "Partition Variable")
}

# If for the general model
if (file.exists(general_input_name)){
  general_input_file_tidyr = tidyr::pivot_longer(general_input_file, cols=c("PopulationSize", "TotalFarmers", "TotalHGs", "RatioFarmertoHG"),
                                                names_to = "Variable")
}

# If the model was run on a square, plot the data
if (file.exists(square_input_name)){
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot
  # ************************************************************************************************
  
  # Push ratio data to its own data frame
  ratio_dat = data.frame(square_input_file$Year,
                         square_input_file$RatioFarmerToHG, square_input_file$RatioFarmerToHG_Partition1,
                         square_input_file$RatioFarmerToHG_Partition2, square_input_file$RatioFarmerToHG_Partition3, 
                         square_input_file$RatioFarmerToHG_Partition4, square_input_file$RatioFarmerToHG_Partition5, 
                         square_input_file$RatioFarmerToHG_Partition6, square_input_file$RatioFarmerToHG_Partition7, 
                         square_input_file$RatioFarmerToHG_Partition8, square_input_file$RatioFarmerToHG_Partition9, 
                         square_input_file$RatioFarmerToHG_Partition10)
  
  
  # Name columns
  colnames(ratio_dat) = c("Year", "RatioFarmerToHG_Total", "RatioFarmerToHG_Partition1", "RatioFarmerToHG_Partition2", 
                           "RatioFarmerToHG_Partition3", "RatioFarmerToHG_Partition4", "RatioFarmerToHG_Partition5", 
                           "RatioFarmerToHG_Partition6", "RatioFarmerToHG_Partition7", "RatioFarmerToHG_Partition8", 
                           "RatioFarmerToHG_Partition9", "RatioFarmerToHG_PartitionTen")
  
  
  # Remove non-informative ratio data
  head_ratio_dat = head(ratio_dat, trimming_range)
  
  # Convert to tidy format
  ratio_dat_tidyr = tidyr::pivot_longer(head_ratio_dat, cols=c("RatioFarmerToHG_Total", "RatioFarmerToHG_Partition1", "RatioFarmerToHG_Partition2", 
                                                               "RatioFarmerToHG_Partition3", "RatioFarmerToHG_Partition4", "RatioFarmerToHG_Partition5", 
                                                               "RatioFarmerToHG_Partition6", "RatioFarmerToHG_Partition7", "RatioFarmerToHG_Partition8", 
                                                               "RatioFarmerToHG_Partition9", "RatioFarmerToHG_PartitionTen"),
                                        names_to = "variable")
  
  # Plot
  ratio_plot = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=variable)) + geom_col(width = 1)
  
  # Create output file name
  ratio_plot_out = paste("ratio_partition_plot", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out, plot = ratio_plot, units = "in", width = 10, height = 5, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot2
  # ************************************************************************************************
  
  ratio_plot2 = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=variable)) + geom_col(width = 1) + facet_grid(variable~.)
  
  # Create output file name
  ratio_plot_out2 = paste("ratio_partition_plot2", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out2, plot = ratio_plot2, units = "in", width = 10, height = 12, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot3
  # ************************************************************************************************
  
  ratio_plot3 = ggplot(ratio_dat_tidyr[ratio_dat_tidyr$Year%%100 == 0 & ratio_dat_tidyr$variable != "RatioFarmerToHG_Total",], aes(x=variable, y=value, col=factor(Year), group=Year)) + geom_line() 
  
  # Create output file name
  ratio_plot_out3 = paste("ratio_partition_plot3", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out3, plot = ratio_plot3, units = "in", width = 10, height = 5, device="tiff", dpi=1000)
  
    
  # ************************************************************************************************
  # Farmer Population Size Partition Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  population_partition_farmer_dat = data.frame(square_input_file$Farmers_in_Partition1, square_input_file$Farmers_in_Partition2,
                                   square_input_file$Farmers_in_Partition3, square_input_file$Farmers_in_Partition4, 
                                   square_input_file$Farmers_in_Partition5, square_input_file$Farmers_in_Partition6, 
                                   square_input_file$Farmers_in_Partition7, square_input_file$Farmers_in_Partition8, 
                                   square_input_file$Farmers_in_Partition9, square_input_file$Farmers_in_Partition10)
  
  # Name columns
  colnames(population_partition_farmer_dat) = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                                    "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                                    "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                                    "Farmers_in_Partition10")
  
  # Remove non-informative data
  head_population_partition_farmer_dat = head(population_partition_farmer_dat, trimming_range)
  
  # Create output file name
  part_farmer_plot_out = paste("population_partition_farmer_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_farmer_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_farmer_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"))
  
  # Create Plot Legend
  legend("bottomright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"),
          legend = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                      "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                      "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                      "Farmers_in_Partition10"))
  
  dev.off()
  
  # ************************************************************************************************
  # HG Population Size Partition Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  population_partition_HG_dat = data.frame(square_input_file$HGs_in_Partition1, square_input_file$HGs_in_Partition2,
                                          square_input_file$HGs_in_Partition3, square_input_file$HGs_in_Partition4, 
                                          square_input_file$HGs_in_Partition5, square_input_file$HGs_in_Partition6, 
                                          square_input_file$HGs_in_Partition7, square_input_file$HGs_in_Partition8, 
                                          square_input_file$HGs_in_Partition9, square_input_file$HGs_in_Partition10)
  
  # Name columns
  colnames(population_partition_HG_dat) = c("HGs_in_Partition1", "HGs_in_Partition2", "HGs_in_Partition3", 
                                           "HGs_in_Partition4", "HGs_in_Partition5", "HGs_in_Partition6", 
                                           "HGs_in_Partition7", "HGs_in_Partition8", "HGs_in_Partition9", 
                                           "HGs_in_Partition10")
  
  # Remove non-informative data
  head_population_partition_HG_dat = head(population_partition_HG_dat, trimming_range)
  
  # Create output file name
  part_HG_plot_out = paste("population_partition_HG_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_HG_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_HG_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"))
  
  # Create Plot Legend
  legend("topright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"),
         legend = c("HGs_in_Partition1", "HGs_in_Partition2", "HGs_in_Partition3", 
                      "HGs_in_Partition4", "HGs_in_Partition5", "HGs_in_Partition6", 
                      "HGs_in_Partition7", "HGs_in_Partition8", "HGs_in_Partition9", 
                      "HGs_in_Partition10"))
  
  dev.off()
  
  # ************************************************************************************************
  # Combined Population Size Partition Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  population_partition_combined_dat = data.frame(square_input_file$Farmers_in_Partition1, square_input_file$Farmers_in_Partition2,
                                            square_input_file$Farmers_in_Partition3, square_input_file$Farmers_in_Partition4, 
                                            square_input_file$Farmers_in_Partition5, square_input_file$Farmers_in_Partition6, 
                                            square_input_file$Farmers_in_Partition7, square_input_file$Farmers_in_Partition8, 
                                            square_input_file$Farmers_in_Partition9, square_input_file$Farmers_in_Partition10,
                                            square_input_file$HGs_in_Partition1, square_input_file$HGs_in_Partition2,
                                            square_input_file$HGs_in_Partition3, square_input_file$HGs_in_Partition4, 
                                            square_input_file$HGs_in_Partition5, square_input_file$HGs_in_Partition6, 
                                            square_input_file$HGs_in_Partition7, square_input_file$HGs_in_Partition8, 
                                            square_input_file$HGs_in_Partition9, square_input_file$HGs_in_Partition10)
  
  # Name columns
  colnames(population_partition_combined_dat) = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                                             "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                                             "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                                             "Farmers_in_Partition10", "HGs_in_Partition1", "HGs_in_Partition2", 
                                             "HGs_in_Partition3", "HGs_in_Partition4", "HGs_in_Partition5", 
                                             "HGs_in_Partition6", "HGs_in_Partition7", "HGs_in_Partition8", 
                                             "HGs_in_Partition9", "HGs_in_Partition10")
  
  # Remove non-informative data
  head_population_partition_combined_dat = head(population_partition_combined_dat, trimming_range)
  
  # Create output file name
  part_combined_plot_out = paste("population_partition_combined_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_combined_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_combined_dat, xlab = "Year", ylab = "Population Size", pch = 20, 
          lwd = c(7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"))
  
  # Create Plot Legend
  legend("bottomright", pch = 20, y.intersp=1, cex = 0.65,
         lwd = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"),
         legend = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                    "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                    "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                    "Farmers_in_Partition10", "HGs_in_Partition1", "HGs_in_Partition2", 
                    "HGs_in_Partition3", "HGs_in_Partition4", "HGs_in_Partition5", 
                    "HGs_in_Partition6", "HGs_in_Partition7", "HGs_in_Partition8", 
                    "HGs_in_Partition9", "HGs_in_Partition10"))
  
  dev.off()
  
}

if (file.exists(general_input_name)){
  
  # ************************************************************************************************
  # Population Size Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  population_dat = data.frame(general_input_file$PopulationSize, 
                              general_input_file$TotalFarmers, general_input_file$TotalHGs)
  
  # Name columns
  colnames(population_dat) = c("PopulationSize", "TotalFarmers", "TotalHGs")
  
  # Remove non-informative data
  head_population_dat = head(population_dat, trimming_range)
  
  # Create output file name
  pop_plot_out = paste("population_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(pop_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_dat, col = c("black", "blue", "red"), xlab = "Year", ylab = "Population Size", pch = 19)
  
  # Create Plot Legend
  legend("right", legend = c("Total Population", "Total Farmers", "Total Hunter Gatherers"), pch = 19, col = c("black", "blue", "red"))
  
  dev.off()
  
  # ************************************************************************************************
  # Population Ratio Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  RatioFarmertoHG = data.frame(general_input_file$RatioFarmertoHG)
  
  # Name columns
  colnames(RatioFarmertoHG) = c("RatioFarmertoHG")
  
  # Remove non-informative data
  head_RatioFarmertoHG = head(RatioFarmertoHG, trimming_range)
  
  # Create output file name
  RatioFarmertoHG_plot_out = paste("RatioFarmertoHG_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(RatioFarmertoHG_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_RatioFarmertoHG, col = c("blue"), xlab = "Year", ylab = "Ratio Farmers:HGs", pch = 19)
  
  dev.off()
  
  # ************************************************************************************************
  # Birth Rate Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  birthrate_dat = data.frame(general_input_file$NewBirths)
  
  # Name columns
  colnames(birthrate_dat) = "New Births"
  
  # Remove non-informative data
  head_birthrate_dat = head(birthrate_dat, trimming_range)
  
  # Create output file name
  birthrate_plot_out = paste("birth_rate_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(birthrate_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_birthrate_dat, col = "blue", xlab = "Year", ylab = "Number of New Births", pch = 20)
  
  dev.off()
  
  # ************************************************************************************************
  # Plot of Reproductive Aged Individuals
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  repro_dat = data.frame(general_input_file$PopulationSize, 
                              general_input_file$Num_Repro_Age_Inds)
  
  # Name columns
  colnames(repro_dat) = c("TotalPopulation", "Num_Repro_Age_Inds")
  
  # Remove non-informative data
  head_repro_dat = head(repro_dat, trimming_range)
  
  # Create output file name
  repro_plot_out = paste("repro_age_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(repro_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_repro_dat, col = c("black", "blue", "red"), xlab = "Year", ylab = "Population Size", pch = 19)
  
  # Create Plot Legend
  legend("right", legend = c("TotalPopulation", "Reproductive Age Individuals"), pch = 19, col = c("black", "blue"))
  
  dev.off()

} 

if (file.exists(general_input_with_deaths_name)){
  
  
  # ************************************************************************************************
  # Birth and Death Rate Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  birthdeathrate_dat = data.frame(general_input_file$NewBirths, general_input_file$NewDeaths)
  
  # Name columns
  colnames(birthdeathrate_dat) = c("NewBirths", "NewDeaths")
  
  # Remove non-informative data
  head_birthdeathrate_dat = head(birthdeathrate_dat, trimming_range)
  
  # Create output file name
  birthrate_plot_out = paste("birth_death_rate_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(birthrate_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_birthdeathrate_dat, col = c("blue", "red"), xlab = "Year", ylab = "Number of New Births/Deaths", pch = 20, lwd = 0.07)
          
  # Create Plot Legend
  legend("topright", legend = c("Births", "Deaths"), pch = 20, lwd = 0.07, col = c("blue", "red"))
          
  dev.off()
  
  # ************************************************************************************************
  # Death Rate Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  deathrate_dat = data.frame(general_input_file$NewDeaths)
  
  # Name columns
  colnames(deathrate_dat) = "New Deaths"
  
  # Remove non-informative data
  head_deathrate_dat = head(deathrate_dat, trimming_range)
  
  # Create output file name
  deathrate_plot_out = paste("death_rate_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(deathrate_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_deathrate_dat, col = "blue", xlab = "Year", ylab = "Number of New Deaths", pch = 20)
  
  dev.off()
  
}

