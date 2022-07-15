# Load Libraries
library(ggplot2)
library(tidyr)

# ********************************************************************************************************************************************
# Set WD, input and output file name additions, as well as trimming range for non-informative data
# ********************************************************************************************************************************************

setwd("/home/tml5905/Documents/HuberLab/HunterGatherFarmerInteractions/longrun_test/7-14")

input_file_extention = "_loop_test"  # Should be the same as the ending of your SLiM output files
                                #   if endings are default this will be an empty string.

output_file_extention = "_loop_test" # Addition to attach to the output files
                                #   Best practice is to keep these the same for simplicity.

trimming_range = -112

# Create file names
square_input_name = paste("sim_square_wave_stats_per_year", input_file_extention, ".txt", sep = "")
general_input_name = paste("sim_pop_stats_per_year", input_file_extention, ".txt", sep = "")
general_input_with_deaths_name = paste("sim_pop_stats_per_year_with_death", input_file_extention, ".csv", sep = "")

# ********************************************************************************************************************************************
# Read in Files
# ********************************************************************************************************************************************

# Check if files exist and read in those that do
if ((file.exists(square_input_name)) && (file.exists(general_input_with_deaths_name))){
  square_input_file = read.csv(square_input_name, sep = ",", header = TRUE)
  general_input_file = read.csv(general_input_with_deaths_name, sep = ",", header = TRUE)
} else if ((file.exists(square_input_name)) && (file.exists(general_input_name))){
  square_input_file = read.csv(square_input_name, sep = ",", header = TRUE)
  general_input_file = read.csv(general_input_name, sep = ",", header = TRUE)
} else if (file.exists(square_input_name)){
  square_input_file = read.csv(square_input_name, sep = ",", header = TRUE)
} else if (file.exists(general_input_with_deaths_name)){
  general_input_file = read.csv(general_input_with_deaths_name, sep = ",", header = TRUE)
} else if (file.exists(general_input_name)){
  general_input_file = read.csv(general_input_name, sep = ",", header = TRUE)
} else{
  print("Requested files do not exist")
}

# ********************************************************************************************************************************************
# Plotting
# ********************************************************************************************************************************************

# If the model was run on a square, plot the data
if (file.exists(square_input_name)){
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot
  # ************************************************************************************************
  
  # Push ratio data to its own data frame
  ratio_dat = data.frame(square_input_file$Year,
                         square_input_file$RatioFarmerToHG_Partition1,
                         square_input_file$RatioFarmerToHG_Partition2, square_input_file$RatioFarmerToHG_Partition3, 
                         square_input_file$RatioFarmerToHG_Partition4, square_input_file$RatioFarmerToHG_Partition5, 
                         square_input_file$RatioFarmerToHG_Partition6, square_input_file$RatioFarmerToHG_Partition7, 
                         square_input_file$RatioFarmerToHG_Partition8, square_input_file$RatioFarmerToHG_Partition9, 
                         square_input_file$RatioFarmerToHG_Partition10, square_input_file$RatioFarmerToHG_Partition11,
                         square_input_file$RatioFarmerToHG_Partition12, square_input_file$RatioFarmerToHG_Partition13, 
                         square_input_file$RatioFarmerToHG_Partition14, square_input_file$RatioFarmerToHG_Partition15, 
                         square_input_file$RatioFarmerToHG_Partition16, square_input_file$RatioFarmerToHG_Partition17, 
                         square_input_file$RatioFarmerToHG_Partition18, square_input_file$RatioFarmerToHG_Partition19, 
                         square_input_file$RatioFarmerToHG_Partition20)
  
  
  # Name columns
  colnames(ratio_dat) = c("Year", "1", "2", 
                           "3", "4", "5", 
                           "6", "7", "8", 
                           "9", "10", "11", "12", 
                          "13", "14", "15", 
                          "16", "17", "18", 
                          "19", "20")
  
  
  # Remove non-informative ratio data
  head_ratio_dat = head(ratio_dat, trimming_range)
  
  # Convert to tidy format
  ratio_dat_tidyr = tidyr::pivot_longer(head_ratio_dat, cols=c("1", "2", 
                                                               "3", "4", "5", 
                                                               "6", "7", "8", 
                                                               "9", "10", "11", "12", 
                                                               "13", "14", "15", 
                                                               "16", "17", "18", 
                                                               "19", "20"),
                                        names_to = "variable")
  
  Partition = as.numeric(ratio_dat_tidyr$variable)
  
  # Plot
  ratio_plot = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=Partition)) + geom_col(width = 1) + labs(y="Partition")
  
  # Create output file name
  ratio_plot_out = paste("ratio_partition_plot", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out, plot = ratio_plot, units = "in", width = 10, height = 5, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot2
  # ************************************************************************************************
  
  ratio_plot2 = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=Partition)) + geom_col(width = 1) + facet_grid(as.numeric(variable)~.)
  
  # Create output file name
  ratio_plot_out2 = paste("ratio_partition_plot2", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out2, plot = ratio_plot2, units = "in", width = 10, height = 12, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer:HG Ratio Partition Plot3
  # ************************************************************************************************
  
  ratio_plot3 = ggplot(ratio_dat_tidyr[ratio_dat_tidyr$Year%%50 == 0 & ratio_dat_tidyr$variable != "Total",], aes(x=as.numeric(variable), y=value, col=factor(Year), group=Year)) + geom_line() 
  
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
                                   square_input_file$Farmers_in_Partition9, square_input_file$Farmers_in_Partition10, 
                                   square_input_file$Farmers_in_Partition11, square_input_file$Farmers_in_Partition12,
                                   square_input_file$Farmers_in_Partition13, square_input_file$Farmers_in_Partition14, 
                                   square_input_file$Farmers_in_Partition15, square_input_file$Farmers_in_Partition16, 
                                   square_input_file$Farmers_in_Partition17, square_input_file$Farmers_in_Partition18, 
                                   square_input_file$Farmers_in_Partition19, square_input_file$Farmers_in_Partition20)
  
  # Name columns
  colnames(population_partition_farmer_dat) = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                                    "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                                    "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                                    "Farmers_in_Partition10", "Farmers_in_Partition11", "Farmers_in_Partition12", "Farmers_in_Partition13", 
                                    "Farmers_in_Partition14", "Farmers_in_Partition15", "Farmers_in_Partition16", 
                                    "Farmers_in_Partition17", "Farmers_in_Partition18", "Farmers_in_Partition19", 
                                    "Farmers_in_Partition20")
  
  # Remove non-informative data
  head_population_partition_farmer_dat = head(population_partition_farmer_dat, trimming_range)
  
  # Create output file name
  part_farmer_plot_out = paste("population_partition_farmer_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_farmer_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_farmer_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                  "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"))
  
  # Create Plot Legend
  legend("bottomright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                 "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"),
          legend = c("Partition1", "Partition2", "Partition3", 
                      "Partition4", "Partition5", "Partition6", 
                      "Partition7", "Partition8", "Partition9", 
                      "Partition10", "Partition11", "Partition12", "Partition13", 
                      "Partition14", "Partition15", "Partition16", 
                      "Partition17", "Partition18", "Partition19", 
                      "Partition20"), cex = 0.65)
  
  dev.off()
  
  
  
  # ************************************************************************************************
  # HG Population Size Partition Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  population_partition_HG_dat = data.frame(square_input_file$HGs_in_Partition1, square_input_file$HGs_in_Partition2,
                                               square_input_file$HGs_in_Partition3, square_input_file$HGs_in_Partition4, 
                                               square_input_file$HGs_in_Partition5, square_input_file$HGs_in_Partition6, 
                                               square_input_file$HGs_in_Partition7, square_input_file$HGs_in_Partition8, 
                                               square_input_file$HGs_in_Partition9, square_input_file$HGs_in_Partition10, 
                                               square_input_file$HGs_in_Partition11, square_input_file$HGs_in_Partition12,
                                               square_input_file$HGs_in_Partition13, square_input_file$HGs_in_Partition14, 
                                               square_input_file$HGs_in_Partition15, square_input_file$HGs_in_Partition16, 
                                               square_input_file$HGs_in_Partition17, square_input_file$HGs_in_Partition18, 
                                               square_input_file$HGs_in_Partition19, square_input_file$HGs_in_Partition20)
  
  # Name columns
  colnames(population_partition_HG_dat) = c("HGs_in_Partition1", "HGs_in_Partition2", "HGs_in_Partition3", 
                                                "HGs_in_Partition4", "HGs_in_Partition5", "HGs_in_Partition6", 
                                                "HGs_in_Partition7", "HGs_in_Partition8", "HGs_in_Partition9", 
                                                "HGs_in_Partition10", "HGs_in_Partition11", "HGs_in_Partition12", "HGs_in_Partition13", 
                                                "HGs_in_Partition14", "HGs_in_Partition15", "HGs_in_Partition16", 
                                                "HGs_in_Partition17", "HGs_in_Partition18", "HGs_in_Partition19", 
                                                "HGs_in_Partition20")
  
  # Remove non-informative data
  head_population_partition_HG_dat = head(population_partition_HG_dat, trimming_range)
  
  # Create output file name
  part_HG_plot_out = paste("population_partition_HG_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_HG_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_HG_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                  "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"))
  
  # Create Plot Legend
  legend("topright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                 "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"),
         legend = c("Partition1", "Partition2", "Partition3", 
                    "Partition4", "Partition5", "Partition6", 
                    "Partition7", "Partition8", "Partition9", 
                    "Partition10", "Partition11", "Partition12", "Partition13", 
                    "Partition14", "Partition15", "Partition16", 
                    "Partition17", "Partition18", "Partition19", 
                    "Partition20"), cex = 0.65)
  
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
                                                 square_input_file$Farmers_in_Partition11, square_input_file$Farmers_in_Partition12,
                                                 square_input_file$Farmers_in_Partition13, square_input_file$Farmers_in_Partition14, 
                                                 square_input_file$Farmers_in_Partition15, square_input_file$Farmers_in_Partition16, 
                                                 square_input_file$Farmers_in_Partition17, square_input_file$Farmers_in_Partition18, 
                                                 square_input_file$Farmers_in_Partition19, square_input_file$Farmers_in_Partition20,
                                                 square_input_file$HGs_in_Partition1, square_input_file$HGs_in_Partition2,
                                                 square_input_file$HGs_in_Partition3, square_input_file$HGs_in_Partition4, 
                                                 square_input_file$HGs_in_Partition5, square_input_file$HGs_in_Partition6, 
                                                 square_input_file$HGs_in_Partition7, square_input_file$HGs_in_Partition8, 
                                                 square_input_file$HGs_in_Partition9, square_input_file$HGs_in_Partition10, 
                                                 square_input_file$HGs_in_Partition11, square_input_file$HGs_in_Partition12,
                                                 square_input_file$HGs_in_Partition13, square_input_file$HGs_in_Partition14, 
                                                 square_input_file$HGs_in_Partition15, square_input_file$HGs_in_Partition16, 
                                                 square_input_file$HGs_in_Partition17, square_input_file$HGs_in_Partition18, 
                                                 square_input_file$HGs_in_Partition19, square_input_file$HGs_in_Partition20)
  
  # Name columns
  colnames(population_partition_combined_dat) = c("Farmers_in_Partition1", "Farmers_in_Partition2", "Farmers_in_Partition3", 
                                             "Farmers_in_Partition4", "Farmers_in_Partition5", "Farmers_in_Partition6", 
                                             "Farmers_in_Partition7", "Farmers_in_Partition8", "Farmers_in_Partition9", 
                                             "Farmers_in_Partition10", "Farmers_in_Partition11", "Farmers_in_Partition12", "Farmers_in_Partition13", 
                                             "Farmers_in_Partition14", "Farmers_in_Partition15", "Farmers_in_Partition16", 
                                             "Farmers_in_Partition17", "Farmers_in_Partition18", "Farmers_in_Partition19", 
                                             "Farmers_in_Partition11", 
                                             "HGs_in_Partition1", "HGs_in_Partition2",
                                             "HGs_in_Partition3", "HGs_in_Partition4", "HGs_in_Partition5", 
                                             "HGs_in_Partition6", "HGs_in_Partition7", "HGs_in_Partition8", 
                                             "HGs_in_Partition9", "HGs_in_Partition10",
                                             "HGs_in_Partition11", "HGs_in_Partition12",
                                             "HGs_in_Partition13", "HGs_in_Partition14", "HGs_in_Partition15", 
                                             "HGs_in_Partition16", "HGs_in_Partition17", "HGs_in_Partition18", 
                                             "HGs_in_Partition19", "HGs_in_Partition20")
  
  # Remove non-informative data
  head_population_partition_combined_dat = head(population_partition_combined_dat, trimming_range)
  
  # Create output file name
  part_combined_plot_out = paste("population_partition_combined_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(part_combined_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_partition_combined_dat, xlab = "Year", ylab = "Population Size", pch = 20, 
          lwd = c(7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                  "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"))
  
  # Create Plot Legend
  legend("right", pch = 20, y.intersp=1, cex = 0.38,
         lwd = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f",
                 "#a19ae3", "#97cc9d", "#feb67a", "#ffb3f8", "#ffff79", "#07bca1", "#4981ff", "#1b8c35", "#febb38", "#a436f9"),
         legend = c("Farmers_Partition1", "Farmers_Partition2", "Farmers_Partition3", 
                    "Farmers_Partition4", "Farmers_Partition5", "Farmers_Partition6", 
                    "Farmers_Partition7", "Farmers_Partition8", "Farmers_Partition9", 
                    "Farmers_Partition10", "Farmers_Partition11", "Farmers_Partition12", 
                    "Farmers_Partition13", "Farmers_Partition14", "Farmers_Partition15", 
                    "Farmers_Partition16", "Farmers_Partition17", "Farmers_Partition18", 
                    "Farmers_Partition19", "Farmers_Partition20", 
                    "HGs_Partition1", "HGs_Partition2", "HGs_Partition3", "HGs_Partition4", 
                    "HGs_Partition5", "HGs_Partition6", "HGs_Partition7", "HGs_Partition8", 
                    "HGs_Partition9", "HGs_Partition10", "HGs_Partition11", "HGs_Partition12",
                    "HGs_Partition13", "HGs_Partition14", "HGs_Partition15", 
                    "HGs_Partition16", "HGs_Partition17", "HGs_Partition18", 
                    "HGs_Partition19", "HGs_Partition20"))
  
  dev.off()
  
  # ************************************************************************************************
  # Farmer Ancestry Partition Plot - Farmers
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  partition_ancestry_dat_farmer = data.frame(square_input_file$Year, square_input_file$Farmer_Ancestry_Partition1_Farmers, square_input_file$Farmer_Ancestry_Partition2_Farmers,
                                                          square_input_file$Farmer_Ancestry_Partition3_Farmers, square_input_file$Farmer_Ancestry_Partition4_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition5_Farmers, square_input_file$Farmer_Ancestry_Partition6_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition7_Farmers, square_input_file$Farmer_Ancestry_Partition8_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition9_Farmers, square_input_file$Farmer_Ancestry_Partition10_Farmers,
                                                          square_input_file$Farmer_Ancestry_Partition11_Farmers, square_input_file$Farmer_Ancestry_Partition12_Farmers,
                                                          square_input_file$Farmer_Ancestry_Partition13_Farmers, square_input_file$Farmer_Ancestry_Partition14_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition15_Farmers, square_input_file$Farmer_Ancestry_Partition16_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition17_Farmers, square_input_file$Farmer_Ancestry_Partition18_Farmers, 
                                                          square_input_file$Farmer_Ancestry_Partition19_Farmers, square_input_file$Farmer_Ancestry_Partition20_Farmers)
  
  # Name columns
  colnames(partition_ancestry_dat_farmer) = c("Year", "1", "2", 
                                              "3", "4", "5", 
                                              "6", "7", "8", 
                                              "9", "10", "11", "12", 
                                              "13", "14", "15", 
                                              "16", "17", "18", 
                                              "19", "20")
  
  # Remove non-informative ratio data
  head_partition_ancestry_dat_farmer = head(partition_ancestry_dat_farmer, trimming_range)
  
  # Convert to tidy format
  ancestry_dat_farmer_tidyr = tidyr::pivot_longer(head_partition_ancestry_dat_farmer, cols=c("1", "2", 
                                                                                             "3", "4", "5", 
                                                                                             "6", "7", "8", 
                                                                                             "9", "10", "11", "12", 
                                                                                             "13", "14", "15", 
                                                                                             "16", "17", "18", 
                                                                                             "19", "20"),
                                        names_to = "variable")
  
  Partition = as.numeric(ancestry_dat_farmer_tidyr$variable)
  
  # Plot
  ancestry_partition_farmer_plot = ggplot(ancestry_dat_farmer_tidyr, aes(x=Year, y=value, fill=Partition)) + geom_col(width = 1) +
                                      facet_grid(as.numeric(variable)~.) + labs(y="Farming Ancestry")
  
  # Create output file name
  ancestry_partition_plot_farmer_out = paste("ancestry_partition_plot_farmers", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_plot_farmer_out, plot = ancestry_partition_farmer_plot, units = "in", width = 10, height = 12, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer Ancestry Partition Plot - HGs
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  partition_ancestry_dat_HG = data.frame(square_input_file$Year, square_input_file$Farmer_Ancestry_Partition1_HGs, square_input_file$Farmer_Ancestry_Partition2HG,
                                         square_input_file$Farmer_Ancestry_Partition3HG, square_input_file$Farmer_Ancestry_Partition4HG, 
                                         square_input_file$Farmer_Ancestry_Partition5HG, square_input_file$Farmer_Ancestry_Partition6HG, 
                                         square_input_file$Farmer_Ancestry_Partition7HG, square_input_file$Farmer_Ancestry_Partition8HG, 
                                         square_input_file$Farmer_Ancestry_Partition9HG, square_input_file$Farmer_Ancestry_Partition10HG,
                                         square_input_file$Farmer_Ancestry_Partition11HG, square_input_file$Farmer_Ancestry_Partition12HG,
                                         square_input_file$Farmer_Ancestry_Partition13HG, square_input_file$Farmer_Ancestry_Partition14HG, 
                                         square_input_file$Farmer_Ancestry_Partition15HG, square_input_file$Farmer_Ancestry_Partition16HG, 
                                         square_input_file$Farmer_Ancestry_Partition17HG, square_input_file$Farmer_Ancestry_Partition18HG, 
                                         square_input_file$Farmer_Ancestry_Partition19HG, square_input_file$Farmer_Ancestry_Partition20HG)
  
  # Name columns
  colnames(partition_ancestry_dat_HG) = c("Year", "1", "2", 
                                            "3", "4", "5", 
                                            "6", "7", "8", 
                                            "9", "10", "11", "12", 
                                            "13", "14", "15", 
                                            "16", "17", "18", 
                                            "19", "20")
  
  # Remove non-informative ratio data
  head_partition_ancestry_dat_HG = head(partition_ancestry_dat_HG, trimming_range)
  
  # Convert to tidy format
  ancestry_dat_HG_tidyr = tidyr::pivot_longer(head_partition_ancestry_dat_HG, cols=c("1", "2", 
                                                                                     "3", "4", "5", 
                                                                                     "6", "7", "8", 
                                                                                     "9", "10", "11", "12", 
                                                                                     "13", "14", "15", 
                                                                                     "16", "17", "18", 
                                                                                     "19", "20"),
                                                  names_to = "variable")
  
  Partition = as.numeric(ancestry_dat_HG_tidyr$variable)
  
  # Plot
  ancestry_partition_HG_plot = ggplot(ancestry_dat_HG_tidyr, aes(x=Year, y=value, fill=Partition)) + geom_col(width = 1) +
    facet_grid(as.numeric(variable)~.) + labs(y="Farming Ancestry")
  
  # Create output file name
  ancestry_partition_plot_HG_out = paste("ancestry_partition_plot_HGs", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_plot_HG_out, plot = ancestry_partition_HG_plot, units = "in", width = 10, height = 12, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer Ancestry Partition Plot - All
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  partition_ancestry_dat_All = data.frame(square_input_file$Year, square_input_file$Farmer_Ancestry_Partition1_All, square_input_file$Farmer_Ancestry_Partition2All,
                                             square_input_file$Farmer_Ancestry_Partition3All, square_input_file$Farmer_Ancestry_Partition4All, 
                                             square_input_file$Farmer_Ancestry_Partition5All, square_input_file$Farmer_Ancestry_Partition6All, 
                                             square_input_file$Farmer_Ancestry_Partition7All, square_input_file$Farmer_Ancestry_Partition8All, 
                                             square_input_file$Farmer_Ancestry_Partition9All, square_input_file$Farmer_Ancestry_Partition10All,
                                             square_input_file$Farmer_Ancestry_Partition11All, square_input_file$Farmer_Ancestry_Partition12All,
                                             square_input_file$Farmer_Ancestry_Partition13All, square_input_file$Farmer_Ancestry_Partition14All, 
                                             square_input_file$Farmer_Ancestry_Partition15All, square_input_file$Farmer_Ancestry_Partition16All, 
                                             square_input_file$Farmer_Ancestry_Partition17All, square_input_file$Farmer_Ancestry_Partition18All, 
                                             square_input_file$Farmer_Ancestry_Partition19All, square_input_file$Farmer_Ancestry_Partition20All)
  
  # Name columns
  colnames(partition_ancestry_dat_All) = c("Year", "1", "2", 
                                              "3", "4", "5", 
                                              "6", "7", "8", 
                                              "9", "10", "11", "12", 
                                              "13", "14", "15", 
                                              "16", "17", "18", 
                                              "19", "20")
  
  # Remove non-informative ratio data
  head_partition_ancestry_dat_All = head(partition_ancestry_dat_All, trimming_range)
  
  # Convert to tidy format
  ancestry_dat_All_tidyr = tidyr::pivot_longer(head_partition_ancestry_dat_All, cols=c("1", "2", 
                                                                                             "3", "4", "5", 
                                                                                             "6", "7", "8", 
                                                                                             "9", "10", "11", "12", 
                                                                                             "13", "14", "15", 
                                                                                             "16", "17", "18", 
                                                                                             "19", "20"),
                                                  names_to = "variable")
  
  Partition = as.numeric(ancestry_dat_All_tidyr$variable)
  
  # Plot
  ancestry_partition_All_plot = ggplot(ancestry_dat_All_tidyr, aes(x=Year, y=value, fill=Partition)) + geom_col(width = 1) +
    facet_grid(as.numeric(variable)~.) + labs(y="Farming Ancestry")
  
  # Create output file name
  ancestry_partition_plot_All_out = paste("ancestry_partition_plot_All", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_plot_All_out, plot = ancestry_partition_All_plot, units = "in", width = 10, height = 12, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer Ancestry Partition Plot2 - Farmers
  # ************************************************************************************************
  
  Partition = as.numeric(ancestry_dat_farmer_tidyr$variable)
  
  ancestry_partition_farmer_plot2 = ggplot(ancestry_dat_farmer_tidyr[ancestry_dat_farmer_tidyr$Year%%50 == 0 & ancestry_dat_farmer_tidyr$variable != "Total",], aes(x=as.numeric(variable), y=value, col=factor(Year), group=Year)) + geom_line() 
  
  # Create output file name
  ancestry_partition_farmer_plot2_out = paste("ancestry_partition_farmer_plot2", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_farmer_plot2_out, plot = ancestry_partition_farmer_plot2, units = "in", width = 10, height = 5, device="tiff", dpi=1000)

  # ************************************************************************************************
  # Farmer Ancestry Partition Plot2 - HGs
  # ************************************************************************************************
  
  ancestry_partition_HG_plot2 = ggplot(ancestry_dat_HG_tidyr[ancestry_dat_HG_tidyr$Year%%50 == 0 & ancestry_dat_HG_tidyr$variable != "Total",], aes(x=as.numeric(variable), y=value, col=factor(Year), group=Year)) + geom_line() 
  
  # Create output file name
  ancestry_partition_HG_plot2_out = paste("ancestry_partition_HG_plot2", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_HG_plot2_out, plot = ancestry_partition_HG_plot2, units = "in", width = 10, height = 5, device="tiff", dpi=1000)
  
  # ************************************************************************************************
  # Farmer Ancestry Partition Plot2 - All
  # ************************************************************************************************
  
  ancestry_partition_all_plot2 = ggplot(ancestry_dat_All_tidyr[ancestry_dat_All_tidyr$Year%%50 == 0 & ancestry_dat_All_tidyr$variable != "Total",], aes(x=as.numeric(variable), y=value, col=factor(Year), group=Year)) + geom_line() 
  
  # Create output file name
  ancestry_partition_all_plot2_out = paste("ancestry_partition_all_plot2", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ancestry_partition_all_plot2_out, plot = ancestry_partition_all_plot2, units = "in", width = 10, height = 5, device="tiff", dpi=1000)
  

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
  birthrate_dat = data.frame(general_input_file$New_Births)
  
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
  legend("bottomright", legend = c("TotalPopulation", "Reproductive Age Individuals"), pch = 19, col = c("black", "blue"))
  
  dev.off()
  
  # ************************************************************************************************
  # Birth Freq Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  birthfreq_dat = data.frame(general_input_file$ReproFreq)
  
  # Name columns
  colnames(birthfreq_dat) = "ReproductionFreq"
  
  # Remove non-informative data
  head_birthfreq_dat = head(birthfreq_dat, trimming_range)
  
  # Create output file name
  birthfreq_plot_out = paste("birth_freq_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(birthfreq_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_birthfreq_dat, col = "blue", xlab = "Year", ylab = "Rate", pch = 20)
  
  dev.off()
  
  # ************************************************************************************************
  # Plot of Farming Ancestry in Population
  # ************************************************************************************************
  
  # Push ancestry size data to its own data frame
  ancestry_dat = data.frame(general_input_file$Farmer_Ancestry_All,
                              general_input_file$Farmer_Ancestry_Farmers, 
                              general_input_file$Farmer_Ancestry_HGs)
  
  # Name columns
  colnames(ancestry_dat) = c("Farmer_Ancestry_All", "Farmer_Ancestry_All_Farmers", "Farmer_Ancestry_All_HGs")
  
  # Remove non-informative data
  head_ancestry_dat = head(ancestry_dat, trimming_range)
  
  # Create output file name
  ancestry_plot_out = paste("ancestry_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(ancestry_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_ancestry_dat, col = c("black", "blue", "red"), xlab = "Year", ylab = "Percent Farming Ancestry", pch = 20)
  
  # Create Plot Legend
  legend("right", legend = c("Farmer Ancestry in All Individuals", "Farmer Ancestry in Farmers", "Farmer Ancestry in HGs"), pch = 20, col = c("black", "blue", "red"))
  
  dev.off()
} 

if (file.exists(general_input_with_deaths_name)){
  
  
  # ************************************************************************************************
  # Birth and Death Rate Plot
  # ************************************************************************************************
  
  # Push population size data to its own data frame
  birthdeathrate_dat = data.frame(general_input_file$New_Births, general_input_file$New_Deaths)
  
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
  deathrate_dat = data.frame(general_input_file$New_Deaths)
  
  # Name columns
  colnames(deathrate_dat) = "NewDeaths"
  
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

