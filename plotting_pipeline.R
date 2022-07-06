# Load Libraries
library(ggplot2)
library(tidyr)

# ****************************************************
# Set WD, input and output file name additions, 
#   as well as trimming range for non-informative data
# ****************************************************

setwd("/home/tml5905/Documents/HuberLab/HunterGatherFarmerInteractions/longrun_test/7-5")

input_file_extention = "_min_mating_age"  # Should be the same as the ending of your SLiM output files
                                #   if endings are default this will be an empty string.

output_file_extention = "_min_mating_age" # Addition to attach to the output files
                                #   Best practice is to keep these the same for simplicity.

# Create file names
square_input_name = paste("sim_square_wave_stats_per_year", input_file_extention, ".txt", sep = "")
general_input_name = paste("sim_pop_stats_per_year", input_file_extention, ".txt", sep = "")

# Check if files exist and read in those that do
if ((file.exists(square_input_name)) && (file.exists(general_input_name))){
  square_input_file = read.csv(square_input_name, sep = " ", header = TRUE)
  general_input_file = read.csv(general_input_name, sep = " ", header = TRUE)
} else if (file.exists(square_input_name)){
  square_input_file = read.csv(square_input_name, sep = " ", header = TRUE)
} else if (file.exists(general_input_name)){
  general_input_file = read.csv(general_input_name, sep = " ", header = TRUE)
} else{
  print("Requested files do not exist")
}

# If the model was run on a square, plot the data
if (file.exists(square_input_name)){
  
  # ****************************************************
  # Farmer:HG Ratio Quadrant Plot
  # ****************************************************
  
  # Push ratio data to its own data frame
  ratio_dat = data.frame(square_input_file$Year,
                         square_input_file$RatioFarmerToHG, square_input_file$RatioFarmerToHG_Quad1,
                         square_input_file$RatioFarmerToHG_Quad2, square_input_file$RatioFarmerToHG_Quad3, 
                         square_input_file$RatioFarmerToHG_Quad4, square_input_file$RatioFarmerToHG_Quad5, 
                         square_input_file$RatioFarmerToHG_Quad6, square_input_file$RatioFarmerToHG_Quad7, 
                         square_input_file$RatioFarmerToHG_Quad8, square_input_file$RatioFarmerToHG_Quad9, 
                         square_input_file$RatioFarmerToHG_Quad10)
  
  
  # Name columns
  colnames(ratio_dat) = c("Year", "Total_RatioFarmerToHG", "RatioFarmerToHG_Quad1", "RatioFarmerToHG_Quad2", 
                           "RatioFarmerToHG_Quad3", "RatioFarmerToHG_Quad4", "RatioFarmerToHG_Quad5", 
                           "RatioFarmerToHG_Quad6", "RatioFarmerToHG_Quad7", "RatioFarmerToHG_Quad8", 
                           "RatioFarmerToHG_Quad9", "RatioFarmerToHG_QuadTen")
  
  # Remove non-informative ratio data
  head_ratio_dat = head(ratio_dat, -4500)
  
  # Convert to tidy format
  ratio_dat_tidyr = tidyr::pivot_longer(head_ratio_dat, cols=c("Total_RatioFarmerToHG", "RatioFarmerToHG_Quad1", "RatioFarmerToHG_Quad2", 
                                                               "RatioFarmerToHG_Quad3", "RatioFarmerToHG_Quad4", "RatioFarmerToHG_Quad5", 
                                                               "RatioFarmerToHG_Quad6", "RatioFarmerToHG_Quad7", "RatioFarmerToHG_Quad8", 
                                                               "RatioFarmerToHG_Quad9", "RatioFarmerToHG_QuadTen"),
                                        names_to = "variable")
  
  # Plot
  ratio_plot = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=variable)) +
    geom_bar(stat="identity", position="dodge") +  geom_col(width = 1)
  
  # Create output file name
  ratio_plot_out = paste("ratio_plot", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(ratio_plot_out, plot = ratio_plot, units = "in", width = 10, height = 5, device="tiff", dpi=1000)

    
  # ****************************************************
  # Farmer Population Size Quadrant Plot
  # ****************************************************
  
  # Push population size data to its own data frame
  population_quad_farmer_dat = data.frame(square_input_file$Farmers_in_Quadrant1, square_input_file$Farmers_in_Quadrant2,
                                   square_input_file$Farmers_in_Quadrant3, square_input_file$Farmers_in_Quadrant4, 
                                   square_input_file$Farmers_in_Quadrant5, square_input_file$Farmers_in_Quadrant6, 
                                   square_input_file$Farmers_in_Quadrant7, square_input_file$Farmers_in_Quadrant8, 
                                   square_input_file$Farmers_in_Quadrant9, square_input_file$Farmers_in_Quadrant10)
  
  # Name columns
  colnames(population_quad_farmer_dat) = c("Farmers_in_Quadrant1", "Farmers_in_Quadrant2", "Farmers_in_Quadrant3", 
                                    "Farmers_in_Quadrant4", "Farmers_in_Quadrant5", "Farmers_in_Quadrant6", 
                                    "Farmers_in_Quadrant7", "Farmers_in_Quadrant8", "Farmers_in_Quadrant9", 
                                    "Farmers_in_Quadrant10")
  
  # Remove non-informative data
  head_population_quad_farmer_dat = head(population_quad_farmer_dat, -4500)
  
  # Create output file name
  quad_farmer_plot_out = paste("population_quad_farmer_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(quad_farmer_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_quad_farmer_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"))
  
  # Create Plot Legend
  legend("bottomright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"),
          legend = c("Farmers_in_Quadrant1", "Farmers_in_Quadrant2", "Farmers_in_Quadrant3", 
                               "Farmers_in_Quadrant4", "Farmers_in_Quadrant5", "Farmers_in_Quadrant6", 
                               "Farmers_in_Quadrant7", "Farmers_in_Quadrant8", "Farmers_in_Quadrant9", 
                               "Farmers_in_Quadrant10"))
  
  dev.off()
  
  # ****************************************************
  # HG Population Size Quadrant Plot
  # ****************************************************
  
  # Push population size data to its own data frame
  population_quad_HG_dat = data.frame(square_input_file$HGs_in_Quadrant1, square_input_file$HGs_in_Quadrant2,
                                          square_input_file$HGs_in_Quadrant3, square_input_file$HGs_in_Quadrant4, 
                                          square_input_file$HGs_in_Quadrant5, square_input_file$HGs_in_Quadrant6, 
                                          square_input_file$HGs_in_Quadrant7, square_input_file$HGs_in_Quadrant8, 
                                          square_input_file$HGs_in_Quadrant9, square_input_file$HGs_in_Quadrant10)
  
  # Name columns
  colnames(population_quad_HG_dat) = c("HGs_in_Quadrant1", "HGs_in_Quadrant2", "HGs_in_Quadrant3", 
                                           "HGs_in_Quadrant4", "HGs_in_Quadrant5", "HGs_in_Quadrant6", 
                                           "HGs_in_Quadrant7", "HGs_in_Quadrant8", "HGs_in_Quadrant9", 
                                           "HGs_in_Quadrant10")
  
  # Remove non-informative data
  head_population_quad_HG_dat = head(population_quad_HG_dat, -4500)
  
  # Create output file name
  quad_HG_plot_out = paste("population_quad_HG_data_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(quad_HG_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_quad_HG_dat, xlab = "Year", ylab = "Population Size", pch = 19, 
          col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"))
  
  # Create Plot Legend
  legend("topright", pch = 19, y.intersp=1,
         col = c("#d00083", "#8e2090", "#5b3bf5", "#1e5187", "#74ffdf", "#b2ceee", "#b3ff80", "#70a71b", "#f8ef1c", "#ff6e2f"),
         legend = c("HGs_in_Quadrant1", "HGs_in_Quadrant2", "HGs_in_Quadrant3", 
                    "HGs_in_Quadrant4", "HGs_in_Quadrant5", "HGs_in_Quadrant6", 
                    "HGs_in_Quadrant7", "HGs_in_Quadrant8", "HGs_in_Quadrant9", 
                    "HGs_in_Quadrant10"))
  
  dev.off()
  
}

if (file.exists(general_input_name)){
  
  # ****************************************************
  # Population Size Plot
  # ****************************************************
  
  # Push population size data to its own data frame
  population_dat = data.frame(general_input_file$PopulationSize, 
                              general_input_file$TotalFarmers, general_input_file$TotalHGs)
  
  # Name columns
  colnames(population_dat) = c("PopulationSize", "TotalFarmers", "TotalHGs")
  
  # Remove non-informative data
  head_population_dat = head(population_dat, -4500)
  
  # Create output file name
  pop_plot_out = paste("population_plot", output_file_extention, ".tiff", sep = "")
  
  # Save pop plot
  tiff(pop_plot_out, units = "in", width = 10, height = 5, res = 1000)
  
  # Plot
  matplot(head_population_dat, col = c("black", "blue", "red"), xlab = "Year", ylab = "Population Size", pch = 19)
  
  # Create Plot Legend
  legend("topleft", legend = c("Total Population", "Total Farmers", "Total Hunter Gatherers"), pch = 19, col = c("black", "blue", "red"))
  
  dev.off()

} 
