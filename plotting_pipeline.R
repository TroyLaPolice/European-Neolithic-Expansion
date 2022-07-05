# Loab Libraries
library(ggplot2)
library(tidyr)

# Set WD and input and output file name additions

setwd("/home/tml5905/Documents/HuberLab/HunterGatherFarmerInteractions")

input_file_extention = "_test"  # Should be the same as the ending of your SLiM output files
                                #   if endings are default this will be an empty string.

output_file_extention = "_test" # Addition to attach to the output files
                                #   Best practice is to keep these the same for simplicity

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
  
  # Plot ratios
  
  # Push ratio data to its own data frame
  ratio_dat = data.frame(square_input_file$Year,
                         square_input_file$RatioFarmerToHG, square_input_file$RatioFarmerToHG_Quad1,
                         square_input_file$RatioFarmerToHG_Quad2, square_input_file$RatioFarmerToHG_Quad3, 
                         square_input_file$RatioFarmerToHG_Quad4, square_input_file$RatioFarmerToHG_Quad5, 
                         square_input_file$RatioFarmerToHG_Quad6, square_input_file$RatioFarmerToHG_Quad7, 
                         square_input_file$RatioFarmerToHG_Quad8, square_input_file$RatioFarmerToHG_Quad9, 
                         square_input_file$RatioFarmerToHG_Quad10)
  
  
  # Name columns
  colnames(ratio_dat) <- c("Year", "Total_RatioFarmerToHG", "RatioFarmerToHG_Quad1", "RatioFarmerToHG_Quad2", 
                           "RatioFarmerToHG_Quad3", "RatioFarmerToHG_Quad4", "RatioFarmerToHG_Quad5", 
                           "RatioFarmerToHG_Quad6", "RatioFarmerToHG_Quad7", "RatioFarmerToHG_Quad8", 
                           "RatioFarmerToHG_Quad9", "RatioFarmerToHG_QuadTen")
  
  # Remove non-informative ratio data
  head_ratio_dat = head(ratio_dat, -7200)
  
  # Convert to tidy format
  ratio_dat_tidyr = tidyr::pivot_longer(head_ratio_dat, cols=c("Total_RatioFarmerToHG", "RatioFarmerToHG_Quad1", "RatioFarmerToHG_Quad2", 
                                                               "RatioFarmerToHG_Quad3", "RatioFarmerToHG_Quad4", "RatioFarmerToHG_Quad5", 
                                                               "RatioFarmerToHG_Quad6", "RatioFarmerToHG_Quad7", "RatioFarmerToHG_Quad8", 
                                                               "RatioFarmerToHG_Quad9", "RatioFarmerToHG_QuadTen"),
                                        names_to='variable', values_to="value")
  
  # Plot
  ratio_plot = ggplot(ratio_dat_tidyr, aes(x=Year, y=value, fill=variable)) +
    geom_bar(stat='identity', position='dodge') +  geom_col(width = 1)
  
  # Create output filename
  plot_out = paste("ratio_plot", output_file_extention, ".tiff", sep = "")
  
  # Save ratio plot
  ggsave(plot_out, plot = last_plot(), units = "in", width = 10, height = 5, device='tiff', dpi=700)
  
}
