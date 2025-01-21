devtools::install_github("ropenscilabs/rnaturalearth")

devtools::install_github("ropenscilabs/rnaturalearthdata")
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")

library(rnaturalearth)
library(sp)
library(ggplot2)
library(data.table)
library(sf)
library(geosphere)
library(scatterpie)
library(ggmap)
library(Rmpfr)

# ----------------------------------------------------------------------------------------------------------------
# Set inputs and read in files to be used for analysis
# ----------------------------------------------------------------------------------------------------------------

map = ne_countries(scale = "medium", returnclass = "sf")

wd = ("/path/to/simulation/data")

setwd(wd)
simulation_data = data.table(read.table("all_sim_data.csv", header = T, sep=","))

simulation_data_filtered = simulation_data[TotalHGs == 0]

setwd("/path/to/ancient/dna/data")

headers = data.table(read.table("qpadm_headers_square_WGM.csv", header = T, sep=","))

# -----------------------------------------
# ANCIENT SWGMPLE FILTERING CODE
# -----------------------------------------

metaD_dir = "/path/to/ancient/dna/data" # directory containing the v62.0_1240k_public.anno.txt file
qpadm_dir = "/path/to/ancient/dna/data" # directory containing the qpAdm_Fixed_3sources_OldSteppe_WHGA_EEF__Summary.rds file

all_ancient_inds = data.table(read.table("LaPolice__SI_Table__qpAdm_estimates.txt", header = T, sep="\t"))
#filtered_inds = data.table(read.table("simulated_qpadm.txt", header = T, sep=","))

#' Keep only plausible model individuals with s.e < 0.022
filtered_inds = all_ancient_inds[se.p3 < 0.022 & p.value >= 0.01 &
                                   (weight.p1 >= 0 & weight.p1 <=1) &
                                   (weight.p2 >= 0 & weight.p2 <=1) &
                                   (weight.p3 >= 0 & weight.p3 <=1)]

# Remove individuals with > 5% Steppe ancestry
filtered_inds <- filtered_inds[weight.p1 <= 0.05]

# Steppe Ancestry Histogram (no ancestry > 5%)
hist(filtered_inds$weight.p1)
# WHG Ancestry Historgram
hist(filtered_inds$weight.p2)
# EEF Ancestry Historgram
hist(filtered_inds$weight.p3)


# Get the distance of each sample frm the origin (Ankara)
library(geosphere)
# Convert 'lat' and 'lon' to numeric
filtered_inds[, `:=`(
  lat = as.numeric(as.character(lat)),
  lon = as.numeric(as.character(lon))
)]
# Define Ankara's coordinates
Ankara_coords <- c(lon = 35, lat = 40)
# Calculate the distance to Ankara in kilometers
# Assuming 'd' is your data.table and has columns 'lat' and 'lon'
filtered_inds[, distFromAnkara_km := distHaversine(cbind(lon, lat), Ankara_coords) / 1000]


# This will give you a table of farming ancestry estimates and the 
# distance of each sample from the farming origin
filtered_inds[,.(distFromAnkara_km, weight.p3)] 

# Read in Param File
#param_file = data.table(read.table("param_inputs_clean_filtered_run_for_fitting_learning.txt", header = T, sep=","))
param_file = data.table(read.table("param_inputs_clean_filtered_WithinGroup_mating.csv", header = T, sep=","))

# ----------------------------------------------------------------------------------------------------------------
#  UNCOMMENT AND COMMENT THE APPLICABLE LINES IF THE FITTING IS FOR WITHIN-GROUP MATING OR LEARNING RATE
# ----------------------------------------------------------------------------------------------------------------

interpolate = function(ancient_x_dist, sim_x_dist, sim_farming_ancestry){
  
  interpolation = approx(sim_x_dist, sim_farming_ancestry, xout=ancient_x_dist, method="linear")
  
  # Pull out the distance from the approximation list
  approximated_ancestry = interpolation[[2]]
  
  View(interpolation)
  
  # Return value
  return(approximated_ancestry)
  
}

main_ancestry_dt = headers

#for (param in param_file$Learning_Prob_n)
for (param in param_file$WithinGroup_Mating_n)
{
  #individual_param_run = simulation_data_filtered[Learning_Prob_n == param]
  individual_param_run = simulation_data_filtered[WithinGroup_Mating_n == param]
  
  sim_farming_ancestry = (individual_param_run$Farmer_Ancestry_Partition_All)
  sim_x_dist = individual_param_run$Mid_Point_km
  ancient_x_dist = filtered_inds$distFromAnkara_km
  
  interpolated = interpolate(ancient_x_dist, sim_x_dist, sim_farming_ancestry)
  
  filtered_inds[, approximated_ancestry := as.numeric(0)]
  filtered_inds[, approximated_ancestry := interpolated]
  
  filtered_inds[, distance_from_sim := as.numeric(0)]
  filtered_inds[, distance_from_sim := (weight.p3 - approximated_ancestry)]
  
  filtered_inds[, distance_from_sim_squared := as.numeric(0)]
  filtered_inds[, distance_from_sim_squared := (distance_from_sim^2)]
  
  filtered_inds[, mean_distance_from_sim_squared := as.numeric(0)]
  filtered_inds[, mean_distance_from_sim_squared := mean(distance_from_sim_squared, na.rm=TRUE)]
  
  filtered_inds[, log_likelihood := as.numeric(0)]
  filtered_inds[, log_likelihood := (dnorm(weight.p3, approximated_ancestry, se.p3, log = T))]
  
  filtered_inds[, sum_log_likelihood := as.numeric(0)]
  filtered_inds[, sum_log_likelihood := sum(log_likelihood, na.rm=TRUE)]
  
  filtered_inds[, exponential_sum_log_likelihood := as.numeric(0)]
  filtered_inds[, exponential_sum_log_likelihood := (exp(sum_log_likelihood))]
  
  #filtered_inds[, Learning_Prob_n := as.numeric(0)]
  #filtered_inds[, Learning_Prob_n := param]
  
  filtered_inds[, WithinGroup_Mating_n := as.numeric(0)]
  filtered_inds[, WithinGroup_Mating_n := param]
  
  main_ancestry_dt = rbind(main_ancestry_dt, filtered_inds)
  
}

#ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))
ggplot(main_ancestry_dt) + geom_point(aes(WithinGroup_Mating_n, sum_log_likelihood))

ggplot(main_ancestry_dt) + geom_point(aes(distFromAnkara_km, weight.p3))

#main_ancestry_dt_best_fit_quad = main_ancestry_dt[Learning_Prob_n <= 0.0025]
main_ancestry_dt_best_fit_quad = main_ancestry_dt[WithinGroup_Mating_n >= 0.9]

#quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)
quadraticModel = lm(sum_log_likelihood ~ WithinGroup_Mating_n + I(WithinGroup_Mating_n^2), data=main_ancestry_dt_best_fit_quad)

#l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])
l_doubleprime = -(quadraticModel[["coefficients"]][["I(WithinGroup_Mating_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

#LearningValues = seq(0, 0.006, 0.0001)
WGMValues = seq(0.9, 1, 0.001)

#sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))
sum_log_likelihoodPredict = predict(quadraticModel,list(WithinGroup_Mating_n=WGMValues, WithinGroup_Mating_n2=WGMValues^2))

# This just gets an approximate best fitting value- which of your tested params fits the best:
           #approximate_prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]
           approximate_prediction = WGMValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]
# To get the actual best fitting prediction put this in the console:
#           quadraticModel["coefficients"]

# Or do this to pull them out as variables

coefficients = quadraticModel$coefficients
coefficients = as.vector(coefficients)
Learning_Prob_n_coefficient = coefficients[2]
ILearning_Prob_n2_coefficient = coefficients[3]

# Then the prediction will equal this equation. Get the Learning_Probs from the coefficients
# prediction = ((-Learning_Prob_n) / (2*[I(Learning_Prob_n^2)]))

# Or do it out with code, not by hand

prediction = ((-as.numeric(Learning_Prob_n_coefficient)) / (2*as.numeric(ILearning_Prob_n2_coefficient)))

confidence_interval_floor = prediction-confidence_interval
confidence_interval_ceiling = prediction+confidence_interval

prediction_rounded <- round(prediction, 6)
confidence_interval_floor_rounded <- round(confidence_interval_floor, 6)
confidence_interval_ceiling_rounded <- round(confidence_interval_ceiling, 6)

formatted_prediction_with_conf = paste(prediction_rounded, sprintf("[%s, %s]", confidence_interval_floor_rounded, confidence_interval_ceiling_rounded))

#wd_refit = paste0(wd, "/refit")

#setwd(wd)
setwd("/home/tml5905/Documents/figures_to_import_to_pptTEMP/latest_round_of_revisions/re_refits")

#if (!file.exists("refit")) {
#dir.create("refit")
#}

#setwd(wd_refit)

#create scatterplot of original data values
png(file="sum_log_likelihoodPredict_WithinGroup_Mating.png", width = 6, height = 6, units = "in", res = 900)
plot(main_ancestry_dt_best_fit_quad$WithinGroup_Mating_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16, xlab="WithinGroup Mating", ylab="Sum Log Likelihood")
#add predicted lines based on quadratic regression model
lines(WGMValues, sum_log_likelihoodPredict, col='blue')
dev.off()


# Fit regression to the ancient DNA data:

#Standard lm                            
empirical_slope_regression = lm(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)
#Robust regression package 1  
#empirical_slope_robust_regression = rlm(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)
#Robust regression package 2  
empirical_slope_robust_regression2 = lmrob(weight.p3 ~ I(distFromAnkara_km), data=filtered_inds)

formatted_prediction_with_conf
