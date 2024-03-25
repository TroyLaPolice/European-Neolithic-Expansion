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

simulation_data_filtered = simulation_data[Learning_Prob_n == 0.0 & TotalHGs == 0]

setwd("/path/to/empirical/data")

ancient_data = data.table(read.table("qpAdm_Fixed__EuroHG_AnatoliaN.txt", header = T))

headers = data.table(read.table("qpadm_headers_square.csv", header = T, sep=","))

ancient_data_filtered = ancient_data[GroupMeanLat < 55 & GroupMeanBP > 4500 & GroupMeanBP < 8000 & GroupNinds > 2]

# Remove greenland individual

ancient_data_filtered = ancient_data_filtered[GroupMeanLong>-20]

# Comp dist from Ankara (29, 41)

ancient_data_filtered[, distFromAnkara := distm(c(GroupMeanLong, GroupMeanLat), rev(c(39.93, 32.85)), fun = distHaversine), target]

ancient_data_filtered[, distFromAnkara_km := distFromAnkara*0.001]

ancient_data_filtered[, time_slice := cut(GroupMeanBP/1000, breaks = seq(12,1,-0.5))]

# Read in Param File
#param_file = data.table(read.table("param_inputs_clean_filtered_fitting.txt", header = T, sep=","))
param_file = data.table(read.table("param_inputs_clean_filtered_assortative_mating.csv", header = T, sep=","))

# ----------------------------------------------------------------------------------------------------------------
# 
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
for (param in param_file$Assortative_Mating_n)
{
  #individual_param_run = simulation_data_filtered[Learning_Prob_n == param]
  individual_param_run = simulation_data_filtered[Assortative_Mating_n == param]

  sim_farming_ancestry = (individual_param_run$Farmer_Ancestry_Bin_All)
  sim_x_dist = individual_param_run$Mid_Point_km
  ancient_x_dist = ancient_data_filtered$distFromAnkara_km
  
  interpolated = interpolate(ancient_x_dist, sim_x_dist, sim_farming_ancestry)
  
  ancient_data_filtered[, approximated_ancestry := as.numeric(0)]
  ancient_data_filtered[, approximated_ancestry := interpolated]
  
  ancient_data_filtered[, distance_from_sim := as.numeric(0)]
  ancient_data_filtered[, distance_from_sim := (scaled_weight_Anatolia_N - approximated_ancestry)]
  
  ancient_data_filtered[, distance_from_sim_squared := as.numeric(0)]
  ancient_data_filtered[, distance_from_sim_squared := (distance_from_sim^2)]
  
  ancient_data_filtered[, mean_distance_from_sim_squared := as.numeric(0)]
  ancient_data_filtered[, mean_distance_from_sim_squared := mean(distance_from_sim_squared)]
  
  ancient_data_filtered[, log_likelihood := as.numeric(0)]
  ancient_data_filtered[, log_likelihood := (dnorm(scaled_weight_Anatolia_N, approximated_ancestry, se.p1, log = T))]
  
  ancient_data_filtered[, sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, sum_log_likelihood := sum(log_likelihood)]
  
  ancient_data_filtered[, exponential_sum_log_likelihood := as.numeric(0)]
  ancient_data_filtered[, exponential_sum_log_likelihood := (exp(sum_log_likelihood))]
  
  ancient_data_filtered[, Learning_Prob_n := as.numeric(0)]
  ancient_data_filtered[, Learning_Prob_n := param]
  
  ancient_data_filtered[, Assortative_Mating_n := as.numeric(0)]
  ancient_data_filtered[, Assortative_Mating_n := param]
  
  main_ancestry_dt = rbind(main_ancestry_dt, ancient_data_filtered)
  
}

#ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))
ggplot(main_ancestry_dt) + geom_point(aes(Assortative_Mating_n, sum_log_likelihood))

main_ancestry_dt_best_fit_quad = main_ancestry_dt[Assortative_Mating_n >= 0.9]

#quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)
quadraticModel = lm(sum_log_likelihood ~ Assortative_Mating_n + I(Assortative_Mating_n^2), data=main_ancestry_dt_best_fit_quad)

#l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])
l_doubleprime = -(quadraticModel[["coefficients"]][["I(Assortative_Mating_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

#LearningValues = seq(0, 0.006, 0.0001)
AMValues = seq(0.9, 1, 0.001)

#sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))
sum_log_likelihoodPredict = predict(quadraticModel,list(Assortative_Mating_n=AMValues, Assortative_Mating_n2=AMValues^2))

# This just gets an approximate best fitting value- which of your tested params fits the best:
#           approximate_prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]

# To get the actual best fitting prediction put this in the console:
#           quadraticModel["coefficients"]

# Then the prediction will equal this equation. Get the Learning_Probs from the coefficients
# prediction = ((-Learning_Prob_n) / (2*[I(Learning_Prob_n^2)]))

#confidence_interval_floor = prediction-confidence_interval
#confidence_interval_ceiling = prediction+confidence_interval

setwd(wd)

#create scatterplot of original data values
pdf(file="sum_log_likelihoodPredict.pdf")
#plot(main_ancestry_dt_best_fit_quad$Learning_Prob_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16, xlab="Learning Rate", ylab="Sum Log Likelihood")
plot(main_ancestry_dt_best_fit_quad$Assortative_Mating_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16, xlab="Assortative Mating", ylab="Sum Log Likelihood", ylim =(c(-100000,-500)))
#add predicted lines based on quadratic regression model
lines(AMValues, sum_log_likelihoodPredict, col='blue')
dev.off()
