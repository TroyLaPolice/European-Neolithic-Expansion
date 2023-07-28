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

# Walking distance function:
# To do

setwd("/home/tml5905/Documents/HunterGatherFarmerInteractions/empirical_data")

ancient_data = data.table(read.table("qpAdm_Fixed__EuroHG_AnatoliaN.txt", header = T))

headers = data.table(read.table("qpadm_headers.csv", header = T, sep=","))

ancient_data_filtered = ancient_data[GroupMeanLat < 55 & GroupMeanBP > 4500 & GroupMeanBP < 8000 & GroupNinds > 2]

simulation_data = data.table(read.table("all_sim_data.csv", header = T, sep=","))

simulation_data_filtered = simulation_data[Replicate_n == 4 & TotalHGs == 0]

# Remove greenland individual

ancient_data_filtered = ancient_data_filtered[GroupMeanLong>-20]

# Comp dist from Ankara (29, 41)

ancient_data_filtered[, distFromAnkara := distm(c(GroupMeanLong, GroupMeanLat), rev(c(39.93, 32.85)), fun = distHaversine), target]

ancient_data_filtered[, distFromAnkara_km := distFromAnkara*0.001]

ancient_data_filtered[, time_slice := cut(GroupMeanBP/1000, breaks = seq(12,1,-0.5))]

# Read in Param File
param_file = data.table(read.table("param_inputs_clean_filtered.txt", header = T, sep=","))

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
  
for (param in param_file$Learning_Prob_n)
{
  individual_param_run = simulation_data_filtered[Learning_Prob_n == param]

  sim_farming_ancestry = (individual_param_run$Farmer_Ancestry_Partition_All)
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
  
  main_ancestry_dt = rbind(main_ancestry_dt, ancient_data_filtered)
  
}

ggplot(main_ancestry_dt) + geom_point(aes(Learning_Prob_n, sum_log_likelihood))

main_ancestry_dt_best_fit_quad = main_ancestry_dt[sum_log_likelihood > -20000]

quadraticModel = lm(sum_log_likelihood ~ Learning_Prob_n + I(Learning_Prob_n^2), data=main_ancestry_dt_best_fit_quad)

l_doubleprime = -(quadraticModel[["coefficients"]][["I(Learning_Prob_n^2)"]])

confidence_interval = 1.96*(1/(sqrt(l_doubleprime)))

LearningValues = seq(0, 0.006, 0.0001)

sum_log_likelihoodPredict = predict(quadraticModel,list(Learning_Prob_n=LearningValues, Learning_Prob_n2=LearningValues^2))

prediction = LearningValues[match(max(sum_log_likelihoodPredict), sum_log_likelihoodPredict)]

#create scatterplot of original data values
plot(main_ancestry_dt_best_fit_quad$Learning_Prob_n, main_ancestry_dt_best_fit_quad$sum_log_likelihood, pch=16)
#add predicted lines based on quadratic regression model
lines(LearningValues, sum_log_likelihoodPredict, col='blue')

# ----------------------------------------------------------------------------------------------------------------
# Plot data
# ----------------------------------------------------------------------------------------------------------------

# Points with color

ggplot(data = map) + geom_sf(aes(), position = "identity", fill = "grey95", col = "grey95") + theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background = element_rect(fill = "#C3D1D2", color = "white", colour = "white")) + 
  coord_sf(ylim=c(min(ancient_data_filtered$GroupMeanLat, na.rm = T),max(ancient_data_filtered$GroupMeanLat, na.rm = T)), xlim=c(-20,70)) +
  geom_point(data = ancient_data_filtered, aes(GroupMeanLong, GroupMeanLat, col = scaled_weight_Anatolia_N), shape = 16, alpha = 0.99, size=3) +
  guides(col="none") + xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  facet_wrap("time_slice") + geom_point(aes(x=29, y=41), col="red")

# Pie charts

ggplot(data = map) + geom_sf(aes(), position = "identity", fill = "grey95", col = "grey95") + theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background = element_rect(fill = "#C3D1D2", color = "white", colour = "white")) + 
  coord_sf(ylim=c(min(ancient_data_filtered$GroupMeanLat, na.rm = T),max(ancient_data_filtered$GroupMeanLat, na.rm = T)), xlim=c(-20,70)) +
  geom_scatterpie(data = d, aes(GroupMeanLong, GroupMeanLat, alpha=1/se.p1), cols=c("scaled_weight_European_HG", "scaled_weight_Anatolia_N"), color=NA, pie_scale=3) +
  guides(col="none") + xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  facet_wrap("time_slice") #+ geom_point(aes(x=29, y=41), col="red")

ancient_data_filtered[, time_slice_2 := cut(GroupMeanBP/1000, breaks = c(0,5,7,8,20))]
ggplot(ancient_data_filtered[GroupMeanLat<55], aes(distFromAnkara/1000, scaled_weight_Anatolia_N/(scaled_weight_Anatolia_N+scaled_weight_European_HG))) + geom_point() + geom_smooth(n=1000, span=10) + facet_wrap("time_slice_2") + coord_cartesian(ylim=c(0,1))

# Plot of farming ancestry after farming expansion and before Steppe expansion

ggplot(data = map) + geom_sf(aes(), position = "identity", fill = "grey95", col = "grey95") + theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background = element_rect(fill = "#C3D1D2", color = "white", colour = "white")) + 
  coord_sf(ylim=c(min(ancient_data_filtered$GroupMeanLat, na.rm = T),max(ancient_data_filtered$GroupMeanLat, na.rm = T)), xlim=c(-20,70)) +
  geom_scatterpie(data = ancient_data_filtered[GroupMeanLat<55 & GroupMeanBP > 4500 & GroupMeanBP < 8000], aes(GroupMeanLong, GroupMeanLat), cols=c("scaled_weight_European_HG", "scaled_weight_Anatolia_N"), color=NA, pie_scale=0.5) +
  guides(col="none") + xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  facet_wrap("time_slice")

ggplot(d[GroupMeanLat<55 & GroupMeanBP > 4500 & GroupMeanBP < 8000 & GroupNinds > 2], aes(distFromAnkara/1000, weight.p2)) + geom_point() + geom_smooth(n=1000, span=10) + coord_cartesian(ylim=c(0,1), xlim=c(0, 4000)) +
  geom_errorbar(aes(ymin=weight.p2-2*se.p1, ymax=weight.p2+2*se.p1), width=.2) 

ggplot() +
  geom_line(data = simulation_data_filtered, aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = factor(Learning_Prob_n)), linewidth = 1) +
  geom_point(data = main_ancestry_dt, aes(distFromAnkara/1000, scaled_weight_Anatolia_N)) + 
  geom_smooth(n=1000, span=10) + coord_cartesian(ylim=c(0,1), xlim=c(0, 4000)) + theme_bw() +
  labs(y = "Remaining Farming Ancestry") + labs(x = "Distance From Origin (km)")

learning_empirical_overlay = ggplot() +
  geom_line(data = simulation_data_filtered, aes(Mid_Point_km, Farmer_Ancestry_Partition_All, col = factor(Learning_Prob_n)), linewidth = 1.25) +
  geom_point(data = main_ancestry_dt, aes(distFromAnkara/1000, scaled_weight_Anatolia_N)) + 
  geom_smooth(n=1000, span=10) + coord_cartesian(ylim=c(0,1), xlim=c(0, 4000)) + theme_bw() +
  labs(y = "Remaining Farming Ancestry") + labs(x = "Distance From Origin (km)")

ggsave("learning_empirical_overlay.png", plot = learning_empirical_overlay, units = "in", width = 14, height = 8, device="png", dpi=700)

learning_empirical_overlay_best_fit = ggplot() +
  geom_line(data = simulation_data_filtered[Learning_Prob_n == 0.0007], aes(Mid_Point_km, Farmer_Ancestry_Partition_All), col = "#CD9600", linewidth = 1.25) +
  geom_point(data = main_ancestry_dt, aes(distFromAnkara/1000, scaled_weight_Anatolia_N)) + 
  geom_smooth(n=1000, span=10) + coord_cartesian(ylim=c(0,1), xlim=c(0, 4000)) + theme_bw() +
  labs(y = "Remaining Farming Ancestry") + labs(x = "Distance From Origin (km)")

ggsave("learning_empirical_overlay_best_fit.png", plot = learning_empirical_overlay_best_fit, units = "in", width = 14, height = 8, device="png", dpi=700)

