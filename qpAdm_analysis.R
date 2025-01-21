#' qpAdm Fixed analysis of European populations to study the expansion
#' of Anatolian-related ancestry through time and space.
#' 3 Source Steppe_EBA_Yamnaya, European_HG & Anatolia_N

#' Load libraries
library(admixtools);library(tidyverse);library(data.table)

#' Set working & output directory
setwd("/path/to/files")
outDIR = ("/path/to/files")
#' Set eigenstrat data prefix 
prefix = "/path/to/files/v62.0_AADR/v62.0_1240k_public__qpadm_ind"
#' Read in metadata & eigenstrat .ind file
ind <- fread(paste0(prefix, ".ind"), header= FALSE)


#' Set Right populations
rightPops <- c("Russia_Afanasievo","WHGB", "Turkey_N", "OldAfrica")


#' Set Source Populations
sourcePops <- c("OldSteppe", "WHGA", "EEF")

#' Define list of target populations
targetpops <- ind[!(is.na(V3))]$V3

targetpops = targetpops[!(targetpops %in% c("Russia_Afanasievo","WHGB", "Turkey_N", "OldAfrica", "OldSteppe", "WHGA", "EEF"))]
#targetpops = targetpops[2:5]


# Load necessary packages
library(foreach)
library(doParallel)

# Set up the parallel backend to use multiple cores
numCores <- parallel::detectCores() - 1  # Use three less than the total cores
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Start the timer
elapsed_time <- system.time({
  
  # Use foreach loop to run iterations in parallel
  results <- foreach(i = 1:length(targetpops), .packages = c("admixtools")) %dopar% {

    # Suppress warnings and messages within each iteration
    suppressWarnings({
      # Run the qpAdm function for the current target population
      qpAdm_result <- qpadm(
        prefix,
        sourcePops,
        rightPops,
        targetpops[i],
        allsnps = TRUE,
        afprod = TRUE
      )
      
      # Create a summary row for the current result
      summary_row <- data.frame(
        target = qpAdm_result$weights$target[1],
        source.p1 = qpAdm_result$weights$left[1],
        source.p2 = qpAdm_result$weights$left[2],
        source.p3 = qpAdm_result$weights$left[3],
        weight.p1 = as.numeric(qpAdm_result$weights$weight[1]),
        weight.p2 = as.numeric(qpAdm_result$weights$weight[2]),
        weight.p3 = as.numeric(qpAdm_result$weights$weight[3]),
        p.value = as.numeric(qpAdm_result$rankdrop$p[1]),
        se.p1 = as.numeric(qpAdm_result$weights$se[1]),
        se.p2 = as.numeric(qpAdm_result$weights$se[2]),
        se.p3 = as.numeric(qpAdm_result$weights$se[3]),
        f4rank = as.numeric(qpAdm_result$rankdrop$f4rank[1]),
        stringsAsFactors = FALSE
      )
      
      # Return both the qpAdm_result and the summary_row as a list
      list(
        qpAdm_result = qpAdm_result,
        summary_row = summary_row
      )
    })
  }
  
})

# After the parallel processing, collect and organize the results
# Since we used a list output in foreach, 'results' is a list of lists

# Separate the qpAdm results and the summary rows
qpAdm_Fixed_3sources_list <- lapply(results, function(x) x$qpAdm_result)
qpAdm_Fixed_3sources_summary <- do.call(rbind, lapply(results, function(x) x$summary_row))

# Extract the elapsed time (in seconds)
elapsed_seconds <- elapsed_time["elapsed"]
print(elapsed_seconds)

# Stop and clean up the parallel cluster
stopCluster(cl)
registerDoSEQ()  # Register sequential backend to prevent issues later


# Generate plausible conditions for 
qpAdm_Fixed_3sources_summary = data.table(qpAdm_Fixed_3sources_summary)

qpAdm_Fixed_3sources_summary[(weight.p1 > 0 & weight.p1 < 1) & (weight.p2 > 0 & weight.p2 < 1) & (weight.p3 > 0 & weight.p3 < 1) := weights_01 = "Y"]
qpAdm_Fixed_3sources_summary[p.value >= 0.01 := p001 = "Y"]
qpAdm_Fixed_3sources_summary[p.value >= 0.05 := p005 = "Y"]
qpAdm_Fixed_3sources_summary[(weight.p1 > 0 & weight.p1 < 1) & (weight.p2 > 0 & weight.p2 < 1) & (weight.p3 > 0 & weight.p3 < 1) & p.value >= 0.05 := weights_01_p005 = "Y"]
qpAdm_Fixed_3sources_summary[(weight.p1 > 0 & weight.p1 < 1) & (weight.p2 > 0 & weight.p2 < 1) & (weight.p3 > 0 & weight.p3 < 1) & p.value >= 0.01 := weights_01_p001 = "Y"]
qpAdm_Fixed_3sources_summary[(weight.p1 > 0 + (2*se.p1) & weight.p1 < 1 - (2*se.p1)) & (weight.p2 > 0 + (2*se.p2) & weight.p2 < 1 - (2*se.p2)) & (weight.p3 > 0 + (2*se.p3) & weight.p3 < 1- (2*se.p3)) & p.value >= 0.01 := weights_01.2SE_p001 = "Y"]
qpAdm_Fixed_3sources_summary[(weight.p1 > 0 + (2*se.p1) & weight.p1 < 1 - (2*se.p1)) & (weight.p2 > 0 + (2*se.p2) & weight.p2 < 1 - (2*se.p2)) & (weight.p3 > 0 + (2*se.p3) & weight.p3 < 1- (2*se.p3)) & p.value >= 0.05 := weights_01.2SE_p005 = "Y"]
qpAdm_Fixed_3sources_summary[p.value >= 0.01 & (se.p1 < 0.022 & se.p2 < 0.022 & se.p3 < 0.022) := SE002_p001 = "Y"]
qpAdm_Fixed_3sources_summary[p.value >= 0.05 & (se.p1 < 0.022 & se.p2 < 0.022 & se.p3 < 0.022) := SE002_p005 = "Y"]
qpAdm_Fixed_3sources_summary[sum(weight.p1, weight.p2, weight.p3) := weights_sum]




# Assuming qpAdm_Fixed_3sources_summary is already a data.table

# Initialize new columns with "N"
qpAdm_Fixed_3sources_summary[, `:=`(
  weights_01 = "N",
  p001 = "N",
  p005 = "N",
  weights_01_p001 = "N",
  weights_01_p005 = "N",
  weights_01_2SE_p001 = "N",
  weights_01_2SE_p005 = "N",
  SE002_p001 = "N",
  SE002_p005 = "N"
)]

# Update columns based on conditions
qpAdm_Fixed_3sources_summary[
  (weight.p1 > 0 & weight.p1 < 1) &
    (weight.p2 > 0 & weight.p2 < 1) &
    (weight.p3 > 0 & weight.p3 < 1),
  weights_01 := "Y"
]

qpAdm_Fixed_3sources_summary[
  p.value >= 0.01,
  p001 := "Y"
]
