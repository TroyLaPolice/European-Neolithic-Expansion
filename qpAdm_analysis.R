#' qpAdm Fixed analysis of European populations to study the expansion
#' of Anatolian-related ancestry through time and space.
#' 2 Source model "European_HG", "Anatolia_N"

#' Load libraries
library(admixtools);library(tidyverse);library(data.table)

#' Set working & output directory
setwd("/path/to/files")
outDIR = ("/path/to/files")
#' Set eigenstrat data prefix 
prefix = "/path/to/files/v52.2_1240K_public/v52.2_1240K_public_qpAdmFarming"
#' Read in metadata & eigenstrat .ind file
meta <- fread("/path/to/files/v52.2_1240K_public/v52.2_1240K_public.anno.txt", header= TRUE)
ind <- fread(paste0(prefix, ".ind"), header= FALSE)

#' Set Right populations
rightPops <- c("Mbuti.DG", "Papuan.DG", "Han.DG", "Karitiana.DG", 
               "Ethiopia_4500BP.SG", "Italy_North_Villabruna_HG",
               "Russia_Ust_Ishim.DG", "Czech_Vestonice", "Russia_MA1_HG.SG",
               "Israel_Natufian")

#' Set Source Populations
#' European_HG == I0585 (LaBrana1), Loschbour.DG, & I1507 (KO1) 
sourcePops <- c("European_HG", "Anatolia_N")

#' Define list of target populations
targetpops <- unique(meta[meta[,FarmingSimsGroups == "y"]]$`Group ID`)
#targetpops <- targetpops[1:5]

#' Check if any popualtions are not in the eigenstrat .ind file
popsToRemove <- list()
for(i in 1:length(targetpops)){
  if(!targetpops[i] %in% ind$V3){
    message(targetpops[i]," population not in the .ind file")
    popsToRemove <- c(popsToRemove, targetpops[i])
  }
}

#' Run qpAdm 
qpAdm_Fixed_2sources_list=list()
qpAdm_Fixed_2sources_summary <- data.frame("target" = as.character(),
                                           "source.p1"= as.character(), "source.p2"= as.character(), 
                                           "weight.p1"= as.numeric(), "weight.p2"= as.numeric(),
                                           "p.value"= as.numeric(), 
                                           "se.p1"= as.numeric(), "se.p2"= as.numeric(),
                                           "f4rank"= as.numeric()) 
elapsed_time <- system.time(for(i in 1:length(targetpops)){
  suppressWarnings({
    message("Running sample  ", targetpops[i]) 
    qpAdm_result <- qpadm(prefix, sourcePops, 
                          rightPops, targetpops[i], allsnps=TRUE, afprod = TRUE, fudge_twice = TRUE)# 
    qpAdm_Fixed_2sources_list[[length(qpAdm_Fixed_2sources_list)+1]] <- qpAdm_result
    qpAdm_Fixed_2sources_summary[nrow(qpAdm_Fixed_2sources_summary)+1, ] <- c(qpAdm_result$weights$target[1],
                                                                              qpAdm_result$weights$left[1],
                                                                              qpAdm_result$weights$left[2],
                                                                              as.numeric(qpAdm_result$weights$weight[1]),
                                                                              as.numeric(qpAdm_result$weights$weight[2]),
                                                                              as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                              as.numeric(qpAdm_result$weights$se[1]),
                                                                              as.numeric(qpAdm_result$weights$se[2]),
                                                                              as.numeric(qpAdm_result$rankdrop$f4rank[1]))
    message(round(i/(length(targetpops)), digits=1) *100, "% Completed")})
})

# Extract the elapsed time (in seconds)
elapsed_seconds <- elapsed_time["elapsed"]
elapsed_seconds
# Save RDS
saveRDS(qpAdm_Fixed_2sources_summary, 
        file = paste0(outDIR, "qpAdm_Fixed_2sources_", "_European_HG_Anatolia_N__Summary", ".rds"))
saveRDS(qpAdm_Fixed_2sources_list, 
        file = paste0(outDIR, "qpAdm_Fixed_2sources_", "_European_HG_Anatolia_N__List", ".rds"))

#' Read RDS
qpAdm_Fixed_2sources_summary = readRDS(file = paste0(outDIR, "qpAdm_Fixed_2sources_", "_European_HG_Anatolia_N__Summary", ".rds"))

#' Fraction of plausible models at 0.05 threshold
nrow(qpAdm_Fixed_2sources_summary %>% filter(p.value >= 0.05 & as.numeric(weight.p1) >= 0 & as.numeric(weight.p1) <= 1 & as.numeric(weight.p2) >= 0 & as.numeric(weight.p2) <= 1)) / nrow(qpAdm_Fixed_2sources_summary)
#' Fraction of plausible models at 0.01 threshold
nrow(qpAdm_Fixed_2sources_summary %>% filter(p.value >= 0.01 & as.numeric(weight.p1) >= 0 & as.numeric(weight.p1) <= 1 & as.numeric(weight.p2) >= 0 & as.numeric(weight.p2) <= 1)) / nrow(qpAdm_Fixed_2sources_summary)

# Only get the individuals included in the qpAdm Analysis
meta_analysis_inds = meta[`Genetic ID` %in% ind[V3 != ""]$V1]


meta_analysis_inds[, GroupMeanBP:= mean(MeanBP), by = qpAdmGroup]
meta_analysis_inds[, GroupMeanLat:= mean(as.numeric(Lat)), by = qpAdmGroup]
meta_analysis_inds[, GroupMeanLong:= mean(as.numeric(Long)), by = qpAdmGroup]
meta_analysis_inds[, GroupNinds:= .N, by = qpAdmGroup]
meta_analysis_inds[, GroupMeanSNPs:= mean(as.numeric(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`)), by = qpAdmGroup]

meta_analysis_inds %>% filter(qpAdmGroup == "Austria_N_LBK") %>% select(GroupMeanBP,GroupMeanLat,GroupMeanLong)


# Merge the two data tables on the ID and sID columns
merged_dt <- unique(merge(qpAdm_Fixed_2sources_summary[, c("target", "source.p1", "source.p2",
                                                    "weight.p1", "weight.p2", "se.p1", "p.value")], 
                          meta_analysis_inds[, c("qpAdmGroup", "GroupMeanBP", "GroupMeanLat", 
                                                 "GroupMeanLong", "GroupNinds", "GroupMeanSNPs")], by.x = "target", by.y = "qpAdmGroup",  all.x = TRUE))
#' filter the .rel inds
#merged_filtered_qpAdmMetaD_2Source = data.table(merged_dt[grepl(".rel.", merged_dt$target), ])
merged_filtered_qpAdmMetaD_2Source = data.table(merged_dt)

#' Define plausible models 
merged_filtered_qpAdmMetaD_2Source[, plausible_p05_weights:= (as.numeric(p.value) >= 0.05 & as.numeric(weight.p1) >= 0 & as.numeric(weight.p1) <= 1 & as.numeric(weight.p2) >= 0 & as.numeric(weight.p2) <= 1)]
merged_filtered_qpAdmMetaD_2Source[, plausible_05p01_weights:= (as.numeric(p.value) >= 0.01 & as.numeric(weight.p1) >= 0 & as.numeric(weight.p1) <= 1 & as.numeric(weight.p2) >= 0 & as.numeric(weight.p2) <= 1)]


#' For non plausilbe models, scaling of admixture weight values to be between 0 an 1
#' First set values < 0  to 0
merged_filtered_qpAdmMetaD_2Source[merged_filtered_qpAdmMetaD_2Source$'weight.p1' < 0]$weight.p1 <- 0
merged_filtered_qpAdmMetaD_2Source[merged_filtered_qpAdmMetaD_2Source$'weight.p2' < 0]$weight.p2 <- 0

merged_filtered_qpAdmMetaD_2Source$weight.p1 <- as.numeric(merged_filtered_qpAdmMetaD_2Source$weight.p1)
merged_filtered_qpAdmMetaD_2Source$weight.p2 <- as.numeric(merged_filtered_qpAdmMetaD_2Source$weight.p2)

#' Generate a normalizing value to scale all weight values to 1: FUN = p * normalizing value
#' normalizing value = (1 / (sum weights))}
merged_filtered_qpAdmMetaD_2Source[,colSum:= ((rowSums(merged_filtered_qpAdmMetaD_2Source[, c('weight.p1', 'weight.p2')], 
                                                                na.rm=TRUE)))] # sum admixture weights for each test
merged_filtered_qpAdmMetaD_2Source[, normalizer:= 1/colSum]
#' Scaling admixture weights from the 2 sources
merged_filtered_qpAdmMetaD_2Source[, scaled_weight_European_HG := mapply(`*`,merged_filtered_qpAdmMetaD_2Source$weight.p1,
                                                                              merged_filtered_qpAdmMetaD_2Source$normalizer)]
merged_filtered_qpAdmMetaD_2Source[, scaled_weight_Anatolia_N := mapply(`*`,merged_filtered_qpAdmMetaD_2Source$weight.p2,
                                                                   merged_filtered_qpAdmMetaD_2Source$normalizer)]


merged_filtered_qpAdmMetaD_2Source$scaled_weight_European_HG <- as.numeric(merged_filtered_qpAdmMetaD_2Source$scaled_weight_European_HG)
merged_filtered_qpAdmMetaD_2Source$scaled_weight_Anatolia_N <- as.numeric(merged_filtered_qpAdmMetaD_2Source$scaled_weight_Anatolia_N)

merged_filtered_qpAdmMetaD_2Source[, scaledWeightsSum := rowSums(.SD), .SDcols = c('scaled_weight_European_HG', 'scaled_weight_Anatolia_N')]
merged_filtered_qpAdmMetaD_2Source$scaledWeightsSum = as.numeric(merged_filtered_qpAdmMetaD_2Source$scaledWeightsSum)



#' Checks & balances
merged_filtered_qpAdmMetaD_2Source %>% filter(plausible_p05_weights == FALSE) %>% dplyr::select(colSum, normalizer, weight.p1, scaled_weight_European_HG, weight.p2, scaled_weight_Anatolia_N)

# Define the file path and name
file_path <- paste0(outDIR, "qpAdm_Fixed__EuroHG_AnatoliaN.txt")

# Use write.table() function to write the data frame to a tab-delimited text file
write.table(merged_filtered_qpAdmMetaD_2Source, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
