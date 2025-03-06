# File: SDOH_Clustering.R
# Description: Preprocesses SDOH data for years 2009, 2014, and 2019,
#              applies transformations and scaling, and performs
#              Agglomerative Hierarchical Clustering.
# Creator: Yuan Meng

# Load required packages
library(cluster)
library(NbClust)
library(scales)

# Define cube root function (defined once for all iterations)
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# Define years to process
years <- c(2009, 2014, 2019)

# Read SDoH categories with directions
categories <- read.csv("AHRQ_directions_12_27.csv")

# Preprocess and scale data for each year
for (yr in years) {
  # Read SDoH indicators
  file_path <- paste0("AHRQ_shared_", yr, ".csv")
  temp <- read.csv(file_path)
  temp <- temp[, 2:285]
  
  # Reverse sign for variables with negative direction
  neg_vars <- categories[categories$Directions == 'Negative', "SDoH.variable"]
  temp[, neg_vars] <- -temp[, neg_vars]
  
  # Apply cube root transformation
  temp <- Math.cbrt(temp)
  
  # Scale data to range [0, 1]
  temp <- data.frame(lapply(temp, rescale, to = c(0, 1)))
  
  # Add county FIPS from the original file
  FIPS <- read.csv(file_path)
  temp$fips <- FIPS$COUNTYFIPS
  
  # Export scaled data
  write.csv(temp, paste0("AHRQ_scale_fips_", yr, ".csv"), row.names = FALSE)
}

# Agglomerative Hierarchical Clustering
AHC_best_num <- data.frame()

for (yr in years) {
  # Read scaled data
  data <- read.csv(paste0("AHRQ_scale_fips_", yr, ".csv"))
  data <- data[, -1]
  data_ahc <- data[, 1:284]
  
  # Determine optimal clusters using NbClust
  Cluster_num <- NbClust(data_ahc, distance = "euclidean", method = "ward.D2")
  best_num <- length(unique(Cluster_num$Best.partition))
  AHC_best_num <- rbind(AHC_best_num, data.frame(Year = yr, Best_Num = best_num))
  
  # Perform Agglomerative Hierarchical Clustering
  ahc <- agnes(data_ahc, method = "ward")
  
  # Save the clustering model
  saveRDS(ahc, paste0("AHC_model_", yr, ".rda"))
  
  # Plot dendrogram and save as a TIFF file (unique file name per year)
  tiff(filename = paste0("dendrogram_", yr, ".tif"), width = 4300, height = 1500, res = 300)
  plot(ahc)
  rect.hclust(as.hclust(ahc), k = 3, border = 2:6)
  dev.off()
  
  # Cut the tree to form clusters based on best_num
  clusters <- cutree(as.hclust(ahc), k = best_num)
  clusters <- data.frame(cluster = clusters)
  
  # Combine scaled data with cluster labels and export
  AHRQ_AHC <- cbind(data, cluster = clusters$cluster)
  write.csv(AHRQ_AHC, paste0("AHC_scale_fips_labels_", yr, ".csv"), row.names = FALSE)
  
  # Calculate and export cluster frequency counts
  label_count <- as.data.frame(table(clusters$cluster))
  write.csv(label_count, paste0("AHC_labels_counts_", yr, ".csv"), row.names = FALSE)
}

# Export best cluster numbers for each year
write.csv(AHC_best_num, "AHC_model_best_num.csv", row.names = FALSE)
