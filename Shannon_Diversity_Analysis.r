# File: Shannon_Diversity_Analysis.R
# Description: Calculates the Shannon Diversity Index for cluster diversity
#              across US states for years 2009, 2014, and 2019.
# Creator: Yuan Meng

# Load required packages
library(sf)
library(dplyr)
library(vegan)
library(scales)

# Define years to process
years <- c(2009, 2014, 2019)

for (yr in years) {
  # Import spatial data for states
  states <- st_read("Administration data.gdb", layer = "tl_2020_us_mainland_state")
  states$state_num <- as.numeric(states$STATEFP)
  
  # Read clustering labels data
  relabel <- read.csv(paste0("AHC_scale_fips_labels_", yr, ".csv"))
  relabel$STATEFP <- as.integer(relabel$fips / 1000)
  
  # Identify unique states and clusters
  states_unique <- unique(relabel$STATEFP)
  cluster_unique <- sort(unique(relabel$cluster))
  
  # Initialize data frame for cluster counts
  cluster_counts <- data.frame(matrix(NA, nrow = length(states_unique), ncol = length(cluster_unique) + 1))
  colnames(cluster_counts) <- c(paste0("Cluster_", cluster_unique), "STATE")
  
  for (j in seq_along(states_unique)) {
    temp <- subset(relabel, STATEFP == states_unique[j])
    for (k in cluster_unique) {
      count_val <- nrow(subset(temp, cluster == k))
      cluster_counts[j, paste0("Cluster_", k)] <- count_val
    }
    cluster_counts[j, "STATE"] <- states_unique[j]
  }
  
  # Calculate Shannon diversity index for each state based on cluster counts
  shannon_div <- diversity(cluster_counts[, paste0("Cluster_", cluster_unique)], index = "shannon")
  states_div <- data.frame(Diversity = shannon_div, States = cluster_counts$STATE)
  
  # Join cluster counts with diversity data
  states_div_all <- left_join(states_div, cluster_counts, by = c("States" = "STATE"))
  
  # Determine dominant cluster with percentage
  dominant_cluster <- apply(states_div_all[, paste0("Cluster_", cluster_unique)], 1, function(x) {
    dom <- names(x)[which.max(x)]
    perc <- percent(max(x) / sum(x), accuracy = 0.1)
    paste0(dom, " (", perc, ")")
  })
  states_div_all$dominant <- dominant_cluster
  
  # Link state codes to state names
  state_code <- read.csv("State_codes.csv")
  states_div_states <- left_join(states_div_all, state_code, by = c("States" = "Code"))
  
  # Sort states by descending diversity and export summary for the year
  states_div_states <- states_div_states[order(states_div_states$Diversity, decreasing = TRUE), ]
  write.csv(states_div_states, paste0("Shannon_summary_", yr, ".csv"), row.names = FALSE)
}

# Optionally, export the Shannon diversity index from the final iteration
write.csv(states_div, "Shannon_summary.csv", row.names = FALSE)
