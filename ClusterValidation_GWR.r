# File: ClusterValidation_GWR_SVI.R
# Description: This script performs data preparation and validation for SDoH clusters.
#              It calculates domain-level AHRQ means, merges spatial and cluster data
#              for internal validation, joins SVI data for external validation, and runs
#              Geographically Weighted Regression (GWR) models for both internal and external validation.
# Creator: Yuan Meng

###############################
# 1. Import Required Packages #
###############################
library(sf)
library(scales)
library(dplyr)
library(spgwr)
library(tmap)
library(HDCI)
library(caret)
library(MASS)

##############################################
# 2. Read Data & Calculate Domain-Level Means #
##############################################
years <- c(2009, 2014, 2019)
# For external SVI join, the layers correspond to these years:
year_SVI <- c("2010", "2014", "2020")
domains_external <- c("Social context", "Economic context", "Physical infrastructure")
domains <- c(
  "Social context", "Economic context", "Education",
  "Physical infrastructure", "Healthcare context", "Environment"
)

# Read SDoH categories and get unique domains
categories <- read.csv("AHRQ_directions_12_27.csv")
category_unique <- sort(unique(categories$Main.SDoH.domain))

# Process each year's data to compute domain-level averages
for (yr in years) {
  temp <- read.csv(paste0("AHRQ_scale_fips_", yr, ".csv"))

  # Average indicators for each SDoH domain
  for (j in seq_along(category_unique)) {
    temp_sel <- temp[, categories[categories$Main.SDoH.domain == category_unique[j], "SDoH.variable"]]
    temp_sel <- rowMeans(temp_sel)
    if (j == 1) {
      temp_domain <- temp_sel
    } else {
      temp_domain <- cbind(temp_domain, temp_sel)
    }
  }
  # Append FIPS (assumed to be the last column in the original data)
  temp_domain <- cbind(temp_domain, temp$fips)
  colnames(temp_domain) <- c(
    "AHRQ_Social_context", "AHRQ_Economic_context", "AHRQ_Education",
    "AHRQ_Physical_infrastructure", "AHRQ_Healthcare_context", "AHRQ_Environment",
    "GEOID"
  )
  temp_domain <- as.data.frame(temp_domain)
  assign(paste0("data_", yr), temp_domain)
}
data_list <- list(data_2009, data_2014, data_2019)

##################################################
# 3. Internal Join: Merge Spatial, AHRQ, and Cluster Data #
##################################################
# Read average cluster values
AHRQ_clusters <- read.csv("Chord_mainDomain_summary.csv")

for (yr in years) {
  temp <- data_list[[which(years == yr)]]

  # Read spatial data (counties) for the given year
  counties <- st_read("Counties.gdb", layer = paste0("counties_", yr))
  counties <- counties[, c("GEOID_num", "Shape")]
  temp <- left_join(counties, temp, by = c("GEOID_num" = "GEOID"))

  # Join cluster labels
  cluster_labels <- read.csv(paste0("AHC_scale_fips_labels_", yr, ".csv"))
  cluster_labels <- cluster_labels[, c("fips", "cluster")]
  cluster_labels$cluster <- as.factor(cluster_labels$cluster)
  levels(cluster_labels$cluster) <- c("Cluster 1", "Cluster 2", "Cluster 3")
  temp <- left_join(temp, cluster_labels, by = c("GEOID_num" = "fips"))

  # Join average cluster values by cluster
  AHRQ_clusters_sel <- AHRQ_clusters[AHRQ_clusters$Year == yr, 2:8]
  temp <- left_join(temp, AHRQ_clusters_sel, by = "cluster")

  assign(paste0("Internal_", yr), temp)
}
Internal_list <- list(Internal_2009, Internal_2014, Internal_2019)

###########################################################
# 4. External Join: Merge SVI Data with Cluster Information #
###########################################################
# Re-read AHRQ_clusters if necessary
AHRQ_clusters <- read.csv("Chord_mainDomain_summary.csv")

for (i in seq_along(year_SVI)) {
  # Read SVI data for the corresponding year
  SVI <- st_read("SVI.gdb", layer = paste0("SVI_us_mainland_", year_SVI[i]))
  SVI <- st_zm(SVI)
  temp <- SVI[, c("GEOID_num", "RPL_THEME1", "RPL_THEME2", "RPL_THEME3", "RPL_THEME4")]

  # Use corresponding internal year for cluster labels: assume the order corresponds (2009, 2014, 2019)
  corresponding_year <- years[i]
  cluster_labels <- read.csv(paste0("AHC_scale_fips_labels_", corresponding_year, ".csv"))
  cluster_labels <- cluster_labels[, c("fips", "cluster")]
  cluster_labels$cluster <- as.factor(cluster_labels$cluster)
  levels(cluster_labels$cluster) <- c("Cluster 1", "Cluster 2", "Cluster 3")
  temp <- left_join(temp, cluster_labels, by = c("GEOID_num" = "fips"))

  # Join average cluster values
  AHRQ_clusters_sel <- AHRQ_clusters[AHRQ_clusters$Year == corresponding_year, 2:8]
  temp <- left_join(temp, AHRQ_clusters_sel, by = "cluster")

  assign(paste0("External_", year_SVI[i]), temp)
}
External_list <- list(External_2010, External_2014, External_2020)

#####################################################
# 5. Validation: GWR Modeling (Internal Validation) #
#####################################################
R2_internal <- data.frame()
for (i in seq_along(years)) {
  data <- Internal_list[[i]]

  for (j in seq_along(domains)) {
    # Select independent and dependent variables from the data
    model_data <- data[, c(j + 1, j + 8)]
    colnames(model_data)[1:2] <- c("independent", "dependent")

    formula <- dependent ~ independent
    # Convert sf object to Spatial object for GWR
    model_data <- sf:::as_Spatial(model_data)

    # Select adaptive bandwidth for GWR
    GWRbandwidth <- gwr.sel(formula, data = model_data, adapt = TRUE)

    # Run the GWR model with Gaussian kernel
    gwr_model <- gwr(formula,
      data = model_data,
      adapt = GWRbandwidth,
      hatmatrix = TRUE,
      se.fit = TRUE
    )

    # Compute quasi-global R2
    qGlobalR2 <- 1 - (gwr_model$results$rss / gwr_model$gTSS)
    R2_temp <- data.frame(
      Quasi_global_R2 = qGlobalR2,
      Domain = domains[j],
      Year = years[i]
    )
    R2_internal <- rbind(R2_internal, R2_temp)

    # Export GWR coefficients for the current domain and year
    write.csv(as.data.frame(gwr_model$SDF),
      paste0("./GWR coefficients/GWR_internal_coefficients_", years[i], "_", domains[j], ".csv"),
      row.names = FALSE
    )
  }
}
write.csv(R2_internal, "GWR_internal_R2.csv", row.names = FALSE)

######################################################
# 6. Validation: GWR Modeling (External Validation)  #
######################################################
# Mapping for external validation:
# "Social context" ~ Racial & Ethnic Minority Status (RPL_THEME3)
# "Economic context" ~ Socioeconomic Status (RPL_THEME1)
# "Physical infrastructure" ~ Housing Type & Transportation (RPL_THEME4)
years_external <- c(2009, 2014, 2019)
domains_external <- c("Social context", "Economic context", "Physical infrastructure")
R2_external <- data.frame()

for (i in seq_along(years_external)) {
  data <- External_list[[i]]

  for (j in seq_along(domains_external)) {
    if (j == 1) {
      model_data <- data[, c(4, 7)]
    } else if (j == 2) {
      model_data <- data[, c(2, 8)]
    } else {
      model_data <- data[, c(5, 10)]
    }
    colnames(model_data)[1:2] <- c("independent", "dependent")
    formula <- dependent ~ independent
    model_data <- sf:::as_Spatial(model_data)

    GWRbandwidth <- gwr.sel(formula, data = model_data, adapt = TRUE)
    gwr_model <- gwr(formula,
      data = model_data,
      adapt = GWRbandwidth,
      hatmatrix = TRUE,
      se.fit = TRUE
    )

    qGlobalR2 <- 1 - (gwr_model$results$rss / gwr_model$gTSS)
    R2_temp <- data.frame(
      Quasi_global_R2 = qGlobalR2,
      Domain = domains_external[j],
      Year = years_external[i]
    )
    R2_external <- rbind(R2_external, R2_temp)

    write.csv(as.data.frame(gwr_model$SDF),
      paste0("GWR_external_coefficients_", years_external[i], "_", domains_external[j], ".csv"),
      row.names = FALSE
    )
  }
}
write.csv(R2_external, "GWR_external_R2.csv", row.names = FALSE)
