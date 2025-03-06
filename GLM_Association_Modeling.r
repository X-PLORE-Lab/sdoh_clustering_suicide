# File: GLM_Association_Modeling.R
# Description: Performs GLM association modeling for SDoH clusters using four models:
#              Poisson, Negative Binomial, Zero Inflated Poisson (ZIP), and
#              Zero Inflated Negative Binomial (ZINB). The script builds coefficient
#              tables, computes AIC/BIC metrics, calculates confidence intervals,
#              and exports results for further plotting.
# Creator: Yuan Meng

#############################
# 1. Import Required Packages
#############################
library(dplyr)
library(MASS)
library(pscl)



#############################
# 2. Build Coefficient Tables
#############################
# For standard GLM models
GLM_coefficients <- data.frame(matrix(NA, 18, 12))
colnames(GLM_coefficients) <- c(
    "Year", "Reference", "Cluster",
    "all_SD", "gender_male_SD", "gender_female_SD",
    "aged_15_24_SD", "aged_25_64_SD", "aged_over_64_SD",
    "race_Hispanic_SD", "race_Black_SD", "race_White_SD"
)
GLM_coefficients$Year <- c(
    rep("2008-2010", 6),
    rep("2013-2015", 6),
    rep("2018-2020", 6)
)
GLM_coefficients$Reference <- c(
    "Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2", "Cluster 3", "Cluster 3",
    "Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2", "Cluster 3", "Cluster 3",
    "Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2", "Cluster 3", "Cluster 3"
)
GLM_coefficients$Cluster <- c(
    "Cluster 2", "Cluster 3",
    "Cluster 1", "Cluster 3",
    "Cluster 1", "Cluster 2",
    "Cluster 2", "Cluster 3",
    "Cluster 1", "Cluster 3",
    "Cluster 1", "Cluster 2",
    "Cluster 2", "Cluster 3",
    "Cluster 1", "Cluster 3",
    "Cluster 1", "Cluster 2"
)

GLM_coefficients_glm <- GLM_coefficients
GLM_coefficients_glm_nb <- GLM_coefficients

# For Zero Inflated models
GLM_coefficients_ZI <- data.frame(matrix(NA, 36, 12))
colnames(GLM_coefficients_ZI) <- c(
    "Year", "Reference", "Cluster",
    "all_SD", "gender_male_SD", "gender_female_SD",
    "aged_15_24_SD", "aged_25_64_SD", "aged_over_64_SD",
    "race_Hispanic_SD", "race_Black_SD", "race_White_SD"
)
GLM_coefficients_ZI$Year <- c(
    rep("2008-2010", 12),
    rep("2013-2015", 12),
    rep("2018-2020", 12)
)
GLM_coefficients_ZI$Reference <- c(
    rep(c("Cluster 1", "Cluster 1", "Cluster 1", "Cluster 1"), 3),
    rep(c("Cluster 2", "Cluster 2", "Cluster 2", "Cluster 2"), 3),
    rep(c("Cluster 3", "Cluster 3", "Cluster 3", "Cluster 3"), 3)
)
GLM_coefficients_ZI$Cluster <- c(
    "Cluster 2 (count)", "Cluster 3 (count)", "Cluster 2 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 3 (count)", "Cluster 1 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 2 (count)", "Cluster 1 (zero)", "Cluster 2 (zero)",
    "Cluster 2 (count)", "Cluster 3 (count)", "Cluster 2 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 3 (count)", "Cluster 1 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 2 (count)", "Cluster 1 (zero)", "Cluster 2 (zero)",
    "Cluster 2 (count)", "Cluster 3 (count)", "Cluster 2 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 3 (count)", "Cluster 1 (zero)", "Cluster 3 (zero)",
    "Cluster 1 (count)", "Cluster 2 (count)", "Cluster 1 (zero)", "Cluster 2 (zero)"
)
GLM_coefficients_glm_zip <- GLM_coefficients_ZI
GLM_coefficients_glm_zinb <- GLM_coefficients_ZI

#############################
# 3. Build AIC/BIC Table and Define Subgroups
#############################
GLM_AIC_BIC <- data.frame()

group_suicide <- c(
    "all_SD", "gender_male_SD", "gender_female_SD",
    "aged_15_24_SD", "aged_25_64_SD", "aged_over_64_SD",
    "race_Hispanic_SD", "race_Black_SD", "race_White_SD"
)
group_pop <- c(
    "all_POP", "gender_male_POP", "gender_female_POP",
    "aged_15_24_POP", "aged_25_64_POP", "aged_over_64_POP",
    "race_Hispanic_POP", "race_Black_POP", "race_White_POP"
)
group_name <- c(
    "All", "Male", "Female",
    "Aged 15-24", "Aged 25-64", "Aged over 65",
    "Hispanic", "Non-Hispanic Black", "Non-Hispanic White"
)
cluster <- c("Cluster 1", "Cluster 2", "Cluster 3")

coefficients_CI_export_glm_all <- data.frame()
coefficients_CI_export_glm_nb_all <- data.frame()
coefficients_CI_export_glm_zip_all <- data.frame()
coefficients_CI_export_glm_zinb_all <- data.frame()

GLM_invalid <- data.frame()

year <- c("2008-2010", "2013-2015", "2018-2020")

#############################
# 4. Loop Over Years and Subgroups to Fit Models
#############################
for (i in 1:length(year)) {
    # Read the data for the given year
    temp_year <- read.csv(paste0("AHRQ_scale_fips_suicides_pop_", year[i], ".csv"))

    for (j in 1:length(group_suicide)) {
        # Select specific subgroup data
        temp <- temp_year[, c(group_suicide[j], group_pop[j], "cluster")]
        temp <- na.omit(temp) # remove NA data
        colnames(temp) <- c("Suicide", "Population", "Cluster")
        temp <- temp[!temp$Population %in% 0, ] # remove 0 populations
        temp$Cluster <- factor(temp$Cluster, levels = c("Cluster 1", "Cluster 2", "Cluster 3"))

        # Check if over 90% of counties have 0 mortality rates; if so, mark as invalid and skip
        if (nrow(temp[temp$Suicide == 0, ]) >= (nrow(temp) * 0.9)) {
            invalid <- c(year[i], group_suicide[j])
            GLM_invalid <- rbind(GLM_invalid, invalid)
            next
        }

        # Loop over reference groups
        for (m in 1:length(cluster)) {
            temp$Cluster <- relevel(temp$Cluster, ref = cluster[m])

            # --- Model Selection ---
            # Poisson model
            glm_cluster <- glm(Suicide ~ offset(log(Population)) + Cluster,
                data = temp, family = poisson
            )
            glm_aic <- AIC(glm_cluster)
            glm_bic <- BIC(glm_cluster)

            # Negative Binomial model
            glm_nb_cluster <- glm.nb(Suicide ~ offset(log(Population)) + Cluster,
                data = temp, link = log
            )
            glm_nb_aic <- AIC(glm_nb_cluster)
            glm_nb_bic <- BIC(glm_nb_cluster)

            # Zero Inflated Poisson (ZIP) model
            glm_zip_cluster <- zeroinfl(Suicide ~ offset(log(Population)) + Cluster,
                data = temp, dist = "poisson"
            )
            glm_zip_aic <- AIC(glm_zip_cluster)
            glm_zip_bic <- BIC(glm_zip_cluster)

            # Zero Inflated Negative Binomial (ZINB) model
            glm_zinb_cluster <- zeroinfl(Suicide ~ offset(log(Population)) + Cluster,
                data = temp, dist = "negbin"
            )
            glm_zinb_aic <- AIC(glm_zinb_cluster)
            glm_zinb_bic <- BIC(glm_zinb_cluster)

            # --- Calculate and Store AIC/BIC ---
            AIC_BIC <- c(
                year[i], group_name[j], cluster[m],
                round(glm_aic, 2), round(glm_nb_aic, 2),
                round(glm_zip_aic, 2), round(glm_zinb_aic, 2),
                round(glm_bic, 2), round(glm_nb_bic, 2),
                round(glm_zip_bic, 2), round(glm_zinb_bic, 2)
            )
            GLM_AIC_BIC <- rbind(GLM_AIC_BIC, AIC_BIC)

            # --- Calculate Coefficients and Confidence Intervals ---
            result_glm <- summary(glm_cluster)
            result_glm_nb <- summary(glm_nb_cluster)
            result_glm_zip <- summary(glm_zip_cluster)
            result_glm_zinb <- summary(glm_zinb_cluster)

            # Standard GLM
            coefficients_glm <- as.data.frame(round(result_glm$coefficients[2:3, "Estimate"], 2))
            CI_glm <- round(confint(glm_cluster, level = 0.95)[2:3, ], 2)
            rownames(coefficients_glm) <- rownames(CI_glm)

            # Negative Binomial
            coefficients_glm_nb <- as.data.frame(round(result_glm_nb$coefficients[2:3, "Estimate"], 2))
            CI_glm_nb <- round(confint(glm_nb_cluster, level = 0.95)[2:3, ], 2)
            rownames(coefficients_glm_nb) <- rownames(CI_glm_nb)

            # ZIP model (combining count and zero parts)
            coefficients_glm_zip <- as.data.frame(round(rbind(result_glm_zip$coefficients[[1]], result_glm_zip$coefficients[[2]])[c(2:3, 5:6), "Estimate"], 2))
            CI_glm_zip <- round(confint(glm_zip_cluster, level = 0.95)[c(2:3, 5:6), ], 2)
            rownames(coefficients_glm_zip) <- rownames(CI_glm_zip)

            # ZINB model
            coefficients_glm_zinb <- as.data.frame(round(rbind(result_glm_zinb$coefficients[[1]], result_glm_zinb$coefficients[[2]])[c(2:3, 5:6), "Estimate"], 2))
            CI_glm_zinb <- round(confint(glm_zinb_cluster, level = 0.95)[c(2:3, 5:6), ], 2)
            rownames(coefficients_glm_zinb) <- rownames(CI_glm_zinb)

            # --- Export Coefficients (with CI) for Plots ---
            coefficients_CI_export_glm <- cbind(coefficients_glm, CI_glm)
            coefficients_CI_export_glm$year <- year[i]
            coefficients_CI_export_glm$subgroup <- group_name[j]
            coefficients_CI_export_glm$reference <- cluster[m]
            coefficients_CI_export_glm$cluster <- rownames(coefficients_CI_export_glm)

            coefficients_CI_export_glm_nb <- cbind(coefficients_glm_nb, CI_glm_nb)
            coefficients_CI_export_glm_nb$year <- year[i]
            coefficients_CI_export_glm_nb$subgroup <- group_name[j]
            coefficients_CI_export_glm_nb$reference <- cluster[m]
            coefficients_CI_export_glm_nb$cluster <- rownames(coefficients_CI_export_glm_nb)

            coefficients_CI_export_glm_zip <- cbind(coefficients_glm_zip, CI_glm_zip)
            coefficients_CI_export_glm_zip$year <- year[i]
            coefficients_CI_export_glm_zip$subgroup <- group_name[j]
            coefficients_CI_export_glm_zip$reference <- cluster[m]
            coefficients_CI_export_glm_zip$cluster <- rownames(coefficients_CI_export_glm_zip)

            coefficients_CI_export_glm_zinb <- cbind(coefficients_glm_zinb, CI_glm_zinb)
            coefficients_CI_export_glm_zinb$year <- year[i]
            coefficients_CI_export_glm_zinb$subgroup <- group_name[j]
            coefficients_CI_export_glm_zinb$reference <- cluster[m]
            coefficients_CI_export_glm_zinb$cluster <- rownames(coefficients_CI_export_glm_zinb)

            coefficients_CI_export_glm_all <- rbind(coefficients_CI_export_glm_all, coefficients_CI_export_glm)
            coefficients_CI_export_glm_nb_all <- rbind(coefficients_CI_export_glm_nb_all, coefficients_CI_export_glm_nb)
            coefficients_CI_export_glm_zip_all <- rbind(coefficients_CI_export_glm_zip_all, coefficients_CI_export_glm_zip)
            coefficients_CI_export_glm_zinb_all <- rbind(coefficients_CI_export_glm_zinb_all, coefficients_CI_export_glm_zinb)

            # --- Calculate and Format P-values ---
            P_values_glm <- as.data.frame(result_glm$coefficients)
            P_values_glm_nb <- as.data.frame(result_glm_nb$coefficients)
            P_values_glm_zip <- as.data.frame(rbind(result_glm_zip$coefficients[[1]], result_glm_zip$coefficients[[2]]))
            P_values_glm_zinb <- as.data.frame(rbind(result_glm_zinb$coefficients[[1]], result_glm_zinb$coefficients[[2]]))

            # Format p-values for GLM
            for (p in 1:nrow(P_values_glm)) {
                if (is.na(P_values_glm[p, "Pr(>|z|)"])) {
                    P_values_glm[p, "P"] <- NA
                } else if (P_values_glm[p, "Pr(>|z|)"] < 0.001) {
                    P_values_glm[p, "P"] <- "<.001"
                } else {
                    P_values_glm[p, "P"] <- round(P_values_glm[p, "Pr(>|z|)"], 3)
                }
            }
            # Format p-values for Negative Binomial
            for (p in 1:nrow(P_values_glm_nb)) {
                if (is.na(P_values_glm_nb[p, "Pr(>|z|)"])) {
                    P_values_glm_nb[p, "P"] <- NA
                } else if (P_values_glm_nb[p, "Pr(>|z|)"] < 0.0001) {
                    P_values_glm_nb[p, "P"] <- "<0.0001"
                } else {
                    P_values_glm_nb[p, "P"] <- format(round(P_values_glm_nb[p, "Pr(>|z|)"], 4), scientific = FALSE)
                }
            }
            # Format p-values for ZIP
            for (p in 1:nrow(P_values_glm_zip)) {
                if (is.na(P_values_glm_zip[p, "Pr(>|z|)"])) {
                    P_values_glm_zip[p, "P"] <- NA
                } else if (P_values_glm_zip[p, "Pr(>|z|)"] < 0.001) {
                    P_values_glm_zip[p, "P"] <- "<.001"
                } else {
                    P_values_glm_zip[p, "P"] <- round(P_values_glm_zip[p, "Pr(>|z|)"], 3)
                }
            }
            # Format p-values for ZINB
            for (p in 1:nrow(P_values_glm_zinb)) {
                if (is.na(P_values_glm_zinb[p, "Pr(>|z|)"])) {
                    P_values_glm_zinb[p, "P"] <- NA
                } else if (P_values_glm_zinb[p, "Pr(>|z|)"] < 0.001) {
                    P_values_glm_zinb[p, "P"] <- "<.001"
                } else {
                    P_values_glm_zinb[p, "P"] <- round(P_values_glm_zinb[p, "Pr(>|z|)"], 3)
                }
            }

            # Append P-values to coefficient estimates
            P_value_cluster <- P_values_glm[2:3, "P"]
            for (k in 1:nrow(CI_glm)) {
                coefficients_glm[k, 1] <- paste0(
                    coefficients_glm[k, 1], "\n(",
                    CI_glm[k, "2.5 %"], " to ", CI_glm[k, "97.5 %"], ")\nP: ", P_value_cluster[k]
                )
            }
            P_value_cluster <- P_values_glm_nb[2:3, "P"]
            for (k in 1:nrow(CI_glm_nb)) {
                coefficients_glm_nb[k, 1] <- paste0(
                    coefficients_glm_nb[k, 1], "\n(",
                    CI_glm_nb[k, "2.5 %"], " to ", CI_glm_nb[k, "97.5 %"], ")\nP: ", P_value_cluster[k]
                )
            }
            P_value_cluster <- P_values_glm_zip[c(2:3, 5:6), "P"]
            for (k in 1:nrow(CI_glm_zip)) {
                coefficients_glm_zip[k, 1] <- paste0(
                    coefficients_glm_zip[k, 1], "\n(",
                    CI_glm_zip[k, "2.5 %"], " to ", CI_glm_zip[k, "97.5 %"], ")\nP: ", P_value_cluster[k]
                )
            }
            P_value_cluster <- P_values_glm_zinb[c(2:3, 5:6), "P"]
            for (k in 1:nrow(CI_glm_zinb)) {
                coefficients_glm_zinb[k, 1] <- paste0(
                    coefficients_glm_zinb[k, 1], "\n(",
                    CI_glm_zinb[k, "2.5 %"], " to ", CI_glm_zinb[k, "97.5 %"], ")\nP: ", P_value_cluster[k]
                )
            }

            # Reorder coefficients by row names (if needed)
            coefficients_glm <- coefficients_glm[order(rownames(coefficients_glm)), , drop = FALSE]
            coefficients_glm_nb <- coefficients_glm_nb[order(rownames(coefficients_glm_nb)), , drop = FALSE]
            coefficients_glm_zip <- coefficients_glm_zip[order(rownames(coefficients_glm_zip)), , drop = FALSE]
            coefficients_glm_zinb <- coefficients_glm_zinb[order(rownames(coefficients_glm_zinb)), , drop = FALSE]

            # Store coefficients in the main tables for the corresponding subgroup and reference
            GLM_coefficients_glm[GLM_coefficients_glm$Year == year[i] & GLM_coefficients_glm$Reference == cluster[m], group_suicide[j]] <- coefficients_glm
            GLM_coefficients_glm_nb[GLM_coefficients_glm_nb$Year == year[i] & GLM_coefficients_glm_nb$Reference == cluster[m], group_suicide[j]] <- coefficients_glm_nb
            GLM_coefficients_glm_zip[GLM_coefficients_glm_zip$Year == year[i] & GLM_coefficients_glm_zip$Reference == cluster[m], group_suicide[j]] <- coefficients_glm_zip
            GLM_coefficients_glm_zinb[GLM_coefficients_glm_zinb$Year == year[i] & GLM_coefficients_glm_zinb$Reference == cluster[m], group_suicide[j]] <- coefficients_glm_zinb
        } # end loop over reference groups
    } # end loop over subgroups
} # end loop over years

#############################
# 5. Export Results
#############################
# Export invalid groups
colnames(GLM_invalid) <- c("Year", "Subgroup")
write.csv(GLM_invalid, "GLM_invalid_subgroups.csv", row.names = FALSE)

# Export AIC/BIC table (only for reference cluster "Cluster 1")
GLM_AIC_BIC <- GLM_AIC_BIC[GLM_AIC_BIC[, 3] %in% "Cluster 1", ]
colnames(GLM_AIC_BIC) <- c(
    "Year", "Subgroup", "Reference cluster",
    "Poisson AIC", "Negative binomial AIC",
    "Zero-inflated poisson AIC", "Zero-inflated negative binomial AIC",
    "Poisson BIC", "Negative binomial BIC",
    "Zero-inflated poisson BIC", "Zero-inflated negative binomial BIC"
)
write.csv(GLM_AIC_BIC, "GLM_clusters_AIC_BIC_4_models.csv", row.names = FALSE)

# Export coefficient tables
write.csv(GLM_coefficients_glm, "GLM_clusters_coefficients_all_glm.csv", row.names = FALSE)
write.csv(GLM_coefficients_glm_nb, "GLM_clusters_coefficients_all_glm_nb.csv", row.names = FALSE)
write.csv(GLM_coefficients_glm_zip, "GLM_clusters_coefficients_all_glm_zip.csv", row.names = FALSE)
write.csv(GLM_coefficients_glm_zinb, "GLM_clusters_coefficients_all_glm_zinb.csv", row.names = FALSE)

# Set column names for coefficients CI export tables
colnames(coefficients_CI_export_glm_all) <- c(
    "coefficients", "2.5 %", "97.5 %",
    "year", "subgroup", "reference", "cluster"
)
colnames(coefficients_CI_export_glm_nb_all) <- colnames(coefficients_CI_export_glm_all)
colnames(coefficients_CI_export_glm_zip_all) <- colnames(coefficients_CI_export_glm_all)
colnames(coefficients_CI_export_glm_zinb_all) <- colnames(coefficients_CI_export_glm_all)

write.csv(coefficients_CI_export_glm_all, "GLM_clusters_coefficients_all_forPlots_glm.csv", row.names = FALSE)
write.csv(coefficients_CI_export_glm_nb_all, "GLM_clusters_coefficients_all_forPlots_glm_nb.csv", row.names = FALSE)
write.csv(coefficients_CI_export_glm_zip_all, "GLM_clusters_coefficients_all_forPlots_glm_zip.csv", row.names = FALSE)
write.csv(coefficients_CI_export_glm_zinb_all, "GLM_clusters_coefficients_all_forPlots_glm_zinb.csv", row.names = FALSE)
