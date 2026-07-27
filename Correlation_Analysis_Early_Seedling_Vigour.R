
###############################################################
# Project : Correlation Analysis of Early Seedling Vigour Traits
# Crop    : Rice (Oryza sativa L.)
#
# Description:
# This script performs Pearson correlation analysis among
# mesocotyl length (ML), plumule length (PL), and root length (RL)
# measured at 7, 10, and 14 days after sowing (DAS).
#
# Outputs:
# 1. Pearson correlation matrices
# 2. Correlation coefficients, p-values and t-statistics
# 3. Pairwise correlation plots
# 4. Correlation heatmaps
# 5. Summary table of correlation coefficients
#
# NOTE:
# Trait measurements were obtained through conventional
# phenotyping and ImageJ-based image analysis.
#
# Author      : Priyanka Negi
# Affiliation : Mahatma Phule Krishi Vidyapeeth (MPKV), Rahuri
#               ICAR-NIASM, Baramati
# Year        : 2025
###############################################################

###############################################################
# 1. Load Required Packages
###############################################################

library(psych)
library(corrplot)
library(PerformanceAnalytics)
library(dplyr)

###############################################################
# 2. Import Dataset
###############################################################

seedling_data <- read.csv(
  "Data/Exp_1_paper_towel_method.csv",
  stringsAsFactors = FALSE
)

str(seedling_data)

###############################################################
# 3. Data Preparation
###############################################################

traits <- c(
  "ML7","PL7","RL7",
  "ML10","PL10","RL10",
  "ML14","PL14","RL14"
)

seedling_data[traits] <- lapply(
  seedling_data[traits],
  function(x) as.numeric(as.character(x))
)

# Check missing observations
colSums(is.na(seedling_data[traits]))

# Trait means
colMeans(seedling_data[traits], na.rm = TRUE)

###############################################################
# 4. Pearson Correlation Analysis
###############################################################

corr_results <- corr.test(
  seedling_data[traits],
  use = "pairwise",
  method = "pearson",
  adjust = "none"
)

corr_results$r
corr_results$t
corr_results$p
corr_results$se

dir.create("Output", showWarnings = FALSE)

sink("Output/correlation_output.txt")
print(corr_results)
sink()

###############################################################
# 5. Pairwise Correlation Plot
###############################################################

png(
  "Output/Pairwise_Correlation_Plots.png",
  width = 3000,
  height = 2500,
  res = 300
)

pairs.panels(
  seedling_data[traits],
  hist.col = "cyan",
  stars = TRUE,
  method = "pearson",
  ellipses = TRUE
)

dev.off()

###############################################################
# 6. Correlation Matrix
###############################################################

tiff(
  "Output/Correlation_Matrix.tiff",
  width = 10,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw"
)

chart.Correlation(
  seedling_data[traits],
  histogram = TRUE,
  pch = 19
)

dev.off()

###############################################################
# 7. Trait-wise Correlation by Day
###############################################################

day7  <- seedling_data[, c("ML7","PL7","RL7")]
day10 <- seedling_data[, c("ML10","PL10","RL10")]
day14 <- seedling_data[, c("ML14","PL14","RL14")]

day7[]  <- lapply(day7, as.numeric)
day10[] <- lapply(day10, as.numeric)
day14[] <- lapply(day14, as.numeric)

cor_7  <- cor(day7, use = "pairwise.complete.obs")
cor_10 <- cor(day10, use = "pairwise.complete.obs")
cor_14 <- cor(day14, use = "pairwise.complete.obs")

###############################################################
# 8. Correlation Heatmaps
###############################################################

png(
  "Output/Correlation_Comparison_Days.png",
  width = 3600,
  height = 1200,
  res = 300
)

par(mfrow = c(1,3))

corrplot(cor_7, method="color", type="upper",
         title="7 DAS", mar=c(0,0,2,0))
corrplot(cor_10, method="color", type="upper",
         title="10 DAS", mar=c(0,0,2,0))
corrplot(cor_14, method="color", type="upper",
         title="14 DAS", mar=c(0,0,2,0))

dev.off()

###############################################################
# 9. Correlation Significance
###############################################################

corr_day7  <- corr.test(day7)
corr_day10 <- corr.test(day10)
corr_day14 <- corr.test(day14)

png(
  "Output/Correlation_Values_Significance.png",
  width=3600,
  height=1200,
  res=300
)

par(mfrow=c(1,3))

corrplot(corr_day7$r, method="number", type="upper",
         p.mat=corr_day7$p, sig.level=0.05,
         insig="blank", number.cex=1.2,
         title="7 DAS")

corrplot(corr_day10$r, method="number", type="upper",
         p.mat=corr_day10$p, sig.level=0.05,
         insig="blank", number.cex=1.2,
         title="10 DAS")

corrplot(corr_day14$r, method="number", type="upper",
         p.mat=corr_day14$p, sig.level=0.05,
         insig="blank", number.cex=1.2,
         title="14 DAS")

dev.off()

###############################################################
# 10. Summary Table
###############################################################

cor_table <- data.frame(
  Day = rep(c("7 DAS","10 DAS","14 DAS"), each = 3),
  Trait_Pair = rep(
    c("Mesocotyl–Plumule",
      "Mesocotyl–Root",
      "Plumule–Root"),
    3
  ),
  Pearson_r = c(
    cor_7[1,2], cor_7[1,3], cor_7[2,3],
    cor_10[1,2], cor_10[1,3], cor_10[2,3],
    cor_14[1,2], cor_14[1,3], cor_14[2,3]
  )
)

print(round(cor_table,3))

write.csv(
  cor_table,
  "Output/Correlation_Values_Table.csv",
  row.names = FALSE
)

###############################################################
# 11. Genotype Mean Correlations
###############################################################

genotype_means <- seedling_data %>%
  group_by(Genotype) %>%
  summarise(across(ML7:RL14, mean, na.rm = TRUE))

cor(genotype_means[,c("ML7","PL7","RL7")])
cor(genotype_means[,c("ML10","PL10","RL10")])
cor(genotype_means[,c("ML14","PL14","RL14")])

###############################################################
# 12. Session Information
###############################################################

sessionInfo()
