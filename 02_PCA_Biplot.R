###############################################################################
# Project      : Rice Early Seedling Phenotyping
# Script       : 02_PCA_Biplot.R
#
# Description  : Performs Principal Component Analysis (PCA) on genotype-wise mean early seedling traits measured using the paper towel method.
#
# Objective
#   - Calculate genotype-wise mean values
#   - Perform Principal Component Analysis (PCA)
#   - Visualize relationships among rice genotypes and seedling traits
#   - Generate publication-quality PCA biplot
#
# Input File: Data/Paper_towel_method_final_sheets_Conventional_method.csv
#
# Output Folder: Output/PCA/PCA_Biplot_PaperTowel.tiff
#
# Author       : Priyanka Negi
# Affiliation  : PMRF Fellow | Department of Agricultural Botany (Plant Physiology)
#                MPKV, Rahuri, Maharashtra, India
#                ICAR-NIASM, Baramati, Maharashtra, India
#
# R Version: 4.3+
###############################################################################

##############################
# Load Required Packages
##############################

library(tidyverse)
library(factoextra)

##############################
# Create Output Folder
##############################

if (!dir.exists("Output")) {
  dir.create("Output")
}

if (!dir.exists("Output/PCA")) {
  dir.create("Output/PCA")
}

##############################
# Import Dataset
##############################

phenotype_data <- read.csv(
  "Data/Paper_towel_method_final_sheets_Conventional_method.csv",
  stringsAsFactors = FALSE
)

phenotype_data$Genotype <- as.factor(
  phenotype_data$Genotype
)

##############################
# Calculate Genotype Means
##############################

genotype_means <- phenotype_data %>%

  group_by(Genotype) %>%

  summarise(

    across(

      where(is.numeric),

      mean,

      na.rm = TRUE

    ),

    .groups = "drop"

  )

##############################
# Prepare Data for PCA
##############################

pca_input <- genotype_means %>%

  select(where(is.numeric))

##############################
# Perform Principal Component Analysis
##############################

pca_result <- prcomp(

  pca_input,

  scale. = TRUE

)

##############################
# Generate PCA Biplot
##############################

pca_plot <-

fviz_pca_biplot(

  pca_result,

  label = "var",

  habillage = genotype_means$Genotype,

  repel = TRUE,

  addEllipses = FALSE,

  geom = c("point", "text"),

  pointsize = 3

) +

labs(

  title = "Principal Component Analysis of Rice Genotypes",

  subtitle = "Early seedling traits measured using the Paper Towel Method",

  x = "Principal Component 1",

  y = "Principal Component 2"

) +

theme_minimal(base_size = 14) +

theme(

  plot.title = element_text(

    face = "bold",

    hjust = 0.5

  ),

  plot.subtitle = element_text(

    hjust = 0.5

  ),

  legend.position = "right"

)

##############################
# Save Figure
##############################

ggsave(

  filename = file.path(

    "Output",

    "PCA",

    "PCA_Biplot_PaperTowel.tiff"

  ),

  plot = pca_plot,

  width = 10,

  height = 8,

  dpi = 600,

  compression = "lzw"

)

##############################
# Completion Message
##############################

cat("PCA analysis completed successfully.\n")

cat("Output saved in Output/PCA/\n")

##############################
# Session Information
##############################

sessionInfo()
