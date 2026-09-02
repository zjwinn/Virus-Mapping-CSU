print("Starting pipeline...")

print("Step 1: Purling main analysis...")
source("purling_script.R")

print("Step 2: Running main analysis...")
source("main_analysis.R")

print("Step 3: Running GWAS in each environment...")
source("gwas_in_each_environment.R")

print("Step 4: Running GWAS interpretation...")
source("gwas_interpretation.R")

print("Step 5: Generating virus frequency images...")
source("virus_frequency_images.R")

print("Step 6: Visualizing PCA for publication...")
source("visualizing_PCA_for_publication.R")

print("Step 7: Running predictive haplotyping...")
source("predictive_haplotyping.R")

if (file.exists(file.path(getwd(), "manuscript"))) {
  print("Step 8: Collecting and formatting images for publication...")
  source(file.path(getwd(), "manuscript/Gathering_Tables_and_Images.R"))
}

print("Pipeline completed successfully!")
