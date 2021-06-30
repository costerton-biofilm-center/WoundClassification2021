#=====================================================================
# This script performs normalization, PCA, k-means clustering  
# and differential expression analysis included in the manuscript, 
# "Transcriptomic Classification of Lower Extremity Wounds".
#
# Written by: Blaine Gabriel Fritz
#
#===========================================================================

#=============================================================================
# Set seed, Data Import, and Formatting of Metadata 
#=============================================================================

# Set Seed ===================================================================

# Generate and Set Seed (Based on system time and process id).
# Seed was selected as 15815 randomly (i.e. it was the seed selected)
# when I decided to run the "final" version of the analysis. 

#seed <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^15)
seed <- 15815
set.seed(seed)

# Import Raw Counts, Metadata, and helper functions  ==========================

# Run source to import helper functions
source("./scripts/analysis_utils.R")

# Import counts, metadata and define relative path to annotation file
counts <- get_counts("./data/counts_newdata/")
#counts <- read.csv("./data/Example_data/Example_counts.csv", row.names = 1)
metadata <- get_metadata("./data/Example_data/TEST_metadata_Ranalysis.tsv")
annotation_path <- c("./data/annotation_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff")
kraken_dir = "./data/kraken_newdata"

# Define data Output_directory 
# out_dir needs to contain empty folder called: 
# "normalized_counts", "DESeq2", and "kmeans"

out_dir <- "./analysis_newSamples/"

# Select columns of interest from metadata ====================================

# Select Metadata Columns of Interest.
# Not sure if the analysis will work with 
# continuous/numeric metadata.

factors<-c("Ulcer_duration_cat", "IDSA_SCORE_1to4",
           "Diabetes_type", "PAD", "Acute0_Chronic1","Type", "is_NEBSmallRNA")

# Convert any numeric to factors
metadata[,factors]<- apply(metadata[,factors], 2, as.factor)


# Filtering ===================================================================

# Filter out samples with less than 1M reads
counts <- counts[,colSums(counts) > 1000000]

# Remove the filtered samples from the metadata set
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

# Batch Normalization ========================================================

# Remove HH5 because it's an extreme outlier - Abnormally high counts for e.g. PADI3.
# Remove MW_CW5 and MW_CW6 because they're acute burn wounds - Added back in for evaluation

counts_filtered<-counts[,colnames(counts)!="HH5"]
metadata_filtered<-subset(metadata, metadata$Sample_ID %in% colnames(counts_filtered))

# Only use genes which were used in testing data

unbiased_genes <- read.delim("./analysis/GO_analysis/counts_batchnorm_genelist.txt", header = F)

counts_batchnorm <- counts_filtered[row.names(counts_filtered) %in% unbiased_genes$V1,]

# Only count mRNA from exonic regions

counts_batchnorm_mRNA <- filterCountsbyGeneType(counts_batchnorm, annotation_path, "exon", "mRNA")


#Normalize the counts 
counts_batchnorm_vst <- DESeq2::vst(as.matrix(counts_batchnorm_mRNA))

#======================================================================
# Analyze Kraken data to determine ratio of Bac/human reads 
#======================================================================

library(tidyverse)

kraken_data <-
  list.files(kraken_dir, full.names = T) %>%
  set_names(., nm = map(.x = ., ~gsub(".kraken.report", "", basename(.x)))) %>%
  map(function(x) {
    x <- read.delim(x, header = F)
    x$V6 <- trimws(x$V6, which = "left")
    x <- subset(x, V6 %in% c("root", "Bacteria", "Homo sapiens"))
    return(x)}) %>%
  map_df(~as.data.frame(.x), .id = "Sample_ID") %>%
  select(Sample_ID, "Percent_assigned" = V1, "Total_assigned" = V2, "Organism" = V6) %>%
  group_by(Sample_ID) %>%
  summarise(
    "Total_reads" = Total_assigned[Organism=="root"],
    "Bac_reads" = Total_assigned[Organism=="Bacteria"],
    "Human_reads" = Total_assigned[Organism=="Homo sapiens"],
    "Bac_prcnt" = Percent_assigned[Organism=="Bacteria"],
    "Human_prcnt" = Percent_assigned[Organism=="Homo sapiens"])


#Add kraken_results to kmeans metadata 
metadata_filtered <- dplyr::left_join(metadata_filtered, kraken_data, by="Sample_ID")

#=================================
# Export the data for validation
#=================================

write.csv(counts_batchnorm_vst, "./data/validation_data/ALL_counts_batchnorm_vst.csv")
write.csv(metadata_filtered, "./data/validation_data/ALL_metadata.csv", row.names = F)
