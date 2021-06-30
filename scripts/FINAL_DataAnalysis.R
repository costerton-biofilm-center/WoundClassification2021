#=====================================================================
# This script performs normalization, PCA, k-means clustering  
# and differential expression analysis included in the manuscript, 
# "Transcriptomic Classification of Lower Extremity Wounds".
#
# Written by: Blaine Gabriel Fritz
#
#=============================================================================

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
#counts <- get_counts("./data/counts/") #Uncomment if importing counts from a folder with raw count files
counts <- read.csv("./data/Example_data/ALL_counts.csv", row.names = 1) #If importing from a count matrix tsv 
metadata <- get_metadata("./data/Example_data/ALL_metadata_Ranalysis.tsv")
annotation_path <- c("./data/annotation_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff")
kraken_dir = "./data/kraken"

# Define data Output_directory 
# out_dir needs to contain empty folder called: 
# "normalized_counts", "DESeq2", "kmeans", and "validation"

out_dir <- "./analysis/"

# Select columns of interest from metadata ====================================

# Select Metadata Columns which might be of interest.
# Not sure if the analysis will work with 
# continuous/numeric metadata.

factors<-c("Ulcer_duration_cat", "IDSA_SCORE_1to4", "not_healed0_healed1", 
           "Diabetes_type", "PAD", "Acute0_Chronic1", "is_NEBSmallRNA")

# Convert any numeric to factors
metadata[,factors]<- apply(metadata[,factors], 2, as.factor)

#Split data and metadata into main and validation(includes validation data)

#Save full data for train/testing
counts_validation <- counts 
metadata_validation <- metadata

#Subset to only use "Main" data set 
counts <- counts[,c(colnames(counts) %in% metadata$Sample_ID[metadata$Data_set == "Main"])]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))


# Filtering ===================================================================

# Filter out samples with less than 1M reads
counts <- counts[,colSums(counts) > 1000000]

# Remove the filtered samples from the metadata set
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

# =============================================================================
# Differential Gene Expression (DESeq2) of metadata factors
# for "CBC" data (Note: this data is referred to as "LHS" in manuscript.)
# =============================================================================


# Get counts for only CBC data
counts_cbconly <- counts[, colnames(counts) %in% metadata$Sample_ID[metadata$Source=="CBC"]]

# Only count mRNA coding regions on exons
counts_cbconly <- filterCountsbyGeneType(counts_cbconly, annotation_path, "exon", "mRNA")

# Generate the DESeq2 Data set for the CBC data

CBC_DEseq2 <- run_DESeq2(metadata = subset(metadata, Source == "CBC"),
                         counts = counts_cbconly,
                         id_colname = "Sample_ID",
                         metadata_vars = c("Ulcer_duration_cat","IDSA_SCORE_1to4","Acute0_Chronic1"),
                         formula = ~ Ulcer_duration_cat + IDSA_SCORE_1to4 + Acute0_Chronic1)

# Test for differential expression
CBC_DEgenes_UlcerDuration<-DESeq2::results(CBC_DEseq2, contrast = c("Ulcer_duration_cat", "2", "0"))
CBC_DEgenes_IDSAScore<-DESeq2::results(CBC_DEseq2, contrast = c("IDSA_SCORE_1to4", "4","2"))
CBC_DEgenes_Acute0_Chronic1<-DESeq2::results(CBC_DEseq2, contrast = c("Acute0_Chronic1", "1", "0"))

# Filter for genes with abs(log2FoldChange > 2) and padj < 0.05
CBC_DEgenes_Acute0_Chronic1_sig <- filter_DESeq(CBC_DEgenes_Acute0_Chronic1, 2, 0.05)
CBC_DEgenes_IDSAScore_sig <- filter_DESeq(CBC_DEgenes_IDSAScore, 2, 0.05)
CBC_DEgenes_UlcerDuration_sig <- filter_DESeq(CBC_DEgenes_UlcerDuration, 2, 0.05)

#==============================================================================
# Batch Normalization and Differential Gene Expression (DESeq2) analysis 
# of metadata factors for "Combined" data 
#==============================================================================

# Batch Normalization ========================================================

# Remove HH5 because it's an extreme outlier - Abnormally high counts for e.g. PADI3.
# Remove MW_CW5 and MW_CW6 because they're acute burn wounds

counts_filtered<-counts[,colnames(counts)!="HH5" & colnames(counts)!="MW_CW5"& colnames(counts)!="MW_CW6"]
metadata_filtered<-subset(metadata, metadata$Sample_ID %in% colnames(counts_filtered))

# Remove genes differentially expressed based on the "Endedness"
# In the metadata, Endedness is either "PE" or "SE" and only the NEB Small RNA kit is "SE"

counts_batchnorm <- remove_batch_effect(counts_filtered, metadata_filtered, ~is_NEBSmallRNA, 1)


## DESeq2 ======================================================================

# Only count mRNA from exonic regions

counts_batchnorm_mRNA <- filterCountsbyGeneType(counts_batchnorm, annotation_path, "exon", "mRNA")

# Generate the DESeq2 Data set for the combined data set (Note: This excludes the MW and KK data bc they dont have
# both Ulcer_duration_cat and IDSA_Score_1to4 metadata)

combined_DESeq2 <- run_DESeq2(metadata = metadata_filtered,
                              counts = counts_batchnorm_mRNA,
                              id_colname = "Sample_ID",
                              metadata_vars = c("Ulcer_duration_cat","IDSA_SCORE_1to4", "Lib_prep"),
                              formula = ~ Ulcer_duration_cat + IDSA_SCORE_1to4 + Lib_prep)

# Test for differential expression between "extreme" values of metadata
combined_DEgenes_UlcerDuration<-DESeq2::results(combined_DESeq2, contrast = c("Ulcer_duration_cat", "2", "0"))
combined_DEgenes_IDSAScore<-DESeq2::results(combined_DESeq2, contrast = c("IDSA_SCORE_1to4", "4","2"))

# Filter for genes with abs(log2FoldChange > 2) and padj < 0.05
combined_DEgenes_UlcerDuration_sig <- filter_DESeq(combined_DEgenes_UlcerDuration, 2, 0.05)
combined_DEgenes_IDSAScore_sig <- filter_DESeq(combined_DEgenes_IDSAScore, 2, 0.05)

#===================================================================================
# Perform K-means clustering analysis 
#===================================================================================

## CBC data ================================================================

# get CBC metadata 
metadata_cbconly <- subset(metadata, Source=="CBC")

# Normalize the counts 
counts_cbconly_vst <- DESeq2::vst(as.matrix(counts_cbconly))

# Remove genes with zero variance across all samples
novar_filter <- apply(counts_cbconly_vst, 1, sd)
novar_filter <- novar_filter == 0
counts_cbconly_vst <- counts_cbconly_vst[!novar_filter, ]

# Perform the kmeans clustering
set.seed(seed) # Set Seed for reproducibility
kmeans_cbc <- kmeans(t(counts_cbconly_vst), centers = 3, nstart = 25)
print(kmeans_cbc$cluster)
# Get the results 

kmeans_groups_cbc <- data.frame(Sample_ID = names(kmeans_cbc$cluster), cluster_res_cbc = as.character(kmeans_cbc$cluster))

## Combined Data =======================================================================

#Normalize the counts 
counts_batchnorm_vst <- DESeq2::vst(as.matrix(counts_batchnorm_mRNA))

#Remove genes with zero variance 
novar_filter <- apply(counts_batchnorm_vst, 1, sd)
novar_filter <- novar_filter == 0
counts_batchnorm_vst <- counts_batchnorm_vst[!novar_filter, ]
set.seed(seed) # Set Seed for reproducibility
kmeans_all <- kmeans(t(counts_batchnorm_vst), centers = 3, nstart = 25)

# Get the results 

kmeans_groups_all <- data.frame(Sample_ID = names(kmeans_all$cluster), cluster_res_all = as.character(kmeans_all$cluster))

####  Fix kmeans cluster names #####################################
kmeans_groups_cbc$cluster_res_cbc<-
  lapply(kmeans_groups_cbc$cluster_res_cbc, function(cluster){
    if(cluster=="1"){
      cluster <- "3"
    }
    else if(cluster=="2"){
      cluster <- "1"
    }
    else if(cluster=="3"){
      cluster<-"2"
    }
    else{
      stop("ERROR CONVERTING KMEANS GROUPS!!")
    }
  })

kmeans_groups_cbc$cluster_res_cbc <- unlist(kmeans_groups_cbc$cluster_res_cbc)

kmeans_groups_all$cluster_res_all<-
  lapply(kmeans_groups_all$cluster_res_all, function(cluster){
    if(cluster=="1"){
      cluster <- "3"
    }
    else if(cluster=="2"){
      cluster <- "1"
    }
    else if(cluster=="3"){
      cluster<-"2"
    }
    else{
      stop("ERROR CONVERTING KMEANS GROUPS!!")
    }
  })

kmeans_groups_all$cluster_res_all <- unlist(kmeans_groups_all$cluster_res_all)

#=======================================================
# DESEQ2 of k-means Results 
#=======================================================

# Add kmeans results to metadata 
metadata_kmeans <- dplyr::left_join(metadata, kmeans_groups_all, by="Sample_ID")
metadata_kmeans <- dplyr::left_join(metadata_kmeans, kmeans_groups_cbc, by="Sample_ID")

# Run DESeq2
dds_kmeans <- run_DESeq2(metadata = metadata_kmeans, 
                         counts = counts_batchnorm_mRNA,
                         id_colname = "Sample_ID",
                         metadata_vars = c("cluster_res_all", "Lib_prep"),
                         formula = ~ cluster_res_all + Lib_prep) 


# Extract DE genes results
DEseq_kmeans_all_3v2 <- DESeq2::results(dds_kmeans, contrast = c("cluster_res_all", "3", "2"))
DEseq_kmeans_all_1v2 <- DESeq2::results(dds_kmeans, contrast = c("cluster_res_all", "1", "2"))
DEseq_kmeans_all_3v1 <- DESeq2::results(dds_kmeans, contrast = c("cluster_res_all", "3", "1"))

# Extract significant genes
DEseq_kmeans_all_3v2_sig <- filter_DESeq(DEseq_kmeans_all_3v2, 2, 0.05)
DEseq_kmeans_all_1v2_sig <- filter_DESeq(DEseq_kmeans_all_1v2, 2, 0.05)
DEseq_kmeans_all_3v1_sig <- filter_DESeq(DEseq_kmeans_all_3v1, 2, 0.05)


#=========================================================
# Count # of differentially expressed genes from DESeq 
#=========================================================

DESeq_summary<-
  summarize_DESeq(list(combined_DEgenes_UlcerDuration, #2vs0
                       combined_DEgenes_IDSAScore, #4vs2
                       DEseq_kmeans_all_3v2,
                       DEseq_kmeans_all_1v2,
                       DEseq_kmeans_all_3v1,
                       CBC_DEgenes_UlcerDuration, #2vs0
                       CBC_DEgenes_IDSAScore, #4vs2 
                       CBC_DEgenes_Acute0_Chronic1
  ),
  names = c("combined_DEgenes_UlcerDuration", #2vs0
            "combined_DEgenes_IDSAScore", #4vs2
            "DEseq_kmeans_all_3v2",
            "DEseq_kmeans_all_1v2",
            "DEseq_kmeans_all_3v1",
            "CBC_DEgenes_UlcerDuration", #2vs0
            "CBC_DEgenes_IDSAScore", #4vs2 
            "CBC_DEgenes_Acute0_Chronic1"))#1vs0

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
metadata_kmeans <- dplyr::left_join(metadata_kmeans, kraken_data, by="Sample_ID")

#=========================
# Analyze Validation Data 
#=========================

# Filtering ====================

# Filter out samples with less than 1M reads
counts_validation <- counts_validation[,colSums(counts_validation) > 1000000]
# Remove the filtered samples from the metadata set
metadata_validation<-subset(metadata_validation, metadata_validation$Sample_ID %in% colnames(counts_validation))

# Remove HH5 because it's an extreme outlier - Abnormally high counts for e.g. PADI3.
# Including MW_CW5 and MW_CW6 because they're acute burn wounds - Added back in for evaluation

counts_validation_filtered<-counts_validation[,colnames(counts_validation)!="HH5"]
metadata_validation_filtered<-subset(metadata_validation, metadata$Sample_ID %in% colnames(counts_validation_filtered))

# Normalization ====================

# Only use genes which were used in main data

unbiased_genes <- row.names(counts_batchnorm)
counts_validation_batchnorm <- counts_validation_filtered[row.names(counts_validation_filtered) %in% unbiased_genes,]

# Only count mRNA from exonic regions
counts_validation_batchnorm_mRNA <- filterCountsbyGeneType(counts_validation_batchnorm, annotation_path, "exon", "mRNA")

#Normalize the counts 
counts_validation_batchnorm_vst <- DESeq2::vst(as.matrix(counts_validation_batchnorm_mRNA))

#Add Kraken Data
metadata_validation_filtered <- dplyr::left_join(metadata_validation_filtered, kraken_data, by="Sample_ID")


#=========================================================================
# Export gene counts data tables and analysis results to output directory
#
# WARNING!!! Will overwrite existing files without asking
#=========================================================================


#Export Seed and Session Info 
writeLines(as.character(seed), paste0(out_dir,"seed.txt"))
writeLines(capture.output(sessionInfo()), paste0(out_dir,"sessionInfo.txt"))

# Export Counts
write.csv(counts_batchnorm_vst, paste0(out_dir,"normalized_counts/Main_counts_batchnorm_vst.csv"))
write.csv(counts_cbconly_vst, paste0(out_dir,"normalized_counts/counts_cbconly_vst.csv"))

# Export DESeq results
write.csv(CBC_DEgenes_Acute0_Chronic1_sig, paste0(out_dir,"DESeq2/CBC_DEgenes_Acute0_Chronic1.csv"))
write.csv(CBC_DEgenes_IDSAScore_sig,paste0(out_dir,"DESeq2/CBC_DEgenes_IDSAScore.csv"))
write.csv(CBC_DEgenes_UlcerDuration_sig,paste0(out_dir,"DESeq2/CBC_DEgenes_UlcerDuration.csv" ))
write.csv(combined_DEgenes_UlcerDuration_sig,paste0(out_dir,"DESeq2/combined_DEgenes_UlcerDuration.csv" ))
write.csv(combined_DEgenes_IDSAScore_sig,paste0(out_dir,"DESeq2/combined_DEgenes_IDSAScore.csv" ))
write.csv(DEseq_kmeans_all_3v2_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_3v2_sig.csv"))
write.csv(DEseq_kmeans_all_1v2_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_1v2_sig.csv"))
write.csv(DEseq_kmeans_all_3v1_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_3v1_sig.csv"))

#Export DESeq summary

write.csv(DESeq_summary, paste0(out_dir,"DESeq2/n_sig_DEgenes.csv"))

# Export kmeans results
write.csv(kmeans_groups_all, paste0(out_dir, "kmeans/kmeans_groups_all.csv"), row.names = FALSE)
write.csv(kmeans_groups_cbc, paste0(out_dir, "kmeans/kmeans_groups_cbc.csv"), row.names = FALSE)

# Export Data for GO_analysis

write.table(row.names(subset(DEseq_kmeans_all_3v2_sig,log2FoldChange > 2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_3v2_sigUP.txt"))
write.table(row.names(subset(DEseq_kmeans_all_3v2_sig,log2FoldChange < -2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_3v2_sigDOWN.txt"))
write.table(row.names(subset(DEseq_kmeans_all_3v1_sig,log2FoldChange > 2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_3v1_sigUP.txt"))
write.table(row.names(subset(DEseq_kmeans_all_3v1_sig,log2FoldChange < -2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_3v1_sigDOWN.txt"))
write.table(row.names(subset(DEseq_kmeans_all_1v2_sig,log2FoldChange > 2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_1v2_sigUP.txt"))
write.table(row.names(subset(DEseq_kmeans_all_1v2_sig,log2FoldChange < -2)),
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_all_1v2_sigDOWN.txt"))
write.table(row.names(counts_batchnorm),
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/counts_batchnorm_genelist.txt"))

#Export validation Data analysis

write.csv(counts_validation_batchnorm_vst, "./data/Example_data/Validation_counts_batchnorm_vst.csv")
write.csv(metadata_validation_filtered, "./data/Example_data/Validation_metadata.csv", row.names = F)






