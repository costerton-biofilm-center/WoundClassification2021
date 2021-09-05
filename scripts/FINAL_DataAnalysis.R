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
counts <- read.csv("./data/Example_data/ALL_counts.csv", row.names = 1, check.names = F, comment.char = "") #If importing from a count matrix tsv 
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

factors<-c("Ulcer_duration_cat", "IDSA_SCORE_1to4", "is_NEBSmallRNA")

# Convert any numeric to factors
metadata[,factors]<- apply(metadata[,factors], 2, as.factor)

# Filtering ===================================================================

# Filter out samples with less than 1M reads and remove from metata
counts <- counts[,colSums(counts) > 1000000]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

# Remove MW_CW5 and MW_CW6, HH5 because it's an extreme outlier - Abnormally high counts for e.g. PADI3.

counts<-counts[,!(colnames(counts) %in% c("HH5","MW_CW5","MW_CW6"))]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

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
metadata <- dplyr::left_join(metadata, kraken_data, by="Sample_ID")

#Add Category for high bacteria to kmeans metadata 

metadata$high_bacteria <- ifelse(metadata$Bac_prcnt>10, "high", "low")

#=========================================================================
# Analyze species and genera abundance from Kraken
#=========================================================================

#Note: Kraken files need to include Eukaryota, Bacteria, Viruses, and Archaea

kraken_bac <- 
  list.files(kraken_dir, full.names = T) %>%
  set_names(., nm = map(.x = ., ~gsub(".kraken.report", "", basename(.x)))) %>%
  map(function(x) {
    #### Read in the files and format a bit 
    x <- read.delim(x, header = F)
    x$V6 <- trimws(x$V6, which = "left")
    x <- x[, c(1,2,4,6)] #"Subset columns of interest"
    colnames(x)<-c("Prcnt_root_taxon", "nr_reads_root_taxon", "Rank_code", "ID")
    
    #### We need to extract the bacteria info, but the order of the domains in the file are not consistent
    #### So here is some ugly code to get only bacterial species
    
    domain_indices <- subset(x, Rank_code == "D") #get indexes of domain ranks 
    domain_indices$index <-c(1:nrow(domain_indices)) #make a vector for indexing
    
    
    bac_index <- domain_indices$index[domain_indices$ID =="Bacteria"] #Find which row is bacteria 
    next_domain <- bac_index+1
    
    bac_index <- as.numeric(row.names(domain_indices)[bac_index])
    
    if(next_domain == nrow(domain_indices)){
      stop_index <- nrow(x)
    }
    else{ 
      stop_index <-as.numeric(row.names(domain_indices)[next_domain])-1
    }
    x <- x[bac_index:stop_index,]
    x <- subset(x, Rank_code == "S")
    return(x)}) %>%
  bind_rows(., .id = "Sample_ID")


kraken_bac_abundance <-
  metadata[,c("Sample_ID","Bac_reads", "Human_reads")] %>%
  left_join(kraken_bac, by = "Sample_ID") %>%
  group_by(Sample_ID) %>%
  mutate("Rel_abund" = nr_reads_root_taxon/sum(nr_reads_root_taxon)*100)



# Get only samples with high bacteria 

high_bac_microbiome <- 
  subset(kraken_bac_abundance, Sample_ID %in% metadata$Sample_ID[metadata$Bac_prcnt>30]) 

high_bac_microbiome_top_species <- 
  high_bac_microbiome %>% 
  group_by(Sample_ID) %>%
  arrange(Sample_ID, desc(Rel_abund)) %>%
  filter(Rel_abund>1)

#=====================================================================================
#Split data and metadata into Main and Validation(Main + Validation)
#=====================================================================================
counts_validation <- counts 
metadata_validation <- metadata

#Subset to only use "Main" data set 
counts <- counts[,c(colnames(counts) %in% metadata$Sample_ID[metadata$Source %in% c("CBC")])]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))


#====================================================================================
# Batch Normalization and Differential Gene Expression (DESeq2) analysis of Metadata
#====================================================================================

# Batch Normalization ========================================================

#counts_batchnorm <- remove_batch_effect(counts, metadata, ~is_NEBSmallRNA, 1)
counts_batchnorm <- counts
# DESeq2 ======================================================================

# Only count mRNA from exonic regions

counts_batchnorm_mRNA <- filterCountsbyGeneType(counts_batchnorm, annotation_path, "exon", c("mRNA"))

# Generate the DESeq2 Data set for the combined data set (Note: This excludes the MW and KK data bc they dont have
# both Ulcer_duration_cat and IDSA_Score_1to4 metadata)

combined_DESeq2 <- run_DESeq2(metadata = metadata,
                              counts = counts_batchnorm_mRNA,
                              id_colname = "Sample_ID",
                              metadata_vars = c("Ulcer_duration_cat","IDSA_SCORE_1to4", "Lib_prep"),
                              formula = ~ Ulcer_duration_cat + IDSA_SCORE_1to4)

# Test for differential expression between "extreme" values of metadata
combined_DEgenes_UlcerDuration<-DESeq2::results(combined_DESeq2, contrast = c("Ulcer_duration_cat", "2", "0"))
combined_DEgenes_IDSAScore<-DESeq2::results(combined_DESeq2, contrast = c("IDSA_SCORE_1to4", "4","2"))

# Filter for genes with abs(log2FoldChange > 2) and padj < 0.05
combined_DEgenes_UlcerDuration_sig <- filter_DESeq(combined_DEgenes_UlcerDuration, 2, 0.05)
combined_DEgenes_IDSAScore_sig <- filter_DESeq(combined_DEgenes_IDSAScore, 2, 0.05)

#==================================================================================
# Quick Heirarchical Clustering
#==================================================================================

dist_mat <- dist(t(counts_batchnorm_vst), method = 'euclidian')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg) #Note: HH28-P509 correspond to the samples with high bacteria! Highlight in figure


#===================================================================================
# Perform K-means clustering analysis 
#===================================================================================

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

#Fix k-means cluster names 

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

# Run DESeq2
dds_kmeans <- run_DESeq2(metadata = metadata_kmeans, 
                         counts = counts_batchnorm_mRNA,
                         id_colname = "Sample_ID",
                         metadata_vars = c("Source","cluster_res_all"),
                         formula = ~ cluster_res_all) 


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
                       DEseq_kmeans_all_3v1),
  names = c("combined_DEgenes_UlcerDuration", #2vs0
            "combined_DEgenes_IDSAScore", #4vs2
            "DEseq_kmeans_all_3v2",
            "DEseq_kmeans_all_1v2",
            "DEseq_kmeans_all_3v1"))

#=============================================================
# Prepare Validation Data 
#=============================================================

# Normalize for Library Prep Kit 

counts_validation_mRNA <- filterCountsbyGeneType(counts_validation, annotation_path, "exon", c("mRNA"))
counts_validation_mRNA_vst <- DESeq2::vst(as.matrix(counts_validation))

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

# Export DESeq results
write.csv(combined_DEgenes_UlcerDuration_sig,paste0(out_dir,"DESeq2/combined_DEgenes_UlcerDuration.csv" ))
write.csv(combined_DEgenes_IDSAScore_sig,paste0(out_dir,"DESeq2/combined_DEgenes_IDSAScore.csv" ))
write.csv(DEseq_kmeans_all_3v2_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_3v2_sig.csv"))
write.csv(DEseq_kmeans_all_1v2_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_1v2_sig.csv"))
write.csv(DEseq_kmeans_all_3v1_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_all_3v1_sig.csv"))

#Export DESeq summary

write.csv(DESeq_summary, paste0(out_dir,"DESeq2/n_sig_DEgenes.csv"))

# Export kmeans results
write.csv(kmeans_groups_all, paste0(out_dir, "kmeans/kmeans_groups_all.csv"), row.names = FALSE)


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

write.csv(counts_validation_vst, "./data/Example_data/Validation_counts_batchnorm_vst.csv")
write.csv(metadata_validation_mRNA_vst, "./data/Example_data/Validation_metadata.csv", row.names = F, quote = TRUE)


