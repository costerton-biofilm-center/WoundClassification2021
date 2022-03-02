#=====================================================================
# This script performs normalization, PCA, k-means clustering  
# and differential expression analysis included in the manuscript, 
# "Transcriptomic Classification of Lower Extremity Wounds".
#
# Written by: Blaine Gabriel Fritz
#
#=============================================================================

#===================================
#Import libraries
#===================================
library(tidyverse)
library(factoextra)

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

# Run source to import helper functions and fix metadata 
source("./scripts/analysis_utils.R")
source("./scripts/Clean_metadata.R")

# Import counts, metadata and define relative path to annotation file
counts <- read.csv("./data/Example_data/Raw_counts.csv", row.names = 1, check.names = F, comment.char = "") #If importing from a count matrix tsv 
metadata <- read.delim("./data/Example_data/Analysis_metadata.csv", sep = ",")
annotation_path <- c("./data/annotation_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff")
kraken_dir = "./data/kraken"

# Define data Output_directory 
# out_dir needs to contain empty folders called: 
# "normalized_counts", "DESeq2", "kmeans", and "validation"

out_dir <- "./analysis/"

# Format Column types 

metadata <-
  metadata %>%
  mutate(across(c("Gender_0M_1F", "Peripheral_Neuropathy", "PAK", "CKD_stage5", 
                  "IHD", "CCF", "Type_of_Diabetes_T1_T2","Ulcer_duration_cat", 
                  "PEDIS_IDSA_1uninfected_2mild_3mod_4severe"), as.factor)) %>%
  mutate(across(c("Age", "Duration_of_Diabetes_years", "Ulcer_duration_weeks", "WCC_10e9perL", "CRP_mgperL", "ESR_mLpermin",
                  "Neutrophils_10e9perL", "HBA1c", "HBA1c_IFCC"), as.numeric))

#=================================================================================
# Quality filtering
#=================================================================================

# Filtering ===================================================================

# Filter out samples with less than 1M reads and remove from metata
counts <- counts[,colSums(counts) > 1000000]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

# Remove HH5 because it's an extreme outlier - Abnormally high counts for e.g. PADI3.
counts<-counts[,!(colnames(counts) %in% c("HH5"))]
metadata<-subset(metadata, metadata$Sample_ID %in% colnames(counts))

# Organize count names 
counts <- counts[,metadata$Sample_ID]

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

#Extract Names - There probably is a better way to do this
kraken_names <- 
  list.files(kraken_dir, full.names = T) %>%
  set_names(., nm = map(.x = ., ~gsub(".kraken.report", "", basename(.x)))) %>%
  map(function(x) {
    x <- read.delim(x, header = F)
    x$V6 <- trimws(x$V6, which = "left")
    x <- subset(x, V6 %in% c("root", "Bacteria", "Homo sapiens"))
    return(x)}) %>%
  map_df(~as.data.frame(.x), .id = "Sample_ID") %>%
  distinct(Sample_ID)

#Add Sample_ID to kraken data
kraken_data$Sample_ID <- kraken_names$Sample_ID


#Add kraken_results to kmeans metadata 
metadata <- dplyr::left_join(metadata, kraken_data, by="Sample_ID")

#Add Category for high bacteria to kmeans metadata 

metadata$high_bacteria_1High0Low <- ifelse(metadata$Bac_prcnt>10, "1", "0") #1 = High, 2=low


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
  subset(kraken_bac_abundance, Sample_ID %in% metadata$Sample_ID[metadata$Bac_prcnt>10]) 

high_bac_microbiome_top_species <- 
  high_bac_microbiome %>% 
  group_by(Sample_ID) %>%
  arrange(Sample_ID, desc(Rel_abund)) %>%
  filter(Rel_abund>1)


#====================================================================================
# Batch Normalization and Differential Gene Expression (DESeq2) analysis of Metadata
#====================================================================================

# Batch Normalization ========================================================

counts_batchnorm <- remove_batch_effect(counts, metadata, ~Source, 1)

# DESeq2 ======================================================================

# Only count mRNA from exonic regions

counts_batchnorm_mRNA <- filterCountsbyGeneType(counts_batchnorm, annotation_path, "exon", c("mRNA"))

# Generate the DESeq2 Data set for the combined data set 

Clinical_DESeq2 <- run_DESeq2(metadata = metadata,
                              counts = counts_batchnorm_mRNA,
                              id_colname = "Sample_ID",
                              metadata_vars = c("Ulcer_duration_cat",
                                                "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", "high_bacteria_1High0Low"),
                              formula = ~ Ulcer_duration_cat + PEDIS_IDSA_1uninfected_2mild_3mod_4severe + high_bacteria_1High0Low)


DEgenes_UlcerDuration<-DESeq2::results(Clinical_DESeq2, contrast = c("Ulcer_duration_cat", "2", "0"))
DEgenes_IDSAScore<-DESeq2::results(Clinical_DESeq2, contrast = c("PEDIS_IDSA_1uninfected_2mild_3mod_4severe", "4","2"))
DEgenes_high_bacteria_1High0Low <- DESeq2::results(Clinical_DESeq2, contrast = c("high_bacteria_1High0Low", "1","0"))

# Filter for genes with abs(log2FoldChange > 2) and padj < 0.05
DEgenes_UlcerDuration_sig <- filter_DESeq(DEgenes_UlcerDuration, 2, 0.05)
DEgenes_IDSAScore_sig <- filter_DESeq(DEgenes_IDSAScore, 2, 0.05)
DEgenes_high_bacteria_1High0Low_sig <- filter_DESeq(DEgenes_high_bacteria_1High0Low, 2, 0.05)

#==================================================================================
# Variance stabilizing transformation
#===================================================================================

counts_batchnorm_vst <- DESeq2::vst(as.matrix(counts_batchnorm_mRNA))

#==================================================================================
# Quick Heirarchical Clustering
#==================================================================================

dist_mat <- dist(t(counts_batchnorm_vst), method = 'euclidian')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg) #Note: HH28-P509 correspond to the samples with high bacteria! Highlight in figure


#===================================================================================
# Perform K-means clustering analysis 
#===================================================================================

#Remove genes with zero variance 
novar_filter <- apply(counts_batchnorm_vst, 1, sd)
novar_filter <- novar_filter == 0
counts_batchnorm_vst_noVar <- counts_batchnorm_vst[!novar_filter, ]

#Optimal number of clusters?
factoextra::fviz_nbclust(t(counts_batchnorm_vst_noVar), kmeans, method = "silhouette")

#Run Kmeans
set.seed(seed) # Set Seed for reproducibility
kmeans_all <- kmeans(t(counts_batchnorm_vst_noVar), centers = 2, nstart = 25)

# Get the results 
set.seed(seed)
kmeans_groups_all <- data.frame(Sample_ID = names(kmeans_all$cluster), cluster_res_all = as.character(kmeans_all$cluster))
kmeans_groups_all$cluster_res_all <- unlist(kmeans_groups_all$cluster_res_all)

# Switch the cluster names for consistency with previous analysis
kmeans_groups_all$cluster_res_all <- ifelse(kmeans_groups_all$cluster_res_all=="1","2","1")

#=======================================================
# DESEQ2 of k-means Results 
#=======================================================

# Add kmeans results to metadata 
metadata <- dplyr::left_join(metadata, kmeans_groups_all, by="Sample_ID")

# Run DESeq2
dds_kmeans <- run_DESeq2(metadata = metadata, 
                         counts = counts_batchnorm_mRNA,
                         id_colname = "Sample_ID",
                         metadata_vars = c("Source","cluster_res_all"),
                         formula = ~ cluster_res_all) 

# Extract DE genes results
DEseq_kmeans_1v2 <- DESeq2::results(dds_kmeans, contrast = c("cluster_res_all", "1", "2"))

# Extract significant genes
DEseq_kmeans_1v2_sig <- filter_DESeq(DEseq_kmeans_1v2, 2, 0.05)

#=========================================================
# Count # of differentially expressed genes from DESeq 
#=========================================================

DESeq_summary<-
  summarize_DESeq(list(DEgenes_UlcerDuration_sig, #2vs0
                       DEgenes_IDSAScore_sig, #4vs2
                       DEseq_kmeans_1v2_sig),
                  names = c("combined_DEgenes_UlcerDuration", #2vs0
                            "combined_DEgenes_IDSAScore", #4vs2
                            "DEseq_kmeans_1v2_sig"))

#=========================================================
# Test Bacterial Load Between clusters
#=========================================================

fit <- glm(Bac_prcnt/100 ~ cluster_res_all, data = metadata, family = binomial())
anova(fit, test = "Chisq")

t.test(metadata$Bac_prcnt[metadata$cluster_res_all=="2"],
       metadata$Bac_prcnt[metadata$cluster_res_all=="1"])

# Do samples in C2 show decreased alpha diversity? 
n_species_1prcntRelAbund<-
  kraken_bac_abundance %>%
  group_by(Sample_ID) %>% 
  filter(Rel_abund >= 5) %>% #remove noise from low abundance bacteria
  arrange(Sample_ID, desc(Rel_abund)) %>%
  mutate(cum_abund = cumsum(Rel_abund)) %>%
  # filter(cum_abund < 75) %>% #Filter those 
  summarize(n_species = n())

n_species_1prcntRelAbund <- 
  n_species_1prcntRelAbund %>%
    left_join(metadata[,c("Sample_ID", "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", "cluster_res_all")], by = "Sample_ID")

t.test(n_species_1prcntRelAbund$n_species[n_species_1prcntRelAbund$cluster_res_all=="2"],
       n_species_1prcntRelAbund$n_species[n_species_1prcntRelAbund$cluster_res_all=="1"])

t.test(n_species_1prcntRelAbund$n_species[n_species_1prcntRelAbund$PEDIS_IDSA_1uninfected_2mild_3mod_4severe=="4"],
       n_species_1prcntRelAbund$n_species[n_species_1prcntRelAbund$PEDIS_IDSA_1uninfected_2mild_3mod_4severe=="2"])


#Are top species present in C2 also present in other samples?

top_speciesC2 <- subset(high_bac_microbiome_top_species, Sample_ID %in% metadata$Sample_ID[metadata$cluster_res_all=="2"] & Rel_abund>10)

top_species_C2inC1<-subset(kraken_bac_abundance, !c(Sample_ID %in% top_speciesC2$Sample_ID) & 
                             ID %in% top_speciesC2$ID & 
                             Rel_abund > 5) %>%
  group_by(ID) %>%
  summarize("mean_abund" = median(Rel_abund), "sd" = median(Rel_abund), n = n())

#=========================================================
# Run PCA analysis
#=========================================================

source("./scripts/PCA_analysis.R")

#===========================================================================
# Differential Expression of Immune Signals between C1/C2
# 
#===========================================================================

immune_genes <- subset(DEseq_kmeans_1v2, grepl("^CXCL|^CCL[0-9]*$|^NFKB|
                                                   ^IL[0-9]*$|^CD[0-9]*$|^TNF$|
                                                   ^IFN", 
                                                   row.names(DEseq_kmeans_1v2)))
immune_genes <- data.frame(immune_genes)

#=========================================================
# Export Data
#=========================================================

#Export Seed and Session Info 
writeLines(as.character(seed), paste0(out_dir,"seed.txt"))
writeLines(capture.output(sessionInfo()), paste0(out_dir,"sessionInfo.txt"))

# Export Counts
write.csv(counts, paste0(out_dir, "counts/Raw_counts.csv"))
write.csv(counts_batchnorm, paste0(out_dir, "counts/Batchnorm_counts.csv"))
write.csv(counts_batchnorm_mRNA, paste0(out_dir, "counts/Batchnorm_mRNA_counts.csv"))
write.csv(counts_batchnorm_vst, paste0(out_dir,"counts/Batchnorm_mRNA_vst.csv"))

#Export Metadata with results
write.csv(metadata, paste0(out_dir, "other/metadata_with_results.csv"))

#Export DESeq summary
write.csv(DESeq_summary, paste0(out_dir,"DESeq2/n_sig_DEgenes.csv"))

# Export DESeq results
write.csv(DEgenes_UlcerDuration_sig,paste0(out_dir,"DESeq2/combined_DEgenes_UlcerDuration.csv" ))
write.csv(DEgenes_IDSAScore_sig,paste0(out_dir,"DESeq2/combined_DEgenes_IDSAScore.csv" ))
write.csv(DEseq_kmeans_1v2_sig, paste0(out_dir,"DESeq2/DEseq_kmeans_1v2_sig.csv"))
write.csv(DEgenes_high_bacteria_1High0Low_sig, paste0(out_dir, "DESeq2/DEseq_HighvsLowBac_sig.csv"))
write.csv(immune_genes, paste0(out_dir,"DESeq2/DEseq_kmeans_1v2_Immune_all.csv" ))

# Export kmeans results
write.csv(kmeans_groups_all, paste0(out_dir, "kmeans/kmeans_groups_all.csv"), row.names = FALSE)

#Export Data for GO_analysis

write.table(row.names(counts_batchnorm_vst),
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/counts_batchnorm_genelist.txt"))

write.table(row.names(subset(DEseq_kmeans_1v2_sig,log2FoldChange > 2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_1v2_sigUP.txt"))

write.table(row.names(subset(DEseq_kmeans_1v2_sig,log2FoldChange < -2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_kmeans_1v2_sigDOWN.txt"))

write.table(row.names(subset(DEgenes_high_bacteria_1High0Low_sig,log2FoldChange > 2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_HighvsLowBac_sigUP.txt"))

write.table(row.names(subset(DEgenes_high_bacteria_1High0Low_sig,log2FoldChange < -2)), 
            row.names = F, sep = "", quote = F, col.names  = F,
            paste0(out_dir,"GO_analysis/DEseq_HighvsLowBac_sigDOWN.txt"))

