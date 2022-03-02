
####################################################################################
##########  Analysis Functions  ####################################################
####################################################################################

get_metadata <- function(metadata_path){
  #Metadata should be a tsv file where the first three columns are:
  #Sample_ID, Read1_name, and Read2_name...Colnames don't have to be exactly this
  #The rest of the columns should be the metadata information
  
  metadata <- lapply(metadata_path, read.delim, stringsAsFactors=F)
  
  if (length(metadata_path)>1){
    require(dplyr)
    metadata <- dplyr::bind_rows(metadata)
  }
  else {
    metadata <- data.frame(metadata)
    return(metadata)
  }
  if (any(is.na(metadata[,1]))){
    stop("NAs in sample_id column. Make sure that the name first column is identical in each metadatafile.")
  }
  
  return(metadata)
}


get_counts <- function(counts_path, min_gene_length=0){#Could also be a vector of folder paths
  file_list <- lapply(counts_path,list.files, pattern = "*.txt", full.names = T)
  file_list <- unlist(file_list)
  counts <- lapply(file_list, read.delim, skip = 1, check.names = F)

  #if using featureCounts, length is not always accurate 
  counts_lengthFiltered <- lapply(counts, function(x) x[x$Length>min_gene_length,] )
  
  #Isolate only the Geneid and the counts
  cleaned<-lapply(counts_lengthFiltered,function(x){
    x_clean<-x[,c(1,ncol(x))]
  })

  gene_ids<-unlist(cleaned[[1]][1]) #Get ids from first list entry
  df_out<-data.frame(c(1:length(gene_ids)), check.names = F)
  rownames(df_out)<-gene_ids
  
  for (sample in 1:length(cleaned)){
    sample_name<-colnames(cleaned[[sample]][2])
    df_out[,sample]<-unlist(cleaned[[sample]][2])
    colnames(df_out)[sample]<-sample_name
  }
  #Clean up the column names....Formatted like "A.BUNCH.OF.JUNK.sample_id.aligned.bam"
  #We only want the sample_id
  
  require(stringr)
  file_names <- basename(colnames(df_out))
  sample_ids <- gsub(pattern = ".aligned.bam", "", file_names)
  colnames(df_out) <- sample_ids
  
  return(df_out)
}

plot_PCA <- function(counts, metadata, PCs, color_variable=NULL, normalize=T){
  require(dplyr)
  require(tibble)
  require(factoextra)
  require(ggplot2)
  
  if(normalize){
    counts<-DESeq2::vst(as.matrix(counts), blind = T, nsub = 1000)
  }
  
  PCA <- prcomp(t(counts), scale. = F)
  data <- data.frame(PCA$x)
  var_explained <- PCA$sdev^2/sum(PCA$sdev^2)
  metadata_id_column_name <- colnames(metadata)[1]
  

  #tidy the data and add the metadata
  plot_data <- 
    data %>%
    rownames_to_column(var=metadata_id_column_name) %>%
    left_join(metadata) 
  
  #scale_color<- scale_color_gradient()
  scale_color <- scale_color_distiller(palette = "Oranges", direction = 1)

  if(is.null(color_variable)){
    plot <- 
      ggplot(data = plot_data, aes_string(x = plot_data[, PCs[1]+1], y= plot_data[, PCs[2]+1]))+
      labs(x=paste0("PC", PCs[1], ": ", round(var_explained[PCs[1]]*100,1),"%"),
           y=paste0("PC", PCs[2], ": ", round(var_explained[PCs[2]]*100,1),"%"))+
      geom_point(aes_string(color=color_variable), na.rm = F)+
      scale_color_gradient(na.value = "grey")+
      scale_x_reverse()
      
      return(plot)
  }
  else{
    #scale_color<- scale_color_gradient()
    scale_color <- scale_color_distiller(palette = "Oranges", direction = 1, na.value = "grey")
    
    if(typeof(plot_data[,color_variable])=="character"){
      scale_color<-scale_colour_brewer(palette="Dark2", na.value = "grey")
      }
  #Make the plot 
  plot <- 
    ggplot(data = plot_data, aes_string(x = plot_data[, PCs[1]+1], y= plot_data[, PCs[2]+1]))+
    labs(x=paste0("PC", PCs[1], ": ", round(var_explained[PCs[1]]*100,1),"%"),
         y=paste0("PC", PCs[2], ": ", round(var_explained[PCs[2]]*100,1),"%"))+
    geom_point(aes_string(color=color_variable), na.rm = F)+
    scale_color+
    theme_classic()+
    scale_x_reverse()
  
  return(plot) 
  }
}

#Clustering needs normalzed counts as input
run_clustering <- function(norm_counts, metadata, n_contrib, PCs, n_clust, plot=F, color_var=NULL){
  require(factoextra)
  require(gridExtra)
  #make some plots
  PCA <- prcomp(t(norm_counts), scale. = F)
  PCA_plot <- plot_PCA(norm_counts, metadata, PCs, color_variable = color_var, normalize=F)
  
  scree_plot <- fviz_screeplot(PCA, addlabels = TRUE, ylim = c(0, 50))
  contrib_1<-fviz_contrib(PCA, choice = "var", axes = PCs[1], top = n_contrib)
  contrib_2<-fviz_contrib(PCA, choice = "var", axes = PCs[2], top = n_contrib)
  
  #Determine the optimal number of clusters
  data_clustering <- PCA$x[,PCs]
  kmeans_res <- kmeans(data_clustering, n_clust, nstart = 25)
  
  if (plot==F){
    return(kmeans_res$cluster)
  }
  
  cluster_plot <- fviz_cluster(kmeans_res, data = data_clustering,
               # = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"),
               ggtheme = theme_minimal(),
               main = "Partitioning Clustering Plot")
  optimal <- fviz_nbclust(data_clustering, kmeans, method = "silhouette")
  grid.arrange(PCA_plot, scree_plot, contrib_1, contrib_2, cluster_plot, optimal)
}


add_clusters_to_metadata <- function(metadata, cluster_res_vector){
  cluster_res <- data.frame(cluster_res_vector)
  
  if (any(row.names(cluster_res)!=metadata[,1])){
    stop("ERROR: Order of cluster results doesn't match order of samples in metadata")
  }
  
  metadata$cluster_res<-as.factor(cluster_res_vector)
  return(metadata)
  
}

remove_batch_effect<-function(counts, metadata, contrast_formula, metadata_id_col) {
  #Trim unused rows from metadata 
  sample_names <- colnames(counts)
  metadata <- metadata[metadata[,1] %in% sample_names,]
  
  
  #Reorder count cols and metadata rows so they match and test
  require(DESeq2)
  row.names(metadata)<-metadata[,metadata_id_col]
  counts_ordered<-counts[,metadata[,metadata_id_col]]
  all(row.names(metadata)==colnames(counts_ordered))
  
  dds <- DESeqDataSetFromMatrix(countData = counts_ordered,
                                colData = metadata,
                                design = contrast_formula)
  dds <- DESeq(dds)
  
  results <- results(dds)
  results_significant<-subset(results,padj<.05)
  print(paste("Filtering ", nrow(results_significant), " genes due to batch effects."))
  batch_genes <- row.names(counts) %in% row.names(results_significant)
  counts_batchnorm<-data.frame(counts[!batch_genes,])
  
  
  return(counts_batchnorm)
}

kraken_bac_human<- function(kraken_dirs){
  files <- list.files(kraken_dirs, full.names=T)
  sample_id <- gsub(".kraken.report", "", basename(files))
  data <- lapply(files, read.delim, header=FALSE)
  names(data)<-sample_id
  
  #Extract Bac and Human Info
  df.bac_human <- lapply(data,function(kraken_report){
    kraken_report$V6<-trimws(kraken_report$V6) #whitespace before and after names
    bac<-subset(kraken_report, V4=="D"&V6=="Bacteria")
    human<-subset(kraken_report, V4=="S"&V6=="Homo sapiens")
    
    df.bac_human <- rbind(human,bac)
    
    return(df.bac_human)
  })
  
  require(dplyr)
  summarized_data <- bind_rows(df.bac_human, .id="Sample_ID")
  
  colnames(summarized_data)<-c("Sample_ID", "prcnt_classified", 
                               "n_subclades","n_clade",
                               "Rank", "tax_id", "name")
  return(summarized_data)
}


generate_summary <- function(counts_df, metadata_df, kraken_df){
  #Returns a data frame with sample_id, total_reads, n_bac_reads, n_human_reads,
  #percent classified_bac, percent_classified_human, ?topbacterialspecies,?%ofbacreads_to_top_species
  total_reads <- colSums(counts_df)
  total_reads <- total_reads[order(names(total_reads))]
  sample_names<-names(total_reads)
  
  bacteria <- subset(kraken_df, name=="Bacteria")
  bacteria <- bacteria[order(bacteria[,1]),c(1,2,3)]
  colnames(bacteria)<-paste("bac",colnames(bacteria))

  human <- subset(kraken_df, name=="Homo sapiens")
  human <- human[order(human[,1]),c(1,2,3)]
  colnames(human)<-paste("human",colnames(human))

  #check that the order is right
  if(all(c(sample_names==human[,1], sample_names==bacteria[,1]))){
    output <- data.frame("Human_reads_exon" = total_reads,
                         bacteria[,2:3], 
                         human[,2:3])
    return(output)
  }

  else{stop("ERROR: Sample Names are not consistent between counts, kraken, and metadata")}
}

plot_PCA_corr <- function(PCA_obj, metadata, PC, corr_variable) {
  metadata_info<-data.frame(metadata[,1],metadata[,corr_variable])
  PCA_info<-data.frame(row.names(PCA_obj$x), PCA_obj$x[,PC])
  
  data<-data.frame(metadata_info[order(metadata_info[,1]),2], PCA_info[order(PCA_info[,1]),2])
  row.names(data)<- PCA_info[order(PCA_info[,1]),1]
  
  require(ggplot2)
  
  ggplot(data=data, aes_string(x=colnames(data)[2], y=colnames(data)[1]))+
    geom_point()
}


run_DESeq2 <- function(metadata, counts, id_colname, formula, metadata_vars){
  row.names(metadata)<-metadata[,id_colname]
  
  #Remove NAs from Metadata
  metadata <- metadata[,metadata_vars, drop = F]
  metadata <- metadata[complete.cases(metadata), 1:ncol(metadata), drop = F]
  
  #Remove counts if not present in metadata
  if(any(FALSE == (colnames(counts) %in% row.names(metadata)))){
    counts<-counts[, colnames(counts) %in% row.names(metadata)]
  }
  
  #Order Sample names in counts and metadata
  metadata <- metadata[order(row.names(metadata)),]
  counts <- counts[,order(colnames(counts))]
  
  
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, metadata, formula)
  dds <- DESeq2::DESeq(dds)
  
  return(dds)

}


filter_DESeq <- function(DEseq_results, LFC, padj_val){
  data_filtered <- subset(DEseq_results, abs(log2FoldChange) > LFC & padj < padj_val)
  return(data_filtered)
}



filterCountsbyGeneType <- function(counts, annotation_path, feature_types, transcript_types){
  #Read in annotations
  annotation <- read.delim(annotation_path,header=FALSE, comment.char="#")
  annotation_featurefilter <- subset(annotation, V3 %in% feature_types)
 
   #Filter V9 to extract the gene_type (mRNA, ncRNA, etc)
  require(stringr)
  annotation_featurefilter$gbkey <- stringr::str_extract(annotation_featurefilter$V9, "gbkey=[A-z,_]*")
  annotation_featurefilter$gbkey <- gsub("gbkey=", "", annotation_featurefilter$gbkey)
  annotation_featurefilter$geneid <- stringr::str_extract(annotation_featurefilter$V9, "gene=.*")
  annotation_featurefilter$geneid <- gsub("gene=", "", annotation_featurefilter$geneid)
  annotation_featurefilter$geneid <- gsub(";.*$", "", annotation_featurefilter$geneid)
  
  #Cleaned up info data frame  
  annotation_featurefilter <- annotation_featurefilter[,c("geneid", "gbkey", "V3", "V1")]
  
  #Subset by transcript_types 
  annotation_featurefilter <- subset(annotation_featurefilter, gbkey %in% transcript_types)
  
  #Get a list of the subsetted genes
  selected_genes <- unique(annotation_featurefilter$geneid)
  
  counts_subset <- subset(counts, row.names(counts) %in% selected_genes)
  
  return(counts_subset)
  
}

get_gene_names <- function(gene_list, annotation_path, feature_types, transcript_types){
  #Read in annotations
  #Takes a vector if gene names, gff file path, transcript/feature types 
  #and returns a df with the input list and the full name from the annotation file
  require(stringr)
  require(dplyr)
  #Read in annotations
  annotation <- read.delim(annotation_path,header=FALSE, comment.char="#")
  annotation_featurefilter <- subset(annotation, V3 %in% feature_types)
  
  #Extract the important data 
  data <- data.frame("geneid" = stringr::str_extract(annotation_featurefilter$V9, "gene=.*"),
                     "product" = stringr::str_extract(annotation_featurefilter$V9, "product=.*"),
                     "gbkey" =  stringr::str_extract(annotation_featurefilter$V9, "gbkey=[A-z,_]*"))
  
  data$gbkey <- gsub("gbkey=", "", data$gbkey)
  data$geneid <- gsub("gene=", "", data$geneid)
  data$geneid <- gsub(";.*$", "", data$geneid)
  data$product <- gsub("product=", "", data$product)
  data$product <- gsub(";.*$", "", data$product)
  
  #delete info about transcript variants
  data$product <- gsub("%2C.*$", "", data$product)

  #delete info about transcript variants
  data$product <- gsub("%2C.*$", "", data$product)
  
  #Subset by transcript_types 
  data <- subset(data, gbkey %in% transcript_types)

  #Subset by transcript_types 
  data <- subset(data, gbkey %in% transcript_types)

  data<-
  data %>%
    group_by(geneid, product) %>%
    summarize()
  
  data_out <- data[,c("geneid", "product")]
  data_out <- data_out[data_out$geneid %in% gene_list,]

  colnames(data_out)<-c("Gene", "Product")
  return(data_out) 
}


summarize_DESeq<-
  function(data_list, names){
    data<-
      lapply(data_list, function(x){
        data <-data.frame(x)
        data_filter<-subset(data, padj<0.05)
        n_up <- nrow(data_filter[data_filter$log2FoldChange>2,])
        n_down <- nrow(data_filter[data_filter$log2FoldChange<(-2),])  
        return(list("n_up" = n_up, "n_down" = n_down))
      })
    
    names(data)<-names
    return(data.frame(n_sig_genes = unlist (data)))
  }

# Adapted from from https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot, scale = "discrete", pointSize = 0.5, 
                           textSize = 3, spaceLegend = 0.1) {
  if(scale == "continuous"){
    myPlot +
      guides(shape = guide_colorbar(override.aes = list(size = pointSize)),
             color = guide_colorbar(override.aes = list(size = pointSize)))+
      theme(legend.title = element_text(size = textSize), 
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  } else {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize)))+
               theme(legend.title = element_text(size = textSize), 
                     legend.text  = element_text(size = textSize),
                     legend.key.size = unit(spaceLegend, "lines"))
             }
}

#Remove redundant terms from panther results
reducePanther<- 
  function(panther_csv, revigo_csv){
    panther <- read.delim(panther_csv, skip = 12, stringsAsFactors = F)
    revigo <- read.delim(revigo_csv, sep = ",")
    #GO term needs to be extracted from a string in panther output
    panther$TermID <- str_extract(panther[,1], "GO:[0-9]*")                         
    
    out <- left_join(panther, revigo, by = "TermID")
    #for some reason a space before word in Eliminated 
    out <- subset(out, Eliminated == " False")
    return(out)
  }

#Calculate accuracy

calc_Accuracy<-
  function(data, guess_col, prediction_col){
    if(length(data[,guess_col])!=length(data[,prediction_col])){
      stop("ERROR: Length of guess and prediction not the same")
    }
    data<-data[complete.cases(data),] #Removes any NA values
    n_correct <- length(which(data[,guess_col] == data[,prediction_col]))
    n_incorrect <- length(which(data[,guess_col]!= data[,prediction_col]))
    if((n_correct+n_incorrect)!= length(data[,guess_col])){
      stop("ERROR: Length of n_correct and n_incorrect do not \n 
           sum to input length. Check input data for NAs.")
    }
    
    n_summary <- list("n_correct" = n_correct, "n_incorrect" = n_incorrect, 
                      "n_total" = n_correct + n_incorrect)
    print(n_summary)
    
    accuracy <- n_correct/(n_correct+n_incorrect)
    return(accuracy)
  }

fitmodel<-
  function(gene_counts, IDSA_scores, high_bacteria){
    require(car)
    #Get counts
    gene_counts <- as.numeric(gene_counts)
    IDSA_scores <- factor(IDSA_scores, levels = c("2","3","4"))
    high_bac <- factor(high_bacteria, levels = c("low", "high"))

    model_data <- data.frame(IDSA_scores, high_bac, gene_counts)
    
    #Remove NA
    model_data <- model_data[complete.cases(model_data),]
    
    #fit model
    model <- car::Anova(lm(gene_counts ~ IDSA_scores + high_bac, data = model_data), type = 2)
    
    #Extract parameters
    df_IDSA<-model$Df[1]
    SS_IDSA <- model$`Sum Sq`[1]
    var_IDSA <- SS_IDSA/df_IDSA
    
    df_highbac <- model$Df[2]
    SS_highbac <- model$`Sum Sq`[2]
    var_highbac <- SS_highbac/df_highbac
  

    df_resid <- model$Df[3]
    SS_resid <- model$`Sum Sq`[3]
    var_resid <- SS_resid/df_resid
    
    #calculate eta squared for each parameter
    
    contrib_highbac <- SS_highbac/sum(model$`Sum Sq`)
    contrib_idsa <- SS_IDSA/sum(model$`Sum Sq`)

    output <- list("var_prcnt_highbac" = contrib_highbac, 
                   "var_prcnt_idsa" = contrib_idsa)

    return(output)
  }

getSharedGO<-function(df1, df2){
  data <-
    lapply(list(df1, df2), function(x){ #Make colnames consistent
      data <- x
      colnames(data) #Make colnames consistent
      colnames(data)[2:7]<-c("n_ref_genes", "n_DE_genes", "Expected", 
                             "over.under", "fold.Enrichment", "P.value")
      return(data)
    })
  #Combine the two dfs
  data_combined <- bind_rows(data)
  data_combined <- 
    data_combined[data_combined$TermID %in% intersect(df1$TermID, df2$TermID),] %>%
    group_by(TermID) %>%
    filter(P.value == max(P.value))
  return(data_combined)
}



#################################################
### Plotting Functions ##########################
#################################################

make_Figure1<-function(counts, metadata, variables){
  
  print(head(colnames(metadata)))
  print(variables)
  
  print(head(metadata[,colnames(metadata) %in% variables]))
  
  #metadata[,colnames(metadata) %in% variables]<-lapply(metadata[,colnames(metadata) %in% variables],as.factor)
  
  #Set up the data
  malone_data<-counts[,c(metadata$Sample_ID[metadata$Source=="CBC"])]
  plots_malone<-lapply(variables, plot_PCA, counts=malone_data, metadata=metadata, PCs=c(1,2), normalize=T)
  plots_alldata<-lapply(variables, plot_PCA, counts=counts, metadata=metadata, PCs=c(1,2), normalize=T)

  
  #Multiple PC2 of malone_data by -1 to match alldata
  plots_malone<-lapply(plots_malone, function(plot){
    plot<- plot + scale_y_reverse() 
    return(plot)
  })

  
  require(cowplot)
  #labels<-c("A","E","B","F","C","G","D","H",)
  plot_grid(plotlist=c(plots_malone, plots_alldata), labels="AUTO", ncol=2, byrow = F, align="hv")
  
}

plot_contribs<-function(norm_counts, PCs, n_contrib){
  PCA <- prcomp(t(norm_counts), scale. = F)
  contribs1 <- factoextra::fviz_contrib(PCA, choice = "var", axes = 1, top = n_contrib)+
    ggtitle("Contributions of variables to PC1")
  contribs2 <- factoextra::fviz_contrib(PCA, choice = "var", axes = 2, top = n_contrib)+
    ggtitle("Contributions of variables to PC2")
  

  return(list(contribs1, contribs2))
}


