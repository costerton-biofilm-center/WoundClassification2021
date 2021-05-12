#================================================================
# This script will generate the figures included 
# in the manuscript ""Transcriptomic Classification of Lower 
# Extremity Wounds"
#
# Written by: Blaine Fritz
#
#================================================================

# The script uses the output of FINAL_DataAnalysis.R, so be sure that the script
# has been run or sourced without error. 
# You can source the script by running:
# source("./scripts/FINAL_DataAnalysis.R")

# The source directory for the R session should be the project folder
# NOT the directory containing the scripts

#===============================================================
# Import Required Libraries
#===============================================================
library(ggplot2)
library(cowplot)
library(gridExtra)

#===============================================================
# Figure 1 
#===============================================================
#NOTE: This also generates and saves PCA plots of all the
#      metadata variables in the metadata file.  


#Make and Export PCAs for all the Metadata Variables 
PCA_plots_CBC<-
  lapply(colnames(metadata_kmeans[,5:ncol(metadata_kmeans)]), function(metadata_cat) {
    print(metadata_cat)
    if(metadata_cat %in% c("Sample_ID", "R1_name", "R2_name", "Specific_ID")){} 
    else{
      
      plot<-
        plot_PCA(counts_cbconly_vst, metadata_kmeans, c(1,2), 
                 metadata_cat, normalize = F)
      # ggsave(paste0("./analysis/Figures/PCA_CBC_", metadata_cat, ".png"),
      #        plot = plot,
      #        width = 5, height = 5)
      return(plot+theme(legend.position = "none"))
    }
  })

PCA_plots_all<-
  lapply(colnames(metadata_kmeans[,5:ncol(metadata_kmeans)]), function(metadata_cat) {
    if(metadata_cat %in% c("Sample_ID", "R1_name", "R2_name", "Specific_ID")){} 
    else{
      plot<-
        plot_PCA(counts_batchnorm_vst, metadata_kmeans, c(1,2), 
                 metadata_cat, normalize = F)
      # ggsave(paste0("./analysis/Figures/PCA_All_", metadata_cat, ".png"),
      #        width = 5, height = 5)
      return(plot+theme(legend.position = "none"))
    }
  })

#Assign names to the list outputs for each access
names(PCA_plots_CBC)<-colnames(metadata_kmeans[,5:ncol(metadata_kmeans)])
names(PCA_plots_all)<-colnames(metadata_kmeans[,5:ncol(metadata_kmeans)])


#Make the labels for the PCA plot section

titles <- lapply(c("a. Design",
                   "b. Pre-normalization",
                   "c. Batch Normalized",
                   "d. Infection Severity (IDSA/PEDIS)",
                   "e. Ulcer Duration", 
                   "f. Proportion Bacterial Reads to Human"),
                 function(title){
                   ggdraw()+
                     draw_label(title, 
                          x=0, hjust = 0)+
                     theme(
                       # add margin on the left of the drawing canvas,
                       # so title is aligned with left edge of first plot
                       plot.margin = margin(0, 0, 0, 7)
                     )
                 })

#Get the flowchart and add the title 

flowchart<-
  ggdraw() + 
  draw_image("./analysis/Figures/Sampling_overview.png")+
  theme(plot.margin = margin(0,0,0,0,"cm"))

flowchart <- plot_grid(titles[[1]], flowchart, 
                       ncol = 1, 
                       rel_heights = c(0.2,0.7)
)

#Make PCAs for pre and post normalization 
pre_normalization_pca <- 
  plot_PCA(counts_filtered, metadata_kmeans, c(1,2), 
           "Lib_prep", normalize = T)+
  theme(legend.position = "none")

post_normalization_pca  <- plot_PCA(counts_batchnorm_vst, metadata_kmeans, c(1,2), 
                                    "Lib_prep", normalize = F)+
  theme(legend.position = "none")


# Manually prettying up the PCAs in Figure1 

PCA_plots_all[["Bac_prcnt"]]<- addSmallLegend(PCA_plots_all[["Bac_prcnt"]], pointSize = 2, 
                                            textSize = 5, spaceLegend = 0.4,
                                            scale = "continuous")+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "% Bacterial Reads")

PCA_plots_all[["IDSA_SCORE_1to4"]]<- addSmallLegend(PCA_plots_all[["IDSA_SCORE_1to4"]], pointSize = 2, 
                                              textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "IDSA Score")+
  scale_color_manual(values = c("#84D7E1B2", "#FF6F00B2", "#C71000B2"), na.value = "#3F4041B2")

PCA_plots_all[["Ulcer_duration_cat"]]<- addSmallLegend(PCA_plots_all[["Ulcer_duration_cat"]], pointSize = 2, 
                                                    textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "Ulcer Duration")+
  scale_color_manual(values = c("#ADE2D0B2","#FF95A8B2", "#8A4198B2"), na.value = "#999999")

pre_normalization_pca<-addSmallLegend(pre_normalization_pca, pointSize = 2, 
                                      textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  scale_color_manual(values = c("#3C5488FF","#E64B35FF", "#91D1C2FF"), na.value = "#999999")+
  labs(color = "Library Prep Kit")

post_normalization_pca<-addSmallLegend(post_normalization_pca, pointSize = 2, 
                                      textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  scale_color_manual(values = c("#3C5488FF","#E64B35FF", "#91D1C2FF"), na.value = "#999999")+
  labs(color = "Library Prep Kit")


#Make the plot_grids
pca_plotgrid <- 
  plot_grid(titles[[4]], plot_grid(plotlist = c(PCA_plots_all["IDSA_SCORE_1to4"])),
            titles[[5]], plot_grid(plotlist = c(PCA_plots_all["Ulcer_duration_cat"])),
            titles[[6]], plot_grid(plotlist = c(PCA_plots_all["Bac_prcnt"])),
            rel_heights = (c(0.2,0.8,0.2,0.8,0.2,0.8)),
            ncol = 1,
            align = "hv"
  )

normalization_plotgrid <- 
  plot_grid(titles[[2]], pre_normalization_pca, titles[[3]], post_normalization_pca,NULL,
            ncol = 1,
            nrow = 5,
            byrow = T, 
            rel_heights = c(0.1,0.4, 0.1, 0.4, 0.1),
            vjust = -0.1)

combined<-plot_grid(normalization_plotgrid,NULL, 
                    pca_plotgrid,NULL, ncol = 4, 
                    rel_widths = c(0.5,0.05, 0.5, 0.05),
                    align = "h")

#Make the final plot 

final_plot <- 
  plot_grid(flowchart, combined,
            nrow = 2, rel_heights = c(0.35, 1))

#ggsave("./analysis/Figures/Figure1.png",
ggsave("./analysis/Figures/Figure1.png",      
       final_plot, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)

#Also save as svg for editing 
ggsave("./analysis/Figures/Figure1.pdf",      
       final_plot, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)


#=============================================#
# Figure 2 Kmeans, PCA contibs, GO            #
#=============================================#

# Go Terms - Output frmo Panther, Reduced by Revigo 
GO_3v2UP <- reducePanther(panther = "./analysis/GO_analysis/RESULTS_3v2sigUP_12032021.txt", 
                          revigo = "./analysis/GO_analysis/Revigo_3v2_OUT.csv")
GO_1v2UP <- reducePanther(panther = "./analysis/GO_analysis/RESULTS_1v2sigUP_12032021.txt", 
                          revigo = "./analysis/GO_analysis/Revigo_1v2_OUT.csv")
GO_3v1UP <- reducePanther(panther = "./analysis/GO_analysis/RESULTS_3v1sigUP_12032021.txt", 
                          revigo = "./analysis/GO_analysis/Revigo_3v1_OUT.csv")

#theme 
mytheme <- gridExtra::ttheme_default(core=list(fg_params=list(hjust = 0, x = 0, fontsize = 9)),
                                     colhead=list(fg_params=list(hjust = 0, x = 0, fontsize = 9))
                                     )

plot_data <- lapply(list(GO_3v2UP, GO_1v2UP), function(x){
  data <- x[,c(1,6,7)] #Select name, log fold change, padj
  colnames(data) <- c("GO Term", "Fold Enric.", "adj. p-value")
  newlines <- lapply(data[,1], function(x){
    if(nchar(x)>32){
      x <- gsub('(.{1,32})(\\s)', '\\1\n', x)
    }
    else{
      x <- gsub(pattern = " *\\(GO", replacement = "\n(GO", x)
    }  
  })
  data[,1]<-unlist(newlines)
  #data[,1]<-gsub(pattern = " *\\(GO", replacement = "\n(GO", data[,1]) #Trim GO Term IDs from the names 
  data <- data[order(data[,3], decreasing = F),] #sort by fold enrichment
  plot_out <- tableGrob(head(data), rows = NULL, theme = mytheme)
  plot_out$widths <- unit(c(47,20,20), "mm") 
  return(plot_out)
})

#Get the contributions to the PCA 

PCA_maincontribs<-
  lapply(list(counts_batchnorm_vst, counts_cbconly_vst), plot_contribs, PCs = c(1,2), n_contrib = 20)


plot_grid(PCA_maincontribs[[2]])

#Fix up PCA plot 
PCA_plots_all[[15]]<-
addSmallLegend(PCA_plots_all[[15]], pointSize = 2, 
               textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "K-means Cluster")


#Make the names for the plot sections

titles <- lapply(c("a. Contributions of Variables to PCs",
                   "b. PCA of K-means Results",
                   "c. Enriched Pathways in Cluster 2",
                   "d. Enriched Pathways in Cluster 3"),
                 function(title){
                   ggdraw()+
                     draw_label(title, 
                                x=0, hjust = 0)+
                     theme(
                       # add margin on the left of the drawing canvas,
                       # so title is aligned with left edge of first plot
                       plot.margin = margin(0, 0, 0, 7)
                     )
                 })



#Build the plot 


GO_grid <- plot_grid(plot_data[[1]], plot_data[[2]], nrow = 2, align = "hv", labels = c("c. Enriched Pathways (C3 vs C2)", 
                                                                                        "d. Enriched Pathways (C1 vs C2)"),
                     label_x = -0.4)
PCA_grid <- plot_grid(PCA_maincontribs[[2]], PCA_plots_all[[15]], nrow = 2, rel_heights = c(0.6,0.3),
                      labels = c("a.", "b."), label_y = c(1, 1.1))

FINAL_Fig2<-
plot_grid(PCA_grid, GO_grid,  ncol = 2)


ggsave("./analysis/Figures/Figure2.png",      
       FINAL_Fig2, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)

ggsave("./analysis/Figures/Figure2.pdf",      
       FINAL_Fig2, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)


#==========================================================#
# Figure  3 - Predictor Genes Table and Expression Plots   #
#==========================================================#

### Import results from SVM ################################

SVM_dir <- "./analysis/linearSVM_out/"
SVM_names <- "GoodGenes_*"

# Get the good genes from the classifier
SVM_files<-list.files(SVM_dir, pattern = "GoodGenes.*.csv$", full.names = T)
SVM_genes <- lapply(SVM_files, read.delim)

# Fix names to get metadata cat 
names(SVM_genes) <- list.files(SVM_dir, pattern = "GoodGenes.*.csv$")
names(SVM_genes) <- gsub("GoodGenes_","", names(SVM_genes))
names(SVM_genes) <- gsub("*.csv","", names(SVM_genes))

# Hard-coding the metadata names in so they look nice in the plots. Prints a warnnig 
warning(paste0("Changing: ", names(SVM_genes), " to ", c("K-means Cluster", "IDSA/PEDIS Score", "Ulcer Duration"),
               " in SVM_names_new. Make sure names are correct.", sep = "\n"))

#This will be used for the graphs to have nice labels
SVM_names_new <- c("K-means Cluster", "IDSA/PEDIS Score", "Ulcer Duration")


### Figure 3a #############################################

#Get the annotations for the good genes (takes several minutes to run!)

gene_products<-
lapply(SVM_genes, function(x){
  data <- unlist(x)
  get_gene_names(data, annotation_path, "exon", "mRNA")
})


### Figure 3b ###############################################

# SVM_genes <- 
#   lapply(list.files(SVM_dir, pattern = SVM_names , full.names = T), read.delim, sep = ",")


# For each metadata cat of interest, make a data set for plotting
plot_data <- lapply(seq_along(SVM_genes), function(x){
  expression_values <- subset(counts_batchnorm_vst, row.names(counts_batchnorm_vst) %in% unlist(SVM_genes[x]))
  expression_values <- as.data.frame(expression_values)
  expression_values <-
    expression_values %>%
      rownames_to_column("Gene_ID") %>%
      pivot_longer(cols = -1, names_to = "Sample_ID", values_to = "Norm_expression") %>%
      left_join(metadata_kmeans, by = "Sample_ID")

})

#Make a faceted grid of violin plots for each metadata cat in SVM output
plots <- lapply(seq_along(plot_data), function(x){
  violin_plot<-
  ggplot(plot_data[[x]], aes_string(x=names(SVM_genes)[x], y="Norm_expression"))+
    geom_violin(scale = "width", trim = FALSE, aes_string(fill=names(SVM_genes)[x]))+
    geom_jitter(position = position_jitter(width=0.05, height = 0))+
    xlab(SVM_names_new[x])+
    ylab("Normalized Expression")+
    theme_half_open()+
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
    facet_wrap(~ Gene_ID, nrow = 5, ncol = 4)
  return(violin_plot)
})


### Combine 3a and 3b for Final Figure

lapply(seq_along(plots), function(x){
  name <- names(SVM_genes)[x]
  
  #Get accuracy plot (Exported by ML python script)
  accuracy <- ggdraw() + 
    draw_image(paste0("./analysis/Figures/Accuracy_", name,".png"))+
    theme(plot.margin = margin(0,0,0,0,"cm"))
  
  plot_leftside <- plot_grid(tableGrob(gene_products[[x]], rows = NULL),
                             accuracy, nrow = 2, ncol = 1, labels = c("", "b"),
                             rel_heights = c(0.5,0.4))
  
  
  out_plot <- plot_grid(plot_leftside,
                    plots[[x]], ncol = 2,
                    labels = c("a", "c"), rel_widths = c(0.5,0.5))
  
  
  ggsave(paste0("./analysis/Figures/Figure3_", name, ".png"),
         out_plot, units = "mm",
         width = 300,
         height = 300,
         dpi = 300) 
  
  ggsave(paste0("./analysis/Figures/Figure3_", name, ".pdf"),
         out_plot, units = "mm",
         width = 300,
         height = 300,
         dpi = 300) 
})


######################################
# Figure 4 - SVC Results Figure     ##
######################################

SVC_res<-
lapply(list.files("./analysis/linearSVM_out/", pattern = "COEF*", full.names = T),
       read.csv)

names <- basename(list.files("./analysis/linearSVM_out/"))
names <- gsub(".csv", "", names)
names

SVC_plots<-lapply(SVC_res, function(plot_data){
  require(ggplot2)
  
  ggplot(plot_data, aes(x = Coefficient, y = reorder(Genes, Coefficient)))+
    geom_bar(stat="identity", fill = "#6497B1")+
    ylab("")+
    theme(plot.margin = unit(c(1,0,1,0), "cm"))
  
})


cluster_plots <- plot_grid(plotlist = SVC_plots[c(2,3,1)], 
                           ncol = 3,
                           labels = c("a. Cluster 1", "b. Cluster 2", "c. Cluster 3"),
                           align = "hv")

cluster_plots

ggsave("./analysis/Figures/Figure4.png",      
       cluster_plots, units = "mm", 
       width = 300, 
       height = 150, 
       dpi = 300)

ggsave("./analysis/Figures/Figure4.pdf",      
       cluster_plots, units = "mm", 
       width = 300, 
       height = 150, 
       dpi = 300)



     


