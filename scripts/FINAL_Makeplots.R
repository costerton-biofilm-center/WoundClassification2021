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
library(ggdendro)
library(RColorBrewer)

#===============================================================
# Figure 1 
#===============================================================

PCA_vars = c("Source", "IDSA_SCORE_1to4", "Ulcer_duration_cat", "Bac_prcnt")

PCA_plots_all<-
  lapply(PCA_vars, function(metadata_cat) {
    plot<-
      plot_PCA(counts_batchnorm_vst, metadata, c(1,2), 
               metadata_cat, normalize = F)
    # ggsave(paste0("./analysis/Figures/PCA_All_", metadata_cat, ".png"),
    #        width = 5, height = 5)
    return(plot+theme(legend.position = "none"))
  })

#Assign names to the list outputs for each access
names(PCA_plots_all)<-PCA_vars

# Manually format the PCAs to make pretty 

PCA_plots_all[["Bac_prcnt"]]<- 
  addSmallLegend(PCA_plots_all[["Bac_prcnt"]], 
                 pointSize = 2, textSize = 5, 
                 spaceLegend = 0.4, scale = "continuous")+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "% Bacterial Reads")

PCA_plots_all[["IDSA_SCORE_1to4"]]<- 
  addSmallLegend(PCA_plots_all[["IDSA_SCORE_1to4"]], 
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "IDSA Score")+
  scale_color_manual(values = c("#84D7E1B2", "#FF6F00B2", "#C71000B2"), na.value = "#3F4041B2")

PCA_plots_all[["Ulcer_duration_cat"]]<-
  addSmallLegend(PCA_plots_all[["Ulcer_duration_cat"]],
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.8,0.8))+
  labs(color = "Ulcer Duration")+
  scale_color_manual(values = c("#ADE2D0B2","#FF95A8B2", "#8A4198B2"), na.value = "#999999")

# Import flowchart 

flowchart<-
  ggdraw() + 
  draw_image("./analysis/Figures/Sampling_overview.png")+
  theme(plot.margin = margin(0,0,0,0,"cm"))

# Import hDendogram

dend <-ggdendogram(hculst_avg)

#Build Plot 

PCA_grid <- plot_grid(PCA_plots_all[["IDSA_SCORE_1to4"]],
                      PCA_plots_all[["Ulcer_duration_cat"]],
                      PCA_plots_all[["Bac_prcnt"]], ncol = 3,
                      labels = c("b. Infection Score", 
                                 "c. Ulcer Duration", 
                                 "d. Bac:Human Reads"), 
                      vjust = -0.5, hjust = -0.04
)

Final_Fig1 <- plot_grid(flowchart, PCA_grid, dend, nrow = 3, 
                        labels = c("a. Design Flowchart",
                                   "", "c. Heirarchical Clustering"), 
                        hjust = -0.04, vjust = c(1.5,0,0.5)
                        
)

# Export Plot 

ggsave("./analysis/Figures/Figure1.png",      
       Final_Fig1, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)

#Also save as svg for editing 
ggsave("./analysis/Figures/Figure1.pdf",      
       Final_Fig1, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)

#=============================================
#Figure 2: Bacterial Load  
#=============================================

mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(high_bac_microbiome_top_species$ID)))

plot_high_bac<-
  ggplot(high_bac_microbiome_top_species, aes(x = Sample_ID, y = Rel_abund, fill = fct_reorder(ID, Rel_abund)))+
  geom_bar(stat = "identity", width = 0.9)+
  scale_fill_manual(values = mycolors)+
  labs(x = "Sample", y = "Relative Activity (%)")+
  guides(fill = guide_legend(title = ""))+
  theme(legend.text = element_text(size = 8))



plot_bac_prop<-
ggplot(subset(metadata_kmeans, Bac_prcnt>0), aes(y=Bac_prcnt, x = fct_reorder(Sample_ID, Bac_prcnt, .desc = T)))+
  geom_bar(stat = "identity")+
  labs(x = "Sample", y = "Bacterial:Human Reads (%)")+
  theme(axis.text.x = element_text(angle = 90))


fig2<-
plot_grid(NULL, plot_bac_prop, NULL, plot_high_bac, nrow = 4, 
          labels = c("","a. Percentage of Bac:Human Reads","", 
                     "b. Relative Activity: # of Reads(Species):# of Reads(Bacteria)"),
          rel_heights = c(0.05, 0.5, 0.08, 0.5),
          vjust = c(0, -0.7, 0, -1.8) , hjust = -0.05)

ggsave("./analysis/Figures/Figure2.pdf",      
       fig2, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)

ggsave("./analysis/Figures/Figure2.png",      
       fig2, units = "mm", 
       width = 180, 
       height = 180, 
       dpi = 300)

#=============================================#
# Figure 3 Kmeans, PCA contibs, GO            #
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
PCA_plots_all[[20]]<-
addSmallLegend(PCA_plots_all[[20]], pointSize = 2, 
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
PCA_grid <- plot_grid(PCA_maincontribs[[2]], PCA_plots_all[[20]], nrow = 2, rel_heights = c(0.6,0.3),
                      labels = c("a.", "b."), label_y = c(1, 1.1))

FINAL_Fig3<-
plot_grid(PCA_grid, GO_grid,  ncol = 2)


ggsave("./analysis/Figures/Figure3.png",      
       FINAL_Fig3, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)

ggsave("./analysis/Figures/Figure3.pdf",      
       FINAL_Fig3, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)

#==========================================================#
# Figure  4 - Predictor Genes Table and Expression Plots   #
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


### Figure 4a #############################################

#Get the annotations for the good genes (takes several minutes to run!)

gene_products<-
lapply(SVM_genes, function(x){
  data <- unlist(x)
  get_gene_names(data, annotation_path, "exon", "mRNA")
})


### Figure 4b ###############################################

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


### Combine 4a and 4b for Final Figure

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
  
  
  ggsave(paste0("./analysis/Figures/Figure4_", name, ".png"),
         out_plot, units = "mm",
         width = 300,
         height = 300,
         dpi = 300) 
  
  ggsave(paste0("./analysis/Figures/Figure4_", name, ".pdf"),
         out_plot, units = "mm",
         width = 300,
         height = 300,
         dpi = 300) 
})


######################################
# Figure 5 - SVC Results Figure     ##
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


cluster_plots <- plot_grid(plotlist = SVC_plots[c(1,2,3)], 
                           ncol = 3,
                           labels = c("a. Cluster 1", "b. Cluster 2", "c. Cluster 3"),
                           align = "hv")

cluster_plots

ggsave("./analysis/Figures/Figure5.png",      
       cluster_plots, units = "mm", 
       width = 300, 
       height = 150, 
       dpi = 300)

ggsave("./analysis/Figures/Figure5.pdf",      
       cluster_plots, units = "mm", 
       width = 300, 
       height = 150, 
       dpi = 300)

###############################
#Figure 6 Model Validations
###############################

#Summary Stats
pred_summary<-
  prediction_guesses %>% 
  select(Specific_Type, cluster_guess, prediction, probab_selection) %>%
  group_by(Specific_Type, cluster_guess) %>%
  summarize("n_correct" = sum(cluster_guess == prediction),
            "n_incorrect" = sum(cluster_guess != prediction),
            "prcnt_correct" = round(n_correct/(n_correct + n_incorrect)*100),
            "Mean Assignment \nConfidence" = round(mean(probab_selection), digits = 3),
            "SD" = round(sd(probab_selection), digits = 2))


fig6_theme <- gridExtra::ttheme_default(core=list(fg_params=list(hjust = 0, x = 0.1, fontsize = 9),
                                               padding=unit(c(5, 5), "mm")),
                                     colhead=list(fg_params=list(hjust = 0, x = 0.1, fontsize = 9))
)

#To Grob Table 

fig6a<-tableGrob(pred_summary, theme = fig6_theme, rows = NULL)


#Plot confidence by cluster 
prob_plot<-
ggplot(subset(prediction_guesses), aes(x = as.factor(prediction), y = probab_selection))+
  geom_violin(scale = "count")+
  geom_jitter(aes(shape = Specific_Type, color = Specific_Type), position = position_jitter(width=0.05, height = 0))+
  scale_fill_discrete(name = "Tissue Type")+
  ylim(0.2, 1.0)+
  ggtitle("Prediction Confidence")+
  ylab("Probability")+
  xlab("Cluster")+
  theme_half_open()

FINAL_Fig6<-
plot_grid(fig6a, prob_plot, nrow = 2, rel_heights = c(0.3, 0.4), labels = c("a", "b"))

FINAL_Fig6

ggsave("./analysis/Figures/Figure6.pdf",      
       FINAL_Fig6, units = "mm", 
       width = 180, 
       height = 150, 
       dpi = 300)

ggsave("./analysis/Figures/Figure6.png",      
       FINAL_Fig6, units = "mm", 
       width = 180, 
       height = 150, 
       dpi = 300)
