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

PCA_vars = c("Source", "IDSA_SCORE_1to4", "Ulcer_duration_cat", "Bac_prcnt", "cluster_res_all")

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

# Import Dendogram

dend <-ggdendrogram(hclust_avg)

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
                                   "", "e. Heirarchical Clustering"), 
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
ggplot(subset(metadata, Bac_prcnt>0), aes(y=Bac_prcnt, x = fct_reorder(Sample_ID, Bac_prcnt, .desc = T)))+
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


#theme 
mytheme <- gridExtra::ttheme_default(core=list(fg_params=list(hjust = 0, x = 0, fontsize = 9)),
                                     colhead=list(fg_params=list(hjust = 0, x = 0, fontsize = 9))
)



plot_data <- lapply(list(GO_C3UP_shared, GO_C1UP_shared), function(x){
  data <- as.data.frame(x[,c(1,6,7)]) #Select name, log fold change, padj
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
  data <- data[order(data[,2], decreasing = T),] #sort by fold enrichment
  plot_out <- tableGrob(head(data), rows = NULL, theme = mytheme)
  plot_out$widths <- unit(c(47,20,20), "mm") 
  return(plot_out)
})

#Get the contributions to the PCA 

PCA_maincontribs_plots<-plot_contribs(counts_batchnorm_vst, c(1,2), n_contrib = 20)

PCA_maincontribs <- gridExtra::grid.arrange(PCA_maincontribs_plots[[1]], PCA_maincontribs_plots[[2]])

#Fix up PCA plot 
PCA_plots_all[["cluster_res_all"]]<-
addSmallLegend(PCA_plots_all[["cluster_res_all"]], pointSize = 2, 
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


GO_grid <- plot_grid(NULL, plot_data[[1]], NULL, plot_data[[2]], nrow = 4, align = "hv", labels = c("","c. Enriched Pathways - C3", 
                                                                                        "","d. Enriched Pathways C1"),
                     label_x = -0.4,
                     vjust = c(0, -1.8 , 0, 0, -1), 
                     rel_heights = c(0.05,0.35,0.05,0.4))

PCA_grid <- plot_grid(PCA_maincontribs, PCA_plots_all[["cluster_res_all"]], nrow = 2, rel_heights = c(0.6,0.3),
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
warning(paste0("Changing: ", names(SVM_genes), " to ", c("Infection Score", "Cluster1", "Cluster2", "Cluster3"),
               " in SVM_names_new. Make sure names are correct.", sep = "\n"))

#Make plots of Coefficients

coef_samples<- list.files(SVM_dir, pattern = "COEF.*.csv", full.names = T)
coef_data <- lapply(coef_samples, read.delim, sep = ",", col.names = c("Gene", "Coef"))
names(coef_data) <- basename(coef_samples)

SVM_coefs_plots <- 
  lapply(seq_along(coef_data), function(x){
    data <- coef_data[[x]]
    name <- str_extract(names(coef_data)[x], "isCluster_[1-3]")
    name<-gsub("isCluster_", "Cluster ", name)
    print(name)
    plot<-
    ggplot(data, aes(x = Coef, y = fct_reorder(Gene, Coef)))+
      geom_bar(stat = "identity")+
      xlim(-0.25, 0.25)+
      labs(title = name, x = "Coefficient")+
      theme(axis.text = element_text(size = 5),
            axis.title.y = element_blank())
    return(plot)
    
  })

SVM_expression_data <- lapply(coef_data, function(x){
  data <- x
  pos_coef <- data$Gene[data$Coef>0]
  #Get expression data 
  expression_values <- subset(counts_batchnorm_vst, row.names(counts_batchnorm_vst) %in% pos_coef)
  expression_values <- as.data.frame(expression_values)
  expression_values <-
    expression_values %>%
    rownames_to_column("Gene_ID") %>%
    pivot_longer(cols = -1, names_to = "Sample_ID", values_to = "Norm_expression") %>%
    left_join(metadata, by = "Sample_ID")
})


#Format names for x axis
metadata_vars<-
lapply(names(SVM_expression_data), function(x){
  name <- gsub("COEF_", "", x)
  name <- gsub("_0.csv", "", name)
  name <- gsub(".csv", "", name)
  if(grepl("IDSA", name)){
    name <- "IDSA_SCORE_1to4"
    return(name)
  }
  else{return(name)}
})

plots <- lapply(seq_along(SVM_expression_data), function(x){
  n_genes <- length(unique(SVM_expression_data[[x]]$Gene_ID))
  gene_names <- unique(SVM_expression_data[[x]]$Gene_ID)
  facet_annotations <- data.frame(label = gene_names, Gene_ID = gene_names, x = 0.9, y = 24)
  if(n_genes < 7){
    violin_plot<-
      ggplot(SVM_expression_data[[x]], 
             aes_string(x=metadata_vars[[x]], y="Norm_expression"))+
      geom_violin(scale = "width", trim = FALSE, aes_string(fill=metadata_vars[[x]]))+
      geom_jitter(position = position_jitter(width=0.05, height = 0), size = 1)+
      scale_y_continuous(limits = c(0,25), breaks = c(seq(5,25,5)))+
      theme_half_open()+
      theme(legend.position = "none",
            panel.spacing = unit(0, "lines"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8))+
      ylab("Normalized\nExpression")+
      facet_grid(~ Gene_ID, scales = "free")+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())+
      geom_text(data = facet_annotations, 
                mapping = aes(x = x, 
                              y = y,
                              label = Gene_ID,),
                size = 2)
  }
  else{
    violin_plot<-
      ggplot(SVM_expression_data[[x]], 
             aes_string(x=metadata_vars[[x]], y="Norm_expression"))+
      geom_violin(scale = "width", trim = FALSE, aes_string(fill=metadata_vars[[x]]))+
      geom_jitter(position = position_jitter(width=0.05, height = 0), size = 1)+
      scale_y_continuous(limits = c(0,25), breaks = c(seq(5,25,5)))+
      xlab(metadata_vars[[x]])+
      theme_half_open()+
      theme(legend.position = "none",
            panel.spacing = unit(0, "lines"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8))+
      ylab("Normalized\nExpression")+
      facet_wrap(~ Gene_ID, ncol = 7)+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())+
      geom_text(data = facet_annotations, 
                mapping = aes(x = x, 
                              y = y,
                              label = Gene_ID,),
                size = 2)
  }
})

#Build the plots 


plots_COEFS<-plot_grid(plotlist = SVM_coefs_plots[4:6], ncol = 3)

#plots_test<-as.list(plots[4], NULL)
expression_plots1 <- plot_grid(plotlist = as.list(plots[4], NULL), ncol = 2, rel_widths = c(0.73, 0.27))
expression_plots2 <- plot_grid(expression_plots1, plotlist = plots[5:6], 
                               nrow = 3, rel_heights = c(0.20,0.3,0.3), 
                               labels = c("b.", "c.", "d."), vjust = c(0,0,0))
 

Final_Fig4<-
plot_grid(plots_COEFS, expression_plots2, nrow = 2, rel_heights = c(0.3, 0.7), labels = c("a.") )

ggsave("./analysis/Figures/Figure4.png",      
       Final_Fig4, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)

ggsave("./analysis/Figures/Figure4.pdf",      
       Final_Fig4, units = "mm", 
       width = 180, 
       height = 200, 
       dpi = 300)
  
### Figure 4a #############################################


#============================================================
# Annotation Tables Out 
#===========================================================
#Only include annotations for clusters

#Get the annotations for the good genes (takes several minutes to run!)

gene_products<-
  lapply(SVM_genes, function(x){
    data <- unlist(x)
    get_gene_names(data, annotation_path, "exon", "mRNA")
  })


annotation_tables <- lapply(seq_along(gene_products), function(x){
  name <- names(gene_products)[x]
  data <- gene_products[[x]]
  
  
  
  #Add line break if annotation longer than 45 chars 
  data$Product <- gsub('(.{1,30})(\\s|$)', '\\1\n', data$Product)
  data$Product <- gsub('\n$', '', data$Product)
  table <- tableGrob(data, rows = NULL,
                     theme = ttheme_default(core=list(fg_params=list(cex = 0.8, hjust=0, x=0.1)),
                                            padding = unit(c(0.8,1.3), "mm")))
})


annotation_out<-
plot_grid(plotlist = annotation_tables[2:4], ncol = 3, align = "hv",
          labels = c("a. Cluster 1", "b. Cluster 2","c. Cluster 3"),
          vjust = 1.1, hjust = -0.08)

annotation_out

ggsave("./analysis/Figures/Fig_Annotation.pdf",      
       annotation_out, units = "mm", 
       width = 250, 
       height = 183, 
       dpi = 300)
ggsave("./analysis/Figures/Fig_Annotation.png",      
       annotation_out, units = "mm", 
       width = 250, 
       height = 183, 
       dpi = 300)

