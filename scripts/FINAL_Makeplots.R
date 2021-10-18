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
library(gtable)
library(grid)

#==============================================================
# Import data for GO analysis 
#==============================================================

GO_UP <- read.delim("./analysis/GO_analysis/PANTHER_DEseq_kmeans_1v2_sigUP.txt", skip = 11)
GO_DOWN <- read.delim("./analysis/GO_analysis/PANTHER_DEseq_kmeans_1v2_sigDOWN.txt", skip =11)

#===============================================================
# Figure 1 
#===============================================================

PCA_vars = c("Source", "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", 
             "Ulcer_duration_cat", "Bac_prcnt", "cluster_res_all")

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
  theme(legend.position = c(0.2,0.8))+
  labs(color = "% Bacterial Reads")


PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]]<- 
  addSmallLegend(PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]], 
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.2,0.8))+
  labs(color = "IDSA Score")+
  scale_color_manual(values = c("#84D7E1B2", "#FF6F00B2", "#C71000B2"), na.value = "#3F4041B2")


PCA_plots_all[["Ulcer_duration_cat"]]<-
  addSmallLegend(PCA_plots_all[["Ulcer_duration_cat"]],
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.2,0.8))+
  labs(color = "Ulcer Duration")+
  scale_color_manual(values = c("#ADE2D0B2","#FF95A8B2", "#8A4198B2"), na.value = "#999999")



#Arrange the PCs
PCA_grid <- plot_grid(PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]],
                      PCA_plots_all[["Ulcer_duration_cat"]],NULL,
                      PCA_plots_all[["Bac_prcnt"]], ncol = 2,
                      vjust = -0.5, hjust = -0.04
)

# Import flowchart 

flowchart<-
  ggdraw() + 
  draw_image("./analysis/Figures/Sampling_overview.png")+
  theme(plot.margin = margin(0,0,0,0,"cm"))

#Export PDFs
# 
# pdf("./analysis/Figures/Fig1/Fig1_corr.pdf",)
# plt_corr
# dev.off()
# 
# pdf("./analysis/Figures/Fig1/Fig1_loadings.pdf")
# plt_loadings
# dev.off()
# 
# pdf("./analysis/Figures/Fig1/Fig1_screeplot.pdf")
# plt_screeplot
# dev.off()
# 
# pdf("./analysis/Figures/Fig1/Fig1_PCA_BacPrcnt.pdf", width = 2.5, height = 2.5)
# PCA_plots_all[["Bac_prcnt"]]
# dev.off()
# 
# pdf("./analysis/Figures/Fig1/Fig1_PCA_IDSAPEDIS.pdf", width = 2.5, height = 2.5)
# PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]]
# dev.off()
# 
# pdf("./analysis/Figures/Fig1/Fig1_PCA_Duration.pdf", width = 2.5, height = 2.5)
# PCA_plots_all[["Ulcer_duration_cat"]]
# dev.off()

#Export as SVG

svg("./analysis/Figures/Fig1/Fig1_corr.svg", width = 20, height = 20)
plt_corr
dev.off()

svg("./analysis/Figures/Fig1/Fig1_loadings.svg", width = 35, height = 11)
plt_loadings
dev.off()

svg("./analysis/Figures/Fig1/Fig1_screeplot.svg")
plt_screeplot
dev.off()

svg("./analysis/Figures/Fig1/Fig1_PCA_BacPrcnt.svg", width = 2.5, height = 2.5)
PCA_plots_all[["Bac_prcnt"]]
dev.off()

svg("./analysis/Figures/Fig1/Fig1_PCA_IDSAPEDIS.svg", width = 2.5, height = 2.5)
PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]]
dev.off()

svg("./analysis/Figures/Fig1/Fig1_PCA_Duration.svg", width = 2.5, height = 2.5)
PCA_plots_all[["Ulcer_duration_cat"]]
dev.off()


#import plots as svg to ggplot 

temp_scree <- ggdraw() +
  draw_image("./analysis/Figures/Fig1/Fig1_screeplot.svg")

temp_loadings <- ggdraw() +
  draw_image("./analysis/Figures/Fig1/Fig1_loadings.svg")

temp_corr <- ggdraw() + 
  draw_image("./analysis/Figures/Fig1/Fig1_corr.svg")

#Organize the final plot 

pca_plotgrid <- plot_grid(PCA_plots_all[["PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]], 
                     PCA_plots_all[["Ulcer_duration_cat"]], 
                     PCA_plots_all[["Bac_prcnt"]],
                     ncol = 3, labels = c("a", "b", "c"))

top_row <- plot_grid(pca_plotgrid, temp_loadings, nrow = 2, rel_heights = c(1,1), labels = c("", "d"))

bottom <- plot_grid(temp_corr, temp_scree, ncol = 2, rel_widths = c(1,0.6), labels = c("e", "f"))

Fig1_out <- plot_grid(top_row, bottom, nrow = 2)

ggsave("./analysis/Figures/Fig1/Fig1_out.tiff", width = 180, height = 200, units = "mm")
ggsave("./analysis/Figures/Fig1/Fig1_out.pdf", width = 180, height = 200, units = "mm")

#=============================================
#Figure 2: Bacterial Load
#=============================================

#Define a theme for the plots 

fig2_theme <- theme_bw()+
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 6),
    legend.text = element_text(size = 6)
)

#Define colors for relativ abundance plot
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(high_bac_microbiome_top_species$ID)))

#Plot the "relative abundance" of species in samples with high amounts of bac reads
plt_topSpeciesAbundance<- #For those species with greater than 1% rel_bundance and >10% bac:human reads
  ggplot(high_bac_microbiome_top_species, aes(x = Sample_ID, y = Rel_abund, fill = fct_reorder(ID, Rel_abund)))+
  geom_bar(stat = "identity", width = 0.9)+
  scale_fill_manual(values = mycolors)+
  labs(x = "Sample", y = "Relative Activity (%)")+
  guides(fill = guide_legend(title = NULL))+
  theme(legend.position = "bottom",
        legend.key.size = unit(3,"mm"),
        legend.text = element_text(size = 6))

#Plot Cluster vs proportion of bacterial reads 
plt_bac_prop<- #Plots percent of bacterial reads relative to human
ggplot(subset(metadata, Bac_prcnt>0), aes(y=Bac_prcnt, x = cluster_res_all, color = cluster_res_all))+
  geom_jitter(width = 0.1)+
  labs(x = "Cluster", y = "Bacterial:Human Reads (%)")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("#228B22","#A020F0"))+
  fig2_theme

plt_bac_prop <- 
  addSmallLegend(plt_bac_prop,
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.2,0.8)) 

# plt_nSpeciesVsIDSA<-
# ggplot(n_species_1prcntRelAbund, aes(x = PEDIS_IDSA_1uninfected_2mild_3mod_4severe, y = as.numeric(n_species)))+
#   geom_jitter(width = 0.1)+
#   geom_boxplot(aes(alpha = 1))+
#   fig2_theme
# 
# plt_nSpeciesVscluster<-
# ggplot(n_species_1prcntRelAbund, aes(x = cluster_res_all, y = as.numeric(n_species)))+
#   geom_violin(width = 0.1)+
#   geom_boxplot(aes(alpha = 1))+
#   fig2_theme


# GO table

# Format the GO data 
GO_data <- lapply(list(GO_UP, GO_DOWN), function(x){
  data <- as.data.frame(x[,c(1,6,8)]) #Select name, log fold change, padj
  colnames(data) <- c("GO Term", "Fold Enric.", "FDR")
  data[,2] <- as.numeric(data[,2])
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
  data_out <- data[order(data[,2], decreasing = T),] #sort by fold enrichment
  data_out <- data_out[1:5,]
  # plot_out <- tableGrob(head(data), rows = NULL, theme = mytheme)
  # plot_out$widths <- unit(c(47,20,20), "mm")
  return(data_out)
})

# Theme for the table 
mytheme <- gridExtra::ttheme_default(core=list(fg_params=list(hjust = 0, x = 0.08, fontsize = 6),
                                               bg_params=list(fill=c(rep("#A0D7A0", length.out = 5),
                                                                     rep("#E4D2F0", length.out = 5)))),
                                     colhead=list(fg_params=list(hjust = 0, x = 0.08, fontsize = 6),
                                                  bg_params=list(fill = rep("#FFFFFF",1))))



GO_plot_out <- tableGrob(bind_rows(GO_data), rows = NULL, theme = mytheme)


GO_plot_out <- gtable_add_grob(GO_plot_out,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(GO_plot_out), l = 1, r = ncol(GO_plot_out))
GO_plot_out <- gtable_add_grob(GO_plot_out,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(GO_plot_out))



#Format kmeans plot from PCAtools
plt_pca_kmeans<-
plt_pca_kmeans+
  fig2_theme  

plt_pca_kmeans <- 
  addSmallLegend(plt_pca_kmeans,
                 pointSize = 2,textSize = 5, spaceLegend = 0.4)+
  theme(legend.position = c(0.2,0.8))


# ORganize the plot grid 

fig2_topleft <- plot_grid(plt_pca_kmeans, plt_bac_prop, nrow = 2, rel_heights = c(0.7,0.4), labels = c("a", "c"))
fig2_top <- plot_grid(fig2_topleft, GO_plot_out, ncol = 2, rel_widths = c(0.6, 0.5), labels = c("", "b"))
fig2_out<-plot_grid(fig2_top, plt_topSpeciesAbundance, nrow = 2, rel_heights = c(0.8, 0.5), labels = c("", "d"))



#Export Main Figure 

ggsave("./analysis/Figures/Figure2.pdf",
       fig2_out, units = "mm",
       width = 178,
       height = 210,
       dpi = 300,
       bg = "white")

ggsave("./analysis/Figures/Figure2.tiff",
       fig2_out, units = "mm",
       width = 178,
       height = 210,
       dpi = 300,
       bg = "white")


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
warning(paste0("Changing: ", names(SVM_genes), " to ", c("Cluster", "high_bacteria", "IDSA/PEDIS Score"),
               " in SVM_names_new. Make sure names are correct.", sep = "\n"))

#Make plots of Coefficients

coef_samples<- list.files(SVM_dir, pattern = "COEF.*.csv", full.names = T)
coef_data <- lapply(coef_samples, read.delim, sep = ",", col.names = c("Gene", "Coef"))
names(coef_data) <- basename(coef_samples)

SVM_coefs_plots <-
  lapply(seq_along(coef_data), function(x){
    data <- coef_data[[x]]
    #name <- str_extract(names(coef_data)[x], "isCluster_[1-3]")
    #name<-gsub("isCluster_", "Cluster ", name)
    name <-names(coef_data)[x]
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
  pos_coef <- data$Gene[data$Coef<0]
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
  name <- gsub("_[0-9].csv", "", name)
  name <- gsub(".csv", "", name)

})

plots <- lapply(seq_along(SVM_expression_data), function(x){
  n_genes <- length(unique(SVM_expression_data[[x]]$Gene_ID))
  gene_names <- unique(SVM_expression_data[[x]]$Gene_ID)
  facet_annotations <- data.frame(label = gene_names, Gene_ID = gene_names, x = 1.5, y = 26)
  if(n_genes < 10){
    violin_plot<-
      ggplot(SVM_expression_data[[x]],
             aes_string(x=metadata_vars[[x]], y="Norm_expression"))+
      geom_violin(scale = "width", trim = FALSE, aes_string(fill=metadata_vars[[x]]))+
      geom_jitter(position = position_jitter(width=0.05, height = 0), size = 1)+
      scale_y_continuous(limits = c(0,30), breaks = c(seq(5,25,5)))+
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
      scale_y_continuous(limits = c(0,30), breaks = c(seq(5,25,5)))+
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

plots_COEFS<-plot_grid(plotlist = SVM_coefs_plots[c(1, 3:5)], ncol = 4)
expression_plots <- plot_grid(plotlist = as.list(plots[c(1, 3:5)]), nrow = 4, ncol = 1, labels = c("b.", "c.", "d.", "e."))

# Accuracy

plots_accuracy_cluster <-   
  ggdraw() + 
  draw_image("./analysis/Figures/Accuracy_cluster_res_all.png")#+
#  theme(plot.margin = margin(0,0,0,0,"cm"))

plots_accuracy_IDSA <-
  ggdraw() + 
  draw_image("./analysis/Figures/Accuracy_PEDIS_IDSA_1uninfected_2mild_3mod_4severe.png")#+
#  theme(plot.margin = margin(0,0,0,0,"cm"))

plots_accuracy <- plot_grid(plots_accuracy_cluster, plots_accuracy_IDSA, ncol = 2, labels = c("f.", "g."))

Final_Fig5<-
plot_grid(plots_COEFS, expression_plots, plots_accuracy, nrow = 3, rel_heights = c(0.2, 0.6, 0.2), labels = c("a."))


Final_Fig5

ggsave("./analysis/Figures/Figure5.png",
       Final_Fig5, units = "mm",
       width = 180,
       height = 210,
       dpi = 300,
       bg = "white")

ggsave("./analysis/Figures/Figure5.pdf",
       Final_Fig5, units = "mm",
       width = 180,
       height = 210,
       dpi = 300,
       bg = "white")




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
plot_grid(plotlist = annotation_tables[1], ncol = 1, align = "v",
          labels = c("Infection Fingerprint Genes"),
          vjust = 1.1, hjust = -0.08)

annotation_out

ggsave("./analysis/Figures/Fig4.pdf",
       annotation_out, units = "mm",
       width = 90,
       height = 150,
       dpi = 300,
       bg = "white")

ggsave("./analysis/Figures/Fig4.png",
       annotation_out, units = "mm",
       width = 90,
       height = 150,
       dpi = 300,
       bg = "white")

