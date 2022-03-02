#Plot PCA 

library(PCAtools)

# Run PCA analysis 
PCA_metadata <- metadata 

#Refactor high bacteria category 
PCA_metadata$high_bacteria_1High0Low <- ifelse(PCA_metadata$high_bacteria_1High0Low=="1", ">=10 %", "<10 %")

#Rename Ulcer duration cat
PCA_metadata$Ulcer_duration_cat <- as.character(PCA_metadata$Ulcer_duration_cat)
PCA_metadata$Ulcer_duration_cat[PCA_metadata$Ulcer_duration_cat==0]  <- "<2 weeks"
PCA_metadata$Ulcer_duration_cat[PCA_metadata$Ulcer_duration_cat==1] <- "2 to 6 weeks"
PCA_metadata$Ulcer_duration_cat[PCA_metadata$Ulcer_duration_cat==2] <- "<6 weeks"

#Metadata row names and count col names should match 
row.names(PCA_metadata) <- metadata$Sample_ID

#Fix Long metadata names 
colnames(PCA_metadata)[grep("PEDIS", colnames(PCA_metadata))] <- "PEDIS/IDSA"

#Run the PCA
pca_data <- pca(counts_batchnorm_vst, metadata = PCA_metadata)

#Screeplot
plt_screeplot<-
PCAtools::screeplot(pca_data, components = getComponents(pca_data, 1:20),
                    hline = 90, vline = 16, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = F) +
  geom_label(aes(8, 90, label = '90% explained variation', vjust = -1, size = 8))

#plot Loadings
plt_loadings<-
plotloadings(pca_data, rangeRetain = 0.05, components = 1:3, labSize = 10, legendPosition = "none",
             axisLabSize = 38,
             widthConnectors = 2,
             labvjust = 0.5,
             labhjust = 0,
             lengthConnectors = unit(0.05, "npc"),
             positionConnectors = "right",
             typeConnectors = "open",
             shapeSizeRange = c(20, 20))+
  scale_y_reverse()+
  scale_x_discrete(expand = expansion(add=c(0.4,1)))


#Highlight K-means Clusters
plt_pca_kmeans<-
biplot(pca_data,
       colby = 'cluster_res_all', colkey = c('1' = 'forestgreen', '2' = 'purple'),
       colLegendTitle = 'K-means cluster',
       # encircle config
       encircle = TRUE, encircleFill = FALSE,
       encircleAlpha = 0.5, encircleLineSize = 3,
       hline = 0, vline = c(-25, 0, 25),
       #legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)


#EigenCor PLot 
plt_corr<-
eigencorplot(pca_data,
             components = getComponents(pca_data, 1:6),
             metavars = c("Ulcer_duration_cat", "PEDIS/IDSA", 
                          "WCC_10e9perL", "CRP_mgperL", "ESR_mLpermin","Neutrophils_10e9perL",
                          "HBA1c", "Bac_prcnt"),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             #col = c('darkgreen', 'forestgreen','gold', 'cornsilk1','white'),
             cexCorval = 3,
             cexTitleX = 3,
             cexLabX = 3,
             cexLabY = 3,
             cexMain = 3,
             cexLabColKey = 3,
             fontCorval = 2,
             posLab = 'bottomleft',
             posColKey = 'top',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Spearman ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'spearman',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))


