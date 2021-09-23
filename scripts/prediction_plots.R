
predictions_C3 <- read.csv("./analysis/validation/isCluster_3_predictions.csv")
predictions_C2 <- read.csv("./analysis/validation/isCluster_2_predictions.csv")
predictions_C1 <- read.csv("./analysis/validation/isCluster_1_predictions.csv")


data <- list(predictions_C3, predictions_C2, predictions_C1)
clusters <- c("C3", "C2", "C1")

data_formatted <- lapply(seq_along(data), function(x){
  require(tidyverse)
  
  cluster_name <- clusters[x]
  data_formatted <- data[[x]]

  data_formatted<-
  data_formatted %>% 
    separate(prediction_proba, into = c("not_clust", "is_clust"), sep = ",") %>%
    mutate("is_clust" = gsub("\\]", "", is_clust)) %>%
    mutate("not_clust" = gsub("\\[", "", not_clust))
    
  data_formatted$is_clust <- as.numeric(data_formatted$is_clust)
  data_formatted$not_clust <- as.numeric(data_formatted$not_clust)
  data_formatted$prediction <- as.factor(data_formatted$prediction)
  
  data_formatted$prediction <- recode_factor(data_formatted$prediction, "1" = "is_clust", "0" = "not_clust")
  
  colnames(data_formatted)[c(3,4)]<-
    c(paste0("not_", cluster_name), paste0("is_", cluster_name))

  return(data_formatted)
})

#Combine into a single data set
data_plotting <- 
  data_formatted[[1]] %>%
  left_join(data_formatted[[2]], by = "Sample_ID", suffix = c("C3", "C2")) %>%
  left_join(data_formatted[[3]], by = "Sample_ID")

colnames(data_plotting)[8] <- "predictionC1"


#Add type data 
data_plotting <-
data_plotting %>% 
  left_join(metadata_validation[,c("Sample_ID", "high_bacteria", "Specific_Type", "Bac_prcnt")])

data_plotting$Specific_Type[grep("HH",data_plotting$Sample_ID)] <- "DFU - Clinically Infected"
data_plotting$Specific_Type[grep("SAW_DFU",data_plotting$Sample_ID)] <- "DFU - Clinically Not Infected"
data_plotting$Specific_Type[grep("MW|KK",data_plotting$Sample_ID)] <- "DFU - Infection Unknown"

#Classification to C3 

ggplot(data = data_plotting, aes(x = predictionC3, y = is_C3))+
  geom_jitter(width = 0.05, aes(color = Specific_Type), size = 2.5, alpha = 0.7)+
  geom_violin(scale = "width", width = 0.5, alpha = 0, trim = F)+
  labs(x = "Prediction", y = "Prediction Confidence", title = "Cluster 3 Predictions")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))


ggplot(data = data_plotting, aes(x = predictionC1, y = is_C1))+
  geom_jitter(width = 0.05, aes(color = Specific_Type), size = 2.5, alpha = 0.7)+
  geom_violin(scale = "width", width = 0.5, alpha = 0, trim = F)+
  labs(x = "Prediction", y = "Prediction Confidence", title = "Cluster 1 Predictions")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))


ggplot(data = data_plotting, aes(x = predictionC2, y = is_C2))+
  geom_jitter(width = 0.05, aes(color = "Specific_Type"))+
  geom_violin(scale = "width", width = 0.5, alpha = 0)

ggplot(data = data_plotting, aes(x = predictionC1, y = is_C1))+
  geom_jitter(width = 0.05, aes(color = Specific_Type))+
  geom_violin(scale = "width", width = 0.5, alpha = 0)
