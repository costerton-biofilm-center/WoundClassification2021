##=========================================
## Make predictions based on clinical data
##=========================================
library(dplyr)
library(tidyr)
library(cowplot)
#C1 - > HealthySkin
#C2 - > DFU but <10% bacteria, <4 weeks 
#C3 - > DFU but greater than 10% bac

prediction_guesses <- data.frame("Sample_ID" = metadata_validation_filtered$Sample_ID,
                                 "Train_test" = metadata_validation_filtered$Source,
                                 "Ulcer_duration_cat" = metadata_validation_filtered$Ulcer_duration_cat,
                                 "Bac_prcnt" = metadata_validation_filtered$Bac_prcnt,
                                 "Source" = metadata_validation_filtered$Source,
                                 "Type" = metadata_validation_filtered$Type,
                                 "IDSA_SCORE_1to4" = metadata_validation_filtered$IDSA_SCORE_1to4,
                                 "Specific_Type" = metadata_validation_filtered$Specific_Type,
                                 stringsAsFactors = F)

prediction_guesses<-
prediction_guesses %>% 
  mutate(Train_test = ifelse(Train_test == "CBC", "Train", "Test")) %>%
  mutate( cluster_guess = ifelse(Type == "DFU" & Bac_prcnt > 10, "3", "")) %>%
  mutate( cluster_guess = ifelse(Type == "DFU" & Bac_prcnt < 10, "2", cluster_guess)) %>%
  mutate( cluster_guess = ifelse(Type != "DFU" & Bac_prcnt < 10, "1",  cluster_guess)) %>%
  mutate( Specific_Type = ifelse(Type == "DFU" & Bac_prcnt < 10, "DFU - Not Infected", Specific_Type)) %>%
  mutate( Specific_Type = ifelse(Type == "DFU" & Bac_prcnt > 10, "DFU - Infected",  Specific_Type))



##===========================
## Analyze Prediction CSV
##============================

# Read in csv 

model_predictions <- read.csv("./analysis/validation/predictions.csv")

# Split the prediction probability string
model_predictions<-
model_predictions %>%
  separate(prediction_proba., c("prob_C1", "prob_C2", "prob_C3"), sep = ",") 

# Remove brackets, whitespace, and make numeric
model_predictions[,c("prob_C1", "prob_C2", "prob_C3")]<-
  apply(model_predictions[,c("prob_C1", "prob_C2", "prob_C3")], 2, function(x){
    no_brackets <- gsub("\\[|\\]", "", x)
    no_brack_no_white<- gsub(" ", "", no_brackets)
    out_numeric <- as.numeric(no_brack_no_white)
    return(out_numeric)
  })


#Combine predictions with prediction_guesses 

prediction_guesses<-
  prediction_guesses %>%
  filter(Train_test != "Train" & cluster_guess != "NA") %>% 
  left_join(model_predictions, by = "Sample_ID")



#Calculate accuracy - All Data

accuracy_All <- calc_Accuracy(prediction_guesses, "cluster_guess",
                              "prediction")

#Calculate DFU - C2 accuracy 

accuracy_noninfDFU <- calc_Accuracy(subset(prediction_guesses, cluster_guess == "2"),
                                 "cluster_guess",
                                 "prediction")
                            

#Infected vs non-infected DFU 

accuracy_InfectedvsNonInfectedDFU <- calc_Accuracy(subset(prediction_guesses, 
                                                   cluster_guess == "3"|prediction == "3"),
                                            "cluster_guess", 
                                            "prediction")

#Healthy as C1  

accuracy_Healthy  <- calc_Accuracy(subset(prediction_guesses, 
                                                   Type == "HealthySkin"),
                                                   "cluster_guess", 
                                                   "prediction")

#uninfected DFU as C2

accuracy_noninfectedDFUC2 <- calc_Accuracy(subset(prediction_guesses, 
                                                  cluster_guess == "2"),
                                           "cluster_guess", 
                                           "prediction")


#Get prob for predicted cluster

prediction_guesses$probab_selection <- 
  apply(prediction_guesses[,c("prob_C1", "prob_C2", "prob_C3")], 1, max)

#Create column to show which were predicted accurately 

prediction_guesses$correct <- c(prediction_guesses$cluster_guess==prediction_guesses$prediction)

#Plot confidence by cluster 

ggplot(subset(prediction_guesses), aes(x = as.factor(prediction), y = probab_selection))+
  geom_violin(scale = "count")+
  geom_jitter(aes(shape = Specific_Type, color = Specific_Type), position = position_jitter(width=0.05, height = 0))+
  scale_fill_discrete(name = "Tissue Type")+
  ylim(0.5, 1.0)+
  ggtitle("Prediction Confidence")+
  ylab("Probability")+
  xlab("Cluster")+
  theme_half_open()



#Get Correct/Incorrect
pred_summary<-
prediction_guesses %>% 
  select(Specific_Type, cluster_guess, prediction, probab_selection) %>%
  group_by(Specific_Type, cluster_guess) %>%
  summarize("n_correct" = sum(cluster_guess == prediction),
            "n_incorrect" = sum(cluster_guess != prediction),
            "prcnt_correct" = n_correct/(n_correct + n_incorrect)*100,
            "Mean Assignment Confidence" = mean(probab_selection),
            "SD" = sd(probab_selection))







