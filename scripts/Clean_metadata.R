#=================================
# Script to format the metadata
#=================================

#Note: The same data was used for both Heravi et al(DOI: 10.3389/fmicb.2021.613697)
#      and Radzieta et al (DOI:10.1038/s41522-021-00202-x). Sample IDs are listed differently
#      between these studies. A key is provided in "./data/Example_data/External_data_keys.csv"

#Import libraries
library(tidyverse)

#==========================================
# Combine External and Internal Metadata
#==========================================

# Import the data 

keys <- read.delim("./data/Example_data/External_data_keys.tsv", sep = "\t")
external_metadata <- read.delim("./data/Example_data/External_metadata_Heravi_Radzieta.tsv", sep = "\t")
internal_metadata <- read.delim("./data/Example_data/Internal_metadata.tsv", sep = "\t")
sequence_metadata <- read.delim("./data/Example_data/Sequencing_metadata.csv", sep = ",")

# The first 22 columns are simiar between data sets

external_metadata_matching <- external_metadata[,1:22]
internal_metadata_matching <- internal_metadata[,1:22]

#### Format Column names for consistency  

colnames(external_metadata_matching)
colnames(external_metadata_matching) <- c("Sample_ID",  #1
                                 "DOB", #2 
                                 "Age", #3
                                 "Gender_0M_1F", #4
                                 "Type_of_Diabetes_T1_T2", #5
                                 "Duration_of_Diabetes_years", #6
                                 "Peripheral_Neuropathy", #7
                                 "PAK", #8
                                 "CKD_stage5", #9
                                 "IHD", #10
                                 "CCF", #11
                                 "Ulcer_duration_cat", #12
                                 "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", #13
                                 "WCC_10e9perL", #14
                                 "CRP_mgperL", #15
                                 "ESR_mLpermin", #16
                                 "Neutrophils_10e9perL", #17
                                 "HBA1c", #18
                                 "HBA1c_IFCC", #19
                                 "Culture.results", #20
                                 "WiFI.Score",
                                 "Sample_collection"
                                 )

colnames(internal_metadata_matching) <-  c("Sample_ID",  #1
                                           "DOB", #2 
                                           "Age", #3
                                           "Gender_0M_1F", #4
                                           "Type_of_Diabetes_T1_T2", #5
                                           "Duration_of_Diabetes_years", #6
                                           "Peripheral_Neuropathy", #7
                                           "PAK", #8
                                           "CKD_stage5", #9
                                           "IHD", #10
                                           "CCF", #11
                                           "ABI", #12
                                           "Toe_pressure_mmHg", #13
                                           "Ulcer_duration_weeks", #14
                                           "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", #15
                                           "WCC_10e9perL", #16
                                           "CRP_mgperL", #17
                                           "ESR_mLpermin", #18
                                           "Neutrophils_10e9perL", #19
                                           "HBA1c", #20
                                           "HBA1c_IFCC", #21
                                           "Culture.results") #22


#Convert all columns to characters for now 

internal_metadata_matching <-
  internal_metadata_matching %>%
    mutate(across(everything(), as.character))

external_metadata_matching <-
  external_metadata_matching %>%
  mutate(across(everything(), as.character))

keys <- 
  keys %>%
    mutate(across(everything(), as.character))

# Add Column for source
external_metadata_matching$Source <- "HH" #Heravi et al 
internal_metadata_matching$Source <- "CBC" #Costerton Biofilm Center

# For External Data, filter and add keys 
external_metadata_matching<-
external_metadata_matching %>%
  filter(Sample_ID %in% keys$Radzieta_article_ID) %>%
  left_join(keys, by = c("Sample_ID" = "Radzieta_article_ID")) %>%
  mutate("Sample_ID" = Heravi_article_ID) #Use article IDs

#Add "P" to Sample ID for internal data 

internal_metadata_matching$Sample_ID[internal_metadata_matching$Source == "CBC"] <- 
  paste0("P",internal_metadata_matching$Sample_ID)

# Join the data into a single data set 

combined_metadata <- bind_rows(external_metadata_matching, 
                               internal_metadata_matching)

#==================================
# Formatting Metadata 
#==================================

# Fix gender column
combined_metadata$Gender_0M_1F <-
  if_else(combined_metadata$Gender_0M_1F %in% c("M", "0"), "M", "F")

# Convert days to weeks in Ulcer Duration

combined_metadata$Ulcer_duration_weeks[combined_metadata$Ulcer_duration_weeks=="5 days"]<-as.character(5/7)
combined_metadata$Ulcer_duration_weeks[combined_metadata$Ulcer_duration_weeks=="3 days"]<-as.character(3/7)


# Remove ABI and Toe Pressure columns since not in both data 

combined_metadata <- combined_metadata %>% select(!c("ABI","Toe_pressure_mmHg"))

#Replace Blank Cells with NA 

combined_metadata <- combined_metadata %>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

# Make Ulcer duration cat where 0 = <2 weeks, 1 = 2-4 weeks , and 2 = >4 weeks 

combined_metadata <-
  combined_metadata %>%
    mutate("Ulcer_duration_weeks" = as.numeric(Ulcer_duration_weeks)) %>%
    mutate("Ulcer_duration_cat_CBC" = replace(Ulcer_duration_weeks, Ulcer_duration_weeks < 2, "0")) %>%
    mutate("Ulcer_duration_cat_CBC" = replace(Ulcer_duration_cat_CBC, Ulcer_duration_weeks >= 2 &  Ulcer_duration_weeks <= 6, "1")) %>%
    mutate("Ulcer_duration_cat_CBC" = replace(Ulcer_duration_cat_CBC, Ulcer_duration_weeks > 6, "2"))

combined_metadata$Ulcer_duration_cat <- c(combined_metadata$Ulcer_duration_cat[combined_metadata$Source == "HH"],
                                          combined_metadata$Ulcer_duration_cat_CBC[combined_metadata$Source == "CBC"])

combined_metadata <- combined_metadata %>% select(!Ulcer_duration_cat_CBC)

# Format Column types 

combined_metadata <-
  combined_metadata %>%
    mutate(across(c("Gender_0M_1F", "Peripheral_Neuropathy", "PAK", "CKD_stage5", 
                    "IHD", "CCF", "Type_of_Diabetes_T1_T2","Ulcer_duration_cat", 
                    "PEDIS_IDSA_1uninfected_2mild_3mod_4severe", "Source"), as.factor)) %>%
    mutate(across(c("Age", "Duration_of_Diabetes_years", "Ulcer_duration_weeks", "WCC_10e9perL", "CRP_mgperL", "ESR_mLpermin",
                    "Neutrophils_10e9perL", "HBA1c", "HBA1c_IFCC"), as.numeric))

# Add Sequence Metadata to combined metadata 

combined_metadata<-sequence_metadata %>%
                      left_join(combined_metadata, by = "Sample_ID")

# Write out metadata as csv

write.csv(combined_metadata, "./data/Example_data/Analysis_metadata.csv", row.names = F)










