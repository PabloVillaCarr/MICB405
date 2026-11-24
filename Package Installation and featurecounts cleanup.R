#install all the libraries we need

# copypaste the following into the terminal without the #
# if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")



#install.packages("tidyverse")


#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("DESeq2")

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("biomaRt")

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")


#Open libraries
library(tidyverse)

counts <- read_tsv('read_counts/counts.txt', skip = 1)

#cleaning up the column names and removing unnecessary columns
sample_names <- c(
  "control_High_Light_0h_1", "control_High_Light_0h_2", "control_High_Light_0h_3",
  "control_High_Light_4h_1", "control_High_Light_4h_2", "control_High_Light_4h_3",
  "control_High_Light_8h_1", "control_High_Light_8h_2", "control_High_Light_8h_3",
  "control_Standard_Light_0h_1", "control_Standard_Light_0h_2", "control_Standard_Light_0h_3",
  "control_Standard_Light_4h_1", "control_Standard_Light_4h_2", "control_Standard_Light_4h_3",
  "control_Standard_Light_8h_1", "control_Standard_Light_8h_2", "control_Standard_Light_8h_3",
  "infected_Dark_0h_1", "infected_Dark_0h_2", "infected_Dark_0h_3",
  "infected_Dark_4h_1", "infected_Dark_4h_2", "infected_Dark_4h_3",
  "infected_Dark_8h_1", "infected_Dark_8h_2", "infected_Dark_8h_3",
  "infected_High_Light_0h_1", "infected_High_Light_0h_2", "infected_High_Light_0h_3",
  "infected_High_Light_4h_1", "infected_High_Light_4h_2", "infected_High_Light_4h_3",
  "infected_High_Light_8h_1", "infected_High_Light_8h_2", "infected_High_Light_8h_3",
  "infected_Standard_Light_0h_1", "infected_Standard_Light_0h_2", "infected_Standard_Light_0h_3",
  "infected_Standard_Light_4h_1", "infected_Standard_Light_4h_2", "infected_Standard_Light_4h_3",
  "infected_Standard_Light_8h_1", "infected_Standard_Light_8h_2", "infected_Standard_Light_8h_3"
)

colnames(counts)[7:ncol(counts)] <- sample_names

counts_clean <- counts |>
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) #all samples

counts_clean_nodark <- counts_clean|>
  select(-infected_Dark_0h_1,
         -infected_Dark_0h_2,
         -infected_Dark_0h_3,
         -infected_Dark_4h_1,
         -infected_Dark_4h_2,
         -infected_Dark_4h_3,
         -infected_Dark_8h_1,
         -infected_Dark_8h_2,
         -infected_Dark_8h_3)
counts_clean_dark <- counts_clean|>
  select(Geneid,
         infected_Dark_0h_1,
         infected_Dark_0h_2,
         infected_Dark_0h_3,
         infected_Dark_4h_1,
         infected_Dark_4h_2,
         infected_Dark_4h_3,
         infected_Dark_8h_1,
         infected_Dark_8h_2,
         infected_Dark_8h_3)

#for high light at 4 hours
counts_clean_h4h <- counts_clean|>
  select(Geneid,
         control_High_Light_4h_1,
         control_High_Light_4h_2,
         control_High_Light_4h_3,
         infected_High_Light_4h_1,
         infected_High_Light_4h_2,
         infected_High_Light_4h_3)
#export into a nice csv for downstream use
write_csv(counts_clean, 'read_counts/tidy/count_matrix.csv')

write_csv(counts_clean_nodark, 'read_counts/tidy/count_matrix_nodark.csv')

write_csv(counts_clean_h4h, 'read_counts/tidy/count_matrix_h4h.csv')

write_csv(counts_clean_dark, 'read_counts/tidy/count_matrix_dark.csv')

