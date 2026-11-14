library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(biomaRt)
library(EnhancedVolcano)

#load read count matrix

sample_reads <- read_csv('read_counts/tidy/count_matrix.csv')|>
  column_to_rownames('Geneid')|>
  as.matrix()

#prepare metadata
samples <- tibble(sample = colnames(sample_reads)) #get all sample names

conditions <- c('infection', 'light', 'time', 'rep')  #set all the conditions we need

columndata <- samples |>
  extract(sample,
          into = conditions,
          regex = "^(control|infected)_(Dark|High_Light|Standard_Light)_(\\d+h)_(\\d+)$",
          remove = FALSE)|>
  mutate(
    infection = factor(infection, levels = c('control', 'infected')),
    light = factor (light, levels = c('Standard_Light', 'High_Light', 'Dark')),
    time = factor(time, levels = c('0h', '4h', '8h')),
    rep = factor(rep),
    group = factor(paste(infection, light, time, sep = '_'))
  )|>
  column_to_rownames('sample')

#create DeSeq2 object

dds_matrix <- DESeqDataSetFromMatrix(countData = sample_reads,
                                     colData = columndata,
                                     design = ~ group)

dds_matrix$group <- relevel(dds_matrix$group, ref = "control_Standard_Light_0h") #this sets our "reference" for all groups to be 0h standard light uninfected

#runs the actual DeSeq analysis
dds <- DESeq(dds_matrix)


###At this stage the DeSeq2 part is done

#PCA analysis

