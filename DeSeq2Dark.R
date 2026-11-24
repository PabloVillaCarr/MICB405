library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

#load read count matrix

sample_reads <- read_csv('read_counts/tidy/count_matrix_dark.csv')|>
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
                                     design = ~ time)

dds_matrix$group <- relevel(dds_matrix$group, ref = "infected_Dark_0h") #this sets our "reference" for all groups to be 0h standard light uninfected

#runs the actual DeSeq analysis
dds <- DESeq(dds_matrix)


###At this stage the DeSeq2 part is done

##################################################################################


#PCA analysis
rld <- rlog(dds) #begin by log transforming the count data

PCA_general_Data <- plotPCA(rld,
                            intgroup = c("infection", 'light', 'time'),
                            ntop = 500,
                            returnData = TRUE) #doesn't generate plot, only gives data

#extract percent variation explained of each axis of the PCA data object
percentVar <- round(100 * attr(PCA_general_Data, "percentVar"))
glimpse(percentVar)

#make the actual plot with ggplot

PCA_general_plot <- PCA_general_Data|>
  ggplot(aes(x = PC1, y = PC2, colour = group, shape = infection))+
  geom_point()+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

PCA_general_plot

ggsave('plots/QC/PCA_general_plot_dark.png', plot = PCA_general_plot)



################################################################################

#Make a heatmap of the distances from sample to sample

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Purples")))(255)

# Generate a heatmap using the pheatmap package
sampledist_heatmap <- pheatmap(sample_dist_matrix,
                               clustering_distance_rows = sample_dists,
                               clustering_distance_cols = sample_dists, 
                               col = colours)

ggsave('plots/QC/sampledist_heatmap_dark.png', plot = sampledist_heatmap)


###########################################################################
#Plotting a Z-score heatmap on all genes and on highly variant genes
vsd <- assay(vst(dds)) #normalization and transformation for variance

Z_scores <- t(scale(t(vsd)))|>
  na.omit()#calculate Z scores and omit NA values

#plot pheatmap for all genes
ZscoreHeatmap_general <- pheatmap(Z_scores, main = "Zscore analysis")
ZscoreHeatmap_general

ggsave('plots/Z-scores/Zscore_allgenes_dark.png', plot= ZscoreHeatmap_general)
#extract all the rows with the highest variation
topVarGenes <- rowVars(vsd) |>
  order(decreasing = TRUE) |>
  head(50) #chose the TOP50 most variant genes (roughly 5% of the genome)

#subset the Zscore matrix for the rows we want
Z_TopGenes <- Z_scores[topVarGenes,]

#plot pheatmap of top 50 genes

ZscoreHeatmap_top50 <-pheatmap(Z_TopGenes,
                               main = "Z-score top 50",
                               cluster_cols = FALSE) # I turned off the clustering of groups
ZscoreHeatmap_top50

ggsave('plots/Z-scores/Zscore_top50genes_dark.pdf', plot= ZscoreHeatmap_top50, height = 10)
