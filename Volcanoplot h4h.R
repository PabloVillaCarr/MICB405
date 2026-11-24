library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggrepel)
#load read count matrix

sample_reads <- read_csv('read_counts/tidy/count_matrix_h4h.csv')|>
  column_to_rownames('Geneid')|>
  as.matrix()

GO_names <- read_tsv('prochlorococcus_marinus_gene_ontology.tsv', col_names = c('Geneid', 'GO'))

#prepare metadata
samples <- tibble(sample = colnames(sample_reads)) #get all sample names

conditions <- c('infection', 'light', 'time', 'rep')  #set all the conditions we need

columndata <- samples |>
  extract(sample,
          into = conditions,
          regex = "^(control|infected)_(High_Light|Standard_Light)_(\\d+h)_(\\d+)$",
          remove = FALSE)|>
  mutate(
    infection = factor(infection, levels = c('control', 'infected')),
    light = factor (light, levels = c('Standard_Light', 'High_Light')),
    time = factor(time, levels = c('0h', '4h', '8h')),
    rep = factor(rep),
    group = factor(paste(infection, light, time, sep = '_'))
  )|>
  column_to_rownames('sample')

#create DeSeq2 object

dds_matrix <- DESeqDataSetFromMatrix(countData = sample_reads,
                                     colData = columndata,
                                     design = ~ group)

dds_matrix$group <- relevel(dds_matrix$group, ref = "control_High_Light_4h") #this sets our "reference" for all groups to be 0h standard light uninfected

#runs the actual DeSeq analysis
dds <- DESeq(dds_matrix)

#extract the results names
res <- results(dds, name = 'group_infected_High_Light_4h_vs_control_High_Light_4h')|>
  as.data.frame()

res_no_NA <- res %>% 
  drop_na()

# How many rows did we filter out!?
glimpse(res_no_NA) 

res_filtered <- res_no_NA |>
  rownames_to_column("gene_id")|>
  mutate(up_down = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "UP", # label genes in the top right quadrant of the volcano plot
    padj < 0.05 & log2FoldChange < -1 ~ "DOWN", # label genes in the top left quadrant of the volcano plot
    TRUE ~ "NONE"
  ))

#get the top 10 differentially expressed genes
topGenes <- res_filtered %>%
  arrange(padj) %>%
  head(n = 10) %>%
  dplyr::select(gene_id) %>%
  pull()

GO_names_top <- GO_names|>
  filter(Geneid %in% topGenes)
  
  

dat_annotated <- res_filtered %>%
  mutate(annotation  = case_when(
    gene_id %in% topGenes ~ "GOI",
    TRUE ~ "normal"
  ))

volcano_h4h <- dat_annotated |>
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = annotation)) + # see that we are using the -log10() function to transform our padj values
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")+
  theme_classic() +
  geom_text_repel(data = filter(dat_annotated, annotation == "GOI"),
                  aes(label = gene_id),
                  colour = "black",
                  size = 3,
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  scale_color_manual(values = c("red", "grey"))+
  ggtitle("infected vs control, High Light 4 hours")

volcano_h4h

ggsave('plots/volcano/volcanoplot_h4h_infectedvscontrol.png', plot= volcano_h4h, height = 10)
