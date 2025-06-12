# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Enter the DEseq2 output data
df_volcano <- read.csv("D12_Deseq2_results.csv")

# Create a new column "delabel" to de, that will contain the name of the top 40 differentially expressed genes (NA in case they are not)
df_volcano$delabel <- ifelse(df_volcano$gene_symbol %in% head(df_volcano[order(df_volcano$FDR), "gene_symbol"], 30), df_volcano$gene_symbol, NA)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df_volcano$diffexpressed <- "NO"

# if log2Foldchange > 0 and FDR < 0.05, set as "UP"
df_volcano$diffexpressed[df_volcano$log2FoldChange > 0 & df_volcano$FDR < 0.05] <- "UP"

# if log2Foldchange < 0 and FDR < 0.05, set as "DOWN"
df_volcano$diffexpressed[df_volcano$log2FoldChange < 0 & df_volcano$FDR < 0.05] <- "DOWN"

# Volcano plot 
all_volcanoplot <- ggplot(data = df_volcano, aes(x = log2FoldChange, y = -log10(FDR), 
                                                    col = diffexpressed, label = delabel)) +
  #geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + # the vertical line should be added if a FC threshold is needed
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color = '', #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  # to customise the breaks in the x axis
  ggtitle('D12ko vs non-targeting control') + 
  # Plot title
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  theme_bw()

# Open the file that will contain your plot (the name is up to you)
pdf(file = "D12_Hisat2_Deseq2_volcano_all_data.pdf", width = 8, height = 6) # you can change the size of the output file
# Execute the plot
all_volcanoplot
# Close the file that will contain the plot
dev.off()


#
# Volcano plot of genes that have a FDR lower than a treshold
#

# Set the FDR threshold
threshold <- 50

# Filter out all those genes that have an FDR higher than the treshold
df_volcano_filt <- df_volcano %>% filter(-log10(FDR) < threshold)

# Volcano plot
filt_volcanoplot <- ggplot(data = df_volcano_filt, aes(x = log2FoldChange, y = -log10(FDR), 
                                       col = diffexpressed, label = delabel)) +
  #geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + # the vertical line should be added if a FC threshold is needed
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-12, 12)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = '', 
       #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  scale_x_continuous(breaks = seq(-12, 12, 2)) + 
  # to customise the breaks in the x axis
  ggtitle('D12ko vs non-targeting control') + 
  # Plot title
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  theme_bw()

# Open the file that will contain your plot (the name is up to you)
pdf(file = "D12_Hisat2_Deseq2_volcano_cropped.pdf", width = 8, height = 6) # you can change the size of the output file
# Execute the plot
filt_volcanoplot
# Close the file that will contain the plot
dev.off()




