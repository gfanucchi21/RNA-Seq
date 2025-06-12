library(goseq)

# Import Deseq2 results
Deseq2_results <- read.csv("Deseq2_results.csv")

# Import the GO table for all genes in drosophila
go_all_genes_dm6 <- read.delim("C:/Users/user/WORK/Refs/go_all_genes_dm6.tsv", header=FALSE)

# These are the commands for setting the goseq parameters
supportedOrganisms()
supportedGenomes()
supportedGeneIDs()

# Set 1 for DE genes and 0 for not DE genes
Deseq2_results$DE <- 0
Deseq2_results$DE[Deseq2_results$FDR < 0.05] <- 1

# Filter the results for up and down regulated genes
Deseq2_up <- Deseq2_results %>% filter(log2FoldChange > 0)
Deseq2_down <- Deseq2_results %>% filter(log2FoldChange < 0)

#
## GO for up-regulated genes
#

# Create a named vector for up-regulated genes.
DEgenes_up <- Deseq2_up$DE
names(DEgenes_up) <- Deseq2_up$gene_id

# Create a vector for the gene lengths.
lenghts_up <- Deseq2_up$Length

# Calculate the probability Weighting function for the up-regulated genes.
pwf_up <- nullp(DEgenes_up, bias.data = lenghts_up)

# Perform the GO analysis on up-regulated genes.
pvals_up <- goseq(
  pwf_up,
  gene2cat = go_all_genes_dm6,
  test.cats = c("GO:CC", "GO:BP", "GO:MF"),
  method = "Wallenius",
  use_genes_without_cat = FALSE
)

# Compute the FDR.
pvals_up$FDR <- p.adjust(pvals_up$over_represented_pvalue, method="BH")

# Print the top 10 GO terms.
head(pvals_up)

#
# GO for down-regulated genes
#

# Create a named vector for down-regulated genes
DEgenes_down <- Deseq2_down$DE
names(DEgenes_down) <- Deseq2_down$gene_id

# Create a vector for the gene lenghts
lenghts_down <- Deseq2_down$Length

# Calculate the probability Weighting function for the down-regulated genes
pwf_down <- nullp(DEgenes_down, bias.data = lenghts_down)

# Perform the GO analysis on up-regulated genes
pvals_down <- goseq(
  pwf_down,
  gene2cat = go_all_genes_dm6,
  test.cats = c("GO:CC", "GO:BP", "GO:MF"),
  method = "Wallenius",
  use_genes_without_cat = FALSE
)


# Compute the FDR
pvals_down$FDR <- p.adjust(pvals_down$over_represented_pvalue, method="BH")

# Separate the three ontologies into three tables and save them
ontology <- c("CC", "BP", "MF")

p_vals_ontology <- lapply(ontology, function(x){
  pvals_filt <- pvals_down %>% filter(ontology == x)
  return(pvals_filt)
}
)

names(p_vals_ontology) <- ontology

lapply(seq_along(p_vals_ontology), function(i) {
  write.csv(p_vals_ontology[[i]], file = paste0("GO_", names(p_vals_ontology)[i], "_downreg.csv"), 
            row.names = FALSE, quote = FALSE)
})

#
### Bubble plots
#

## Bubble plot for down-regulated genes

# Select only the top 10 terms for each ontology and prepare them for plotting
p_vals_ontology_top10 <- lapply(p_vals_ontology, function(x){
  df <- head(x, n = 10)
  df$term <- factor(df$term, levels = df$term)
  df <- df[order(df$over_represented_pvalue), ]
  return(df)
})

library(ggplot2)
library(viridisLite)

# Create bubble plots
bubble_plots <- lapply(seq_along(p_vals_ontology_top10), function(i){
  plot <- ggplot(p_vals_ontology_top10[[i]], aes(x = (numDEInCat/numInCat)*100 , y = term, size = numDEInCat, color = -log10(over_represented_pvalue))) +
    geom_point(alpha = 0.8) +
    scale_size(range = c(3, 15), name = "Count") +
    scale_color_viridis_c(option = "D", name = "-Log10(P value)", direction = -1 ) +
    scale_y_discrete(limits = rev) +
    labs(
      title = paste0("Top over-represented categories in ", names(p_vals_ontology_top10)[i], " for down-regulated genes"),
      subtitle = "Wallenius method",
      x = "% DE in category",
      y = "Category"
    ) +
    theme_bw()
  })

# Plot for Cellular component ontology
pdf(file = "top10_go_downreg_CC.pdf", width = 8, height = 6) # you can change the size of the output file
# Execute the plot
bubble_plots[[1]]
# Close the file that will contain the plot
dev.off()

# Plot for Biological Processes ontology
pdf(file = "top10_go_downreg_BP.pdf", width = 10, height = 6) # you can change the size of the output file
# Execute the plot
bubble_plots[[2]]
# Close the file that will contain the plot
dev.off()

# Plot for molecular Function ontology
pdf(file = "top10_go_downreg_MF.pdf", width = 12, height = 6) # you can change the size of the output file
# Execute the plot
bubble_plots[[3]]
# Close the file that will contain the plot
dev.off()


## Bubble plot using FDR < 0.05
# Select only the terms with FDR < 0.05 for each ontology 
FDR_ontology_top <- lapply(p_vals_ontology, function(x){
  df <- filter(x, FDR < 0.05)
  df$term <- factor(df$term, levels = df$term)
  df <- df[order(df$over_represented_pvalue), ]
  return(df)
})

FDR_total <- rbind(FDR_ontology_top[["CC"]], FDR_ontology_top[["BP"]], FDR_ontology_top[["MF"]])



## Bubble plot for up-regulated genes
df <- pvals_up[1:10,]

df$term <- as.factor(df$term)

df <- df[order(df$over_represented_pvalue), ]

df$term <- factor(df$term, levels = df$term)

library(ggplot2)
library(viridisLite)

up_plot <- ggplot(df, aes(x = (numDEInCat/numInCat)*100 , y = term, size = numDEInCat, color = -log10(over_represented_pvalue))) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 15), name = "Count") +
  scale_color_viridis_c(option = "D", name = "-Log10(P value)", direction = -1 ) +
  scale_y_discrete(limits = rev) +
  xlim(45, 100) +
  labs(
    title = "Top over-represented categories in CC,BP,MF for up-regulated genes",
    subtitle = "Wallenius method",
    x = "% DE in category",
    y = "Category"
  ) +
  theme_bw()

pdf(file = "D12_Hisat2_Deseq2_top10_go_upreg_CC_BP_MF.pdf", width = 10, height = 6) # you can change the size of the output file
# Execute the plot
up_plot
# Close the file that will contain the plot
dev.off()

# Data from PANGEA studying the KEGG pathways
KEGG <- read.csv("PANGEA_KEGG_downreg_D12_Hisat2.csv")
KEGG <- KEGG %>% select(!X.Select.)
KEGG <- head(KEGG, n = 10)
KEGG$Gene.Set.Name <- factor(KEGG$Gene.Set.Name, levels = KEGG$Gene.Set.Name)
KEGG <- KEGG[order(KEGG$P.value), ]

KEGG_plot <- ggplot(KEGG, aes(x = (Count.Overlap.Gene / Gene.Set.Size)*100, y = Gene.Set.Name, size = Count.Overlap.Gene, color = P.value.calc.log_10)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 15), name = "Count") +
  scale_color_viridis_c(option = "D", name = "-Log10(P value)", direction = -1 ) +
  scale_y_discrete(limits = rev) +
  #xlim(45, 100) +
  labs(
    title = "Top over-represented categories in KEGG for down-regulated genes",
    x = "% DE in category",
    y = "Category"
  ) +
  theme_bw()

# Plot for Cellular component ontology
pdf(file = "D12_Hisat2_Deseq2_top10_go_downreg_KEGG_PANGEA.pdf", width = 8, height = 6) # you can change the size of the output file
# Execute the plot
KEGG_plot
# Close the file that will contain the plot
dev.off()
