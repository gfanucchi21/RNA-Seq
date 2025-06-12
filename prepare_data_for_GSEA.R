# Import DEseq2 results
Deseq2_results <- read.csv("Deseq2_results.csv")

# Write a table containing all FlyBase IDs contained in the Deseq2_results.csv file
write.table(Deseq2_results$gene_id, "Deseq2_all_FB_gene_id.txt", quote = F, row.names = F, col.names = F, sep = "\t" )

# As the gene sets downloaded use the CG form of gene ID, 
# the gene ids present in Deseq2_results need to be converted online on https://www.biotools.fr/drosophila/fbgn_converter 
# using the Fbgn>>CG option and keeping original IDs in the output
# and saved in file named Deseq2_all_CG_gene_id.txt

Deseq2_all_CG_gene_id <- read.delim("Deseq2_all_CG_gene_id.txt", header=FALSE)

# The CG gene IDs are added to the Deseq2_results table
Deseq2_results <- left_join(Deseq2_results, Deseq2_all_CG_gene_id, by = join_by(gene_id == V1))

# Remove the rows that do not have a CG gene ID
Deseq2_results <- Deseq2_results %>% 
  filter(V2 != "") %>%
  dplyr::select(V2, bam.N1_md_filt_chr.bam, bam.N2_md_filt_chr.bam, bam.N3_md_filt_chr.bam,
         bam.D1_md_filt_chr.bam, bam.D2_md_filt_chr.bam, bam.D3_md_filt_chr.bam)

# Add a column named "description" filled with NAs. It is needed for the GSEA
Deseq2_results <- Deseq2_results %>% mutate(description = "NA", .after = V2)

# Set column names
colnames(Deseq2_results) <- c("Name", "Description", "N1",  "N2", "N3", "D1", "D2", "D3")

# Select only the needed columns
df <- Deseq2_results[, c("Name", "Description", "D1", "D2", "D3", "N1", "N2", "N3")]

# Save the table ready for GSEA
write.table(df, "Deseq2_D12_counts_for_GSEA.txt", quote = F, row.names = F, col.names = T, sep = "\t" )
