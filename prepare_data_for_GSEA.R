Deseq2_results <- read.csv("D12_Deseq2_results.csv")

write.table(Deseq2_results$gene_id, "Deseq2_all_FB_gene_id.txt", quote = F, row.names = F, col.names = F, sep = "\t" )

# As the gene sets downloaded use the CG form of gene id, 
# the gene ids present in D12_Deseq2_results need to be converted online on https://www.biotools.fr/drosophila/fbgn_converter 
# using the Fbgn>>CG option and keeping original IDs in the output

Deseq2_all_CG_gene_id <- read.delim("Deseq2_all_CG_gene_id.txt", header=FALSE)

# The converted ids are added to D12_Deseq2_results using the script prepare_data_for_GSEA.R on R studio
Deseq2_results <- left_join(Deseq2_results, Deseq2_all_CG_gene_id, by = join_by(gene_id == V1))

Deseq2_results <- Deseq2_results %>% 
  filter(V2 != "") %>%
  dplyr::select(V2, bam.N1_md_filt_chr.bam, bam.N2_md_filt_chr.bam, bam.N3_md_filt_chr.bam,
         bam.D1_md_filt_chr.bam, bam.D2_md_filt_chr.bam, bam.D3_md_filt_chr.bam)

Deseq2_results <- Deseq2_results %>% mutate(description = "NA", .after = V2)

colnames(Deseq2_results) <- c("Name", "Description", "N1",  "N2", "N3", "D1", "D2", "D3")

df <- Deseq2_results[, c("Name", "Description", "D1", "D2", "D3", "N1", "N2", "N3")]

write.table(df, "Deseq2_D12_counts_for_GSEA.txt", quote = F, row.names = F, col.names = T, sep = "\t" )
