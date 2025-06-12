# Import human genes symbols
human_genes <- read.table("C:/Users/gabri/WORK/Chip-seq_zzz3/zzz3-bound_downreg_genes.txt", quote="\"", comment.char="")
library(biomaRt)

# Connect to the Asia mirror
listMarts(host="https://asia.ensembl.org")

ENSEMBL_MART_ENSEMBL = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://asia.ensembl.org")
head(listDatasets(ENSEMBL_MART_ENSEMBL))

ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="https://asia.ensembl.org")

# The "listFilters" function will give you the list of available filters for a given mart and species:
head(listFilters(ensembl))

# The "listAttributes" function will give you the list of the available attributes for a given mart and species:
head(listAttributes(ensembl))

# Extract for each gene symbol, the respective ensemble gene id
human_gene_id <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = human_genes,
  mart = ensembl,
  uniqueRows = TRUE
)

# Save ensemble gene id in a vector
human_ensembl_id <- human_gene_id$ensembl_gene_id


### Find drosophila orthologs
# The use of the dec2021 archive is required to avoid errors when using getLDS()
listEnsemblArchives()
human.mart <- useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
dm.mart <-  useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="dmelanogaster_gene_ensembl")

# Find drosophila ortholgs from human enseble gene id
dm_orthologs <- getLDS(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = human_ensembl_id,
  mart = human.mart,
  attributesL = c("ensembl_gene_id", "external_gene_name"),
  martL = dm.mart # Specify drosophila
)

# save orthologs in a vector
dm_orthologs_name <- dm_orthologs$Gene.name.1


#
### Check which genes are common in both datasets
#

# Import genes down-regulated in YEATS2 kd in drosophila
Deseq2_results <- read.csv("D12_Deseq2_results.csv")
down_reg_symbols <- (Deseq2_results %>% filter((log2FoldChange < 0) & (FDR < 0.05)) %>% select("gene_symbol"))$gene_symbol

# Find common genes
common_genes <- intersect(dm_orthologs_name, down_reg_symbols)

# Run GO analysis on these genes
common_genes_GO <- getBM(
  attributes = c("external_gene_name", "go_id", "name_1006", "namespace_1003" ),
  filters = "external_gene_name",
  values = common_genes,
  mart = dm.mart,
  uniqueRows = TRUE
)

# Print the list of common genes
common_genes

# Save gene names in a file
write(common_genes, "genes.regulated.by.ATAC.in.tumor.and.dm.txt", sep="\t")
