#!/usr/bin/env Rscript
# Differential expression analysis with the DESeq2 package.
# 
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#

# Load the library.
library(DESeq2)

# The name of the file that contains the counts.
counts_file = "featurecounts/counts.txt"

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = "design.csv"

# Specify the final result file name.
output_file = "Deseq2_results.csv"

# Specify the PCA plot file name.
PCA_file = "Deseq2_PCA.pdf"


# Read the sample file.
colData <- read.csv(design_file, stringsAsFactors=F)

# Turn conditions into factors.
colData$condition = factor(colData$condition)

# The first level should correspond to the first entry in the file!
# Required later when building a model.
colData$condition = relevel(colData$condition, toString(colData$condition[1]))

# Isolate the sample names.
sample_names <- colData$sample

# Read the data
df = read.delim(counts_file, header=TRUE, row.names=1, comment.char="#" )

# Remove the repetition of the chromosome and strand indication
# Keep only the start and end coordinate of the gene
df$Chr <- sapply(strsplit(df$Chr, ";"), function(x) x[1])
df$Start <- sapply(strsplit(df$Start, ";"), function(x) x[1])
df$End <- sapply(strsplit(df$End, ";"), function(x) tail(x, 1))
df$Strand <- sapply(strsplit(df$Strand, ";"), function(x) x[1])

# Set the smallest group size
smallestGroupSize <- 3

# Filter out all rows in which at least 3 samples have more or equal than 10 reads
keep <- rowSums(df[,7:12] >= 10) >= smallestGroupSize
df <- df[keep,]

# Created rounded integers for the count data
countData = round(df[, sample_names])

# Other columns in the dataframe that are not sample information. 
otherCols = df[!(names(df) %in% sample_names)]

#
# Running DESeq2
#

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

# Run DESeq
dse = DESeq(dds)

# Format the results.
res = results(dse)

# Create a dataset of counts to which the variance stabilizing transformation has been applied
# It is needed for the quality control visualization
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

# Heatmap of the sample-to-sample distance
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Principal component plot
plotPCA(vsd, intgroup= "condition")
PCA_plot <- plotPCA(vsd, intgroup= "condition")

# Save the PCA plot in a .pdf file
pdf(file = PCA_file, width = 6, height = 6) 
PCA_plot
dev.off()

#
# The rest of the code is about formatting the output dataframe.
#

# Turn the DESeq2 results into a data frame.
data = cbind(otherCols, data.frame(res))

# Create the foldChange column.
data$foldChange = 2 ^ data$log2FoldChange

# Rename columns to better reflect reality.
names(data)[names(data)=="pvalue"] <-"PValue"
names(data)[names(data)=="padj"] <- "FDR"

# Create a real adjusted p-value
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Sort the data by PValue to compute false discovery counts.
data = data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR

# Create the additional columns that we wish to present.
data$baseMeanA = 1
data$baseMeanB = 1

# Get the normalized counts.
normed = counts(dse, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

# Merge the two datasets by row names.
total <- merge(data, normed, by=0)

# Sort again for output.
total = total[with(total, order(PValue, -foldChange)), ]

# Sample names for condition A
col_names_A = data.frame(split(colData, colData$condition)[1])[,1]

# Sample names for condition B
col_names_B = data.frame(split(colData, colData$condition)[2])[,1]

# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])

# Bringing some sanity to numbers. Round columns to fewer digits.
total$foldChange = round(total$foldChange, 3)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
total$falsePos = round(total$falsePos, 0)

# Rename the first column.
colnames(total)[1] <- "gene_id"

# Reorganize columns names to make more sense.
new_cols = c("gene_id", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
             "log2FoldChange","lfcSE","stat","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]

library(dplyr)
diffexp <- total %>% filter(as.numeric(FDR) < 0.05)

n_up_exp <-  nrow(total %>% filter(as.numeric(FDR) < 0.05 & (log2FoldChange > 0)))
print(n_up_exp)

n_down_exp <-  nrow(total %>% filter(as.numeric(FDR) < 0.05 & (log2FoldChange < 0)))
print(n_down_exp)

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)
total$FDR = formatC(total$FDR, format = "e", digits = 1)

# Write the results
write.csv(total, file=output_file, row.names=FALSE, quote=FALSE)
