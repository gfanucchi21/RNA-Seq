# RNA-Seq
Scripts for RNA-Seq analysis

## Project description
This project has been created to analyse the RNA-Seq data from Drosophila cells

## Analysis steps
- Sequencing reads are aligned to the Drosophila melanogaster genome using Hisat2
- Reads are filtered based on the alignment quality
- Statistical analysis is performed both before and after filtering
- BigWig coverage files are generated
- The quantification matrix is computed using FeatureCounts
- Differential expression analysis is performed using DESeq2
- PCA plot, heatmap, and volcano plots are generated
- Gene Ontology analysis is performed using the goseq package
- Results are prepared to perform GSEA

## Running instructions
master_script.sh contains all the bash commands to be executed through the Linux terminal or WSL.
Inside master_script.sh has been specified when to execute the R scripts (deseq2.r, heatmap.r, volcanoplot.r, goseq.r, prepare_data_for_GSEA.r) on R Studio.
