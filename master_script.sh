#!/bin/bash
# Generate ids file
conda activate bioinfo
cd /mnt/c/Users/user/WORK/RNAseq
parallel -j 1 echo {1}{2} ::: D N ::: 1 2 3 > ids

# Download dm genome from Flybase https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta/dmel-all-aligned-r6.63.fasta.gz
cd /mnt/c/Users/user/Downloads
gzip -d dmel-all-chromosome-r6.63.fasta.gz 
mv /mnt/c/Users/gabri/Downloads/dmel-all-aligned-r6.63.fasta /mnt/c/Users/gabri/WORK/Refs

### Align dm mRNA from Refseq on dm genome to be visualized on IGV
cd /mnt/c/Users/user/WORK/Refs
# The index to the genome
IDX=dmel-all-chromosome-r6.63.fasta
# Build the index
conda activate hisat2
hisat2-build --noauto --bmaxdivn 32 --dcv 4096 $IDX $IDX
# Index the reference genome with samtools.
conda activate bioinfo
samtools faidx $IDX
# Create the transcript alignment BAM file.
hisat2 -x $IDX -f -U mrna.fa -S mrna.sam
conda activate bioinfo
samtools sort mrna.sam > mrna.bam
# Index the BAM file
samtools index mrna.bam

### Basic statistics on the .fastq files
cd /mnt/c/Users/user/WORK/Refs
seqkit stats dmel-all-chromosome-r6.63.fasta

### Read Alignment using hisat2
conda activate hisat2
cd /mnt/c/Users/user/WORK/RNAseq
# Create bam folder
mkdir bam
IDX=/mnt/c/Users/user/WORK/Refs/dmel-all-chromosome-r6.63.fasta
# Align the FASTQ files to the reference genome. 
# parameters taken from https://www.biorxiv.org/content/10.1101/2024.04.10.588867v1.full#F1
cat ids | parallel "hisat2 --min-intronlen 40 --max-intronlen 200000 --rna-strandness RF -x $IDX -1 reads/{}_1.fastq -2 reads/{}_2.fastq -S bam/{}.sam"
# Sort each SAM file in a BAM file
conda activate bioinfo
cat ids | parallel -j 1 "samtools sort bam/{}.sam > bam/{}.bam"
# Index each BAM file.
cat ids | parallel -j 4 "samtools index bam/{}.bam"
# Mark duplicated reads
cat ids | parallel -j 1 "picard MarkDuplicates I=bam/{}.bam O=bam/{}_md.bam M=bam_stats/{}_md_metrics.txt"
# Index each BAM file.
cat ids | parallel -j 4 "samtools index bam/{}_md.bam"
# Run statistics on unfiltered bams
cat ids | parallel -j 4 "samtools stats bam/{}_md.bam > bam_stats/Unfiltered/{}_md.stats"
cat ids | parallel -j 4 "samtools flagstat bam/{}_md.bam > bam_stats/Unfiltered/{}_md.flagstats"
cat ids | parallel -j 4 "samtools idxstat bam/{}_md.bam > bam_stats/Unfiltered/{}_md.idxstats"
multiqc bam_stats/Unfiltered -f --title "unfiltered bam stats" -n Unfiltered/unfiltered_bam.multiqc

### Read Filtering
# Read filtering feeping only uniquely mapping reads and properly paired reads
cat ids | parallel -j 4 "samtools view -b -q 30 -f 3 bam/{}_md.bam >  bam/{}_md_filt.bam"
cat ids | parallel -j 6 "samtools index bam/{}_md_filt.bam"
# Select only the alignemts on the Chr 2 3 4 X and Y
cat ids | parallel -j 6 "samtools view -b -o bam/{}_md_filt_chr.bam bam/{}_md_filt.bam 2L 2R 3L 3R 4 X Y "
cat ids | parallel -j 6 "samtools index bam/{}_md_filt_chr.bam"
# Run statistics on filtered bams
cat ids | parallel -j 4 "samtools stats bam/{}_md_filt_chr.bam > bam_stats/Filtered/{}_md_filt.stats"
cat ids | parallel -j 4 "samtools flagstat bam/{}_md_filt_chr.bam > bam_stats/Filtered/{}_md_filt.flagstats"
cat ids | parallel -j 4 "samtools idxstat bam/{}_md_filt_chr.bam > bam_stats/Filtered/{}_md_filt.idxstats"
multiqc bam_stats/Filtered -f --title "filtered bam stats" -n Filtered/filtered_bam.multiqc

### Generate BigWig coverage files
mkdir bigwig
# Turn each BAM file into bedGraph coverage
cat ids | parallel -j 4 "bedtools genomecov -ibam bam/{}_md_filt_chr.bam -split -bg > bigwig/{}_md_filt.bg"
# Sort bedGraph files
for id in $(cat ids)
do 
sort -k1,1 -k2,2n -i bigwig/${id}_md_filt.bg > bigwig/${id}_md_filt_sorted.bg
done
# Convert each bedGraph coverage into bigWig coverage
conda activate bedgraphtobigwig
cat ids | parallel -j 4 "bedGraphToBigWig bigwig/{}_md_filt_sorted.bg $IDX.fai bigwig/{}_md_filt.bw"

### Quantification matrix
conda activate featurecounts
# Download genome annotation file from https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/gtf/dmel-all-r6.63.gtf.gz
# Run the featureCounts program to summarize the number of reads that overlap with features.
REF=/mnt/c/Users/user/WORK/Refs/dmel-all-r6.63.gtf
# -s 2 for reverse strandness -p paired-end
featureCounts -T 6 -p -a ${REF} -s 2 -t exon -g gene_id --extraAttributes gene_symbol -o featurecounts/counts.txt bam/N?_md_filt_chr.bam bam/D?_md_filt_chr.bam
# Run summary statistics on featureCounts
conda activate multiqc
multiqc ./featurecounts -f -n featurecounts/D12_featurecounts.multiqc


### Run Deseq2.r script on R studio
### Run Heatmap.r scripts on R studio
### Run Volcanoplot.r on R studio


### Gene Ontology analysis
# Download flybase current gene ontology https://current.geneontology.org/annotations/fb.gaf.gz
# Prepare the file to be used in goseq in R
tail -n +34 fb.gaf | cut -f 2,5 - > go_all_genes_dm6.tsv
## Run goseq.R on R studio
## Gene ontology for KEGG patways was performed on PANGEA using the list of downregulated genes


### Gene Set Enrichment Analysis (GSEA)
# Gene sets specific for drosophila for KEGG and REACTOME were downloaded from blob:https://github.com/77b0bdac-b5c3-4ee6-b576-4f82943ff1ec blob:https://github.com/6c4303c2-d0db-405b-8148-3c1b481440d7
# Expression data in prepared using the script prepare_data_for_GSEA.r on R studio
# The file Deseq2_D12_counts_for_GSEA.txt was used together with the downloaded gene set for GSEA analysis on GSEA 4.4.0
## Configuration:
# Collapse/remap -> No collapse
# Permutation type -> gene set
# Other options -> default
