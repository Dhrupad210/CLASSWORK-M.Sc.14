
# Bioinformatics Workflows

This repository documents major bioinformatics workflows such as Microarray, RNA-Seq, Alternative Splicing, Genome Assembly, Variant Calling, QIIME2, Bismark, MAFFT, and ChIP-Seq.  
Each section includes placeholders where you can later add your own commands, R code, and notes.

---

## Table of Contents
1. [Microarray Data Analysis (Affymetrix)](#microarray-data-analysis-affymetrix)
2. [RNA-Seq Pipeline](#rna-seq-pipeline)
3. [Alternative Splicing Analysis](#alternative-splicing-analysis)
4. [Genome Assembly and Annotation](#genome-assembly-and-annotation)
5. [Variant Calling and Annotation](#variant-calling-and-annotation)
6. [QIIME2 Metagenomics Workflow](#qiime2-metagenomics-workflow)
7. [Bismark Methylation Analysis](#bismark-methylation-analysis)
8. [ChIP-Seq Peak Calling and Analysis](#chip-seq-peak-calling-and-analysis)

---

## MICROARRAY DATA ANALYSIS (Affymetrix)
```R
## Load necessary libraries
library(affy)       # Core package for Affymetrix data
library(limma)      # For differential expression analysis
library(ggplot2)    # For advanced plotting
library(dplyr)      # For data manipulation
library(ggrepel)    # For non-overlapping plot labels

# --- Step 2: Loading Raw .CEL File Data ---
# List all .CEL files in the working directory
cel_files <- list.celfiles(full.names = TRUE)

# Read raw probe intensity data into an AffyBatch object
raw_affy_data <- ReadAffy(filenames = cel_files)

# --- Step 3: Exploration of Raw Probe-Level Data (PM & MM) ---
# Extract Perfect Match (PM) and Mismatch (MM) probe intensities
pm_probes <- as.data.frame(pm(raw_affy_data))
mm_probes <- as.data.frame(mm(raw_affy_data))

# --- Visualization of Raw Data ---
# Set up plotting grid (1 row, 2 columns)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Histogram of raw PM intensities for the first sample (WT1)
hist(pm_probes$EWE11465_0001_WT1.CEL,
     main = "Raw PM Intensities (WT1)",
     xlab = "Intensity",
     col = "steelblue",
     breaks = 100)

# Histogram of log2-transformed PM intensities for WT1
hist(log2(pm_probes$EWE11465_0001_WT1.CEL),
     main = "Log2 Transformed PM (WT1)",
     xlab = "Log2 Intensity",
     col = "salmon",
     breaks = 100)

# Set up plotting grid for scatter plots
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Scatter plot of raw PM values (WT1 vs MUT1)
plot(pm_probes$EWE11465_0001_WT1.CEL ~ pm_probes$EWE11465_0004_drna1.CEL,
     main = "Raw PM: WT1 vs MUT1",
     xlab = "MUT1 Raw Intensity",
     ylab = "WT1 Raw Intensity",
     pch = 19,
     col = rgb(0.2, 0.4, 0.6, 0.3))

# Scatter plot of log2-transformed PM values (WT1 vs MUT1)
plot(log2(pm_probes$EWE11465_0001_WT1.CEL) ~ log2(pm_probes$EWE11465_0004_drna1.CEL),
     main = "Log2 PM: WT1 vs MUT1",
     xlab = "MUT1 Log2 Intensity",
     ylab = "WT1 Log2 Intensity",
     pch = 19,
     col = rgb(0.8, 0.3, 0.3, 0.3))

# Reset plotting window to a single panel
par(mfrow = c(1, 1))

# --- Step 4: Pre-Normalization Quality Control ---
# Boxplot of raw log-intensities for all samples
boxplot(raw_affy_data,
        main = "Boxplot of Raw Log-Intensities",
        ylab = "Log2 Intensity",
        col = "#69b3a2",
        las = 2)

# --- Step 5: Normalization ---
# Perform RMA normalization
normalized_eset <- rma(raw_affy_data)

# Extract the matrix of normalized expression values
normalized_exprs <- exprs(normalized_eset)

# Boxplot of RMA normalized intensities
boxplot(normalized_exprs,
        main = "Boxplot of RMA Normalized Intensities",
        ylab = "Normalized Log2 Intensity",
        col = "#404080",
        las = 2)

# --- Step 6: Post-Normalization QC (PCA and Density Plots) ---
# Perform PCA on transposed normalized data
pca_results <- prcomp(t(exprs(normalized_eset)), scale. = TRUE)

# Define experimental groups (ensure order matches CEL files)
exp_groups <- factor(c("WT", "WT", "WT", "MUT", "MUT", "MUT"))

# Create data frame for PCA plotting
pca_plot_data <- data.frame(
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  Group = exp_groups,
  Sample = colnames(normalized_eset)
)

# Calculate variance explained by PC1 and PC2
variance_explained <- summary(pca_results)$importance[2, 1:2]
pc1_var <- round(variance_explained[1] * 100, 1)
pc2_var <- round(variance_explained[2] * 100, 1)

# Generate PCA plot using ggplot2
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), box.padding = 0.5) +
  scale_color_manual(values = c("WT" = "#0072B2", "MUT" = "#D55E00")) +
  labs(
    title = "PCA of Microarray Samples",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Experimental Group",
    shape = "Experimental Group"
  ) +
  theme_bw() +
  coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5))

# --- Step 6b: Density Plots (Before vs After Normalization) ---
raw_exprs_matrix <- exprs(raw_affy_data)
normalized_exprs_matrix <- normalized_exprs
sample_names <- colnames(raw_exprs_matrix)

par(mfrow = c(2, 6), mar = c(4, 4, 3, 1))
for (i in 1:ncol(raw_exprs_matrix)) {
  plot(density(raw_exprs_matrix[, i]),
       main = paste(sample_names[i], "\n(Raw)"),
       xlab = "Log2 Intensity",
       col = "#D55E00",
       lwd = 2,
       cex.main = 0.9)
}
for (i in 1:ncol(normalized_exprs_matrix)) {
  plot(density(normalized_exprs_matrix[, i]),
       main = paste(sample_names[i], "\n(Normalized)"),
       xlab = "Log2 Intensity",
       col = "#0072B2",
       lwd = 2,
       cex.main = 0.9)
}
par(mfrow = c(1, 1))

# --- Step 7: Differential Expression Analysis with limma ---
design_matrix <- model.matrix(~0 + exp_groups)
colnames(design_matrix) <- c("WT", "MUT")
contrast_matrix <- makeContrasts(MUT - WT, levels = design_matrix)
fit_model <- lmFit(normalized_eset, design_matrix)
fit_contrast <- contrasts.fit(fit_model, contrast_matrix)
fit_bayes <- eBayes(fit_contrast)
de_results <- topTable(fit_bayes, number = Inf, adjust.method = "BH")

# --- Step 8: Visualization of DE Results ---
results_df <- de_results %>%
  mutate(
    significance_level = case_when(
      logFC > 1 & adj.P.Val < 0.05 ~ "Upregulated",
      logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    gene_label = ifelse(adj.P.Val < head(sort(adj.P.Val), 10), rownames(.), NA)
  )

# Volcano plot
ggplot(results_df, aes(x = logFC, y = -log10(adj.P.Val), color = significance_level)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  labs(title = "Volcano Plot: MUT vs. WT",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Regulation") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# MA plot
ggplot(results_df, aes(x = AveExpr, y = logFC, color = significance_level)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#e41a1c", "Downregulated" = "#377eb8", "Not Significant" = "grey")) +
  geom_hline(yintercept = 0, color = "blue", linetype = "solid") +
  labs(
    title = "MA Plot: MUT vs. WT",
    x = "Average Log2 Expression",
    y = "Log2 Fold Change",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

# --- Step 9: Save Final Results ---
write.csv(de_results, file = "differential_expression_results.csv") [Add your shell or R commands here]
```

---

## RNA-SEQ PIPELINE
```bash
# --- Step 1: Initial Quality Control with FastQC ---
echo "Running initial FastQC..."
/usr/lib/jvm/java-17-openjdk-amd64/bin/java -jar /usr/share/java/fastqc.jar \
conA_rep1.fq conA_rep2.fq conB_rep1.fq conB_rep2.fq

# Comment:
# FastQC detects strong base composition bias in the first 12 bases across all reads.
# The first 12 bases should be trimmed to remove bias.

# --- Step 2: Trim Reads with Trimmomatic ---
echo "Trimming first 12 bases from all reads..."
java -jar /usr/share/java/trimmomatic-0.39.jar SE -phred33 conA_rep1.fq conA_rep1_trimmed.fq HEADCROP:12
java -jar /usr/share/java/trimmomatic-0.39.jar SE -phred33 conA_rep2.fq conA_rep2_trimmed.fq HEADCROP:12
java -jar /usr/share/java/trimmomatic-0.39.jar SE -phred33 conB_rep1.fq conB_rep1_trimmed.fq HEADCROP:12
java -jar /usr/share/java/trimmomatic-0.39.jar SE -phred33 conB_rep2.fq conB_rep2_trimmed.fq HEADCROP:12

# Comment:
# After trimming, FastQC shows uniform per-base content distribution.

# --- Step 3: Verify Trimmed Read Quality ---
echo "Verifying post-trim quality..."
/usr/lib/jvm/java-17-openjdk-amd64/bin/java -jar /usr/share/java/fastqc.jar \
conA_rep1_trimmed.fq conA_rep2_trimmed.fq conB_rep1_trimmed.fq conB_rep2_trimmed.fq

# --- Step 4: Build Genome Index for Bowtie2 ---
mkdir bowtie2_index
bowtie2 --version
bowtie2-build "Copy of GCF_000146045.2_R64_genomic.fna" bowtie2_index/genome_index

# --- Step 5: Align Trimmed Reads to the Genome ---
bowtie2 -x bowtie2_index/genome_index -U conA_rep1_trimmed.fq -S conA_rep1.sam
bowtie2 -x bowtie2_index/genome_index -U conA_rep2_trimmed.fq -S conA_rep2.sam
bowtie2 -x bowtie2_index/genome_index -U conB_rep1_trimmed.fq -S conB_rep1.sam
bowtie2 -x bowtie2_index/genome_index -U conB_rep2_trimmed.fq -S conB_rep2.sam

# --- Step 6: Convert, Sort, and Index BAM Files ---
samtools view -bS conA_rep1.sam > conA_rep1.bam
samtools sort conA_rep1.bam -o conA_rep1_sorted.bam
samtools index conA_rep1_sorted.bam

samtools view -bS conA_rep2.sam > conA_rep2.bam
samtools sort conA_rep2.bam -o conA_rep2_sorted.bam
samtools index conA_rep2_sorted.bam

samtools view -bS conB_rep1.sam > conB_rep1.bam
samtools sort conB_rep1.bam -o conB_rep1_sorted.bam
samtools index conB_rep1_sorted.bam

samtools view -bS conB_rep2.sam > conB_rep2.bam
samtools sort conB_rep2.bam -o conB_rep2_sorted.bam
samtools index conB_rep2_sorted.bam

# --- Step 7: Gene Count Quantification ---
sudo apt install subread
featureCounts -a "Copy of GCF_000146045.2_R64_genomic.gtf" -o gene_counts.txt \
  conA_rep1_sorted.bam conA_rep2_sorted.bam conB_rep1_sorted.bam conB_rep2_sorted.bam
```

```r
# Load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(GEOQuery)
library(ggrepel)

# Import gene count data
counts_data <- read.delim("gene_counts.txt", skip = 1, row.names = 1)
counts_data <- counts_data[, 6:ncol(counts_data)]
colnames(counts_data) <- c("conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2")

# Create metadata for samples
col_data <- data.frame(
  sample = colnames(counts_data),
  condition = factor(c("conA", "conA", "conB", "conB"))
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~ condition
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
results_table <- results(dds)
results_ordered <- results_table[order(results_table$padj), ]
write.csv(as.data.frame(results_ordered), file = "deseq2_results.csv")

# Filter significant genes
significant_genes <- subset(results_table, padj < 0.05 & abs(log2FoldChange) > 1)

# Variance stabilizing transformation for visualization
vst_data <- vst(dds, blind = FALSE)
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

# Heatmap of significant genes
pheatmap(assay(vst_data)[rownames(significant_genes), ],
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_col = annotation_col,
         show_rownames = FALSE)

# Top 20 DE genes heatmap
top_genes <- significant_genes[order(significant_genes$padj), ]
top_20_gene_names <- rownames(head(top_genes, 20))
pheatmap(assay(vst_data)[top_20_gene_names, ],
         main = "Top 20 Differentially Expressed Genes",
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_col = annotation_col,
         show_rownames = TRUE)

# --- Volcano Plot Visualization ---
results_df <- na.omit(as.data.frame(results_table))
results_df$diffexpressed <- "NO"
results_df$diffexpressed[results_df$log2FoldChange > 1.5 & results_df$padj < 0.05] <- "UP"
results_df$diffexpressed[results_df$log2FoldChange < -1.5 & results_df$padj < 0.05] <- "DOWN"

# Identify top genes for labeling
top_genes <- bind_rows(
  results_df %>% filter(diffexpressed == "UP") %>% arrange(padj) %>% head(20),
  results_df %>% filter(diffexpressed == "DOWN") %>% arrange(padj) %>% head(20)
)
results_df$delabel <- NA
results_df$delabel[rownames(results_df) %in% rownames(top_genes)] <-
  rownames(results_df)[rownames(results_df) %in% rownames(top_genes)]

# Volcano plot
ggplot(data = results_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point() +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")

# Volcano plot with gene labels
ggplot(data = results_df,
       aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_point() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot with Gene Labels",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
```

---

## ALTERNATIVE SPLICING ANALYSIS
```bash
# --- Step 1: Build HISAT2 Genome Index ---
hisat2-build 'Copy of GCF_000146045.2_R64_genomic.fna' genome_index

# --- Step 2: Align Reads to Reference Genome ---
hisat2 -x genome_index -U conA_rep1.fq -S conA_rep1.sam
hisat2 -x genome_index -U conA_rep2.fq -S conA_rep2.sam
hisat2 -x genome_index -U conB_rep1.fq -S conB_rep1.sam
hisat2 -x genome_index -U conB_rep2.fq -S conB_rep2.sam

# --- Step 3: Convert SAM to BAM ---
samtools view -bS conA_rep1.sam > conA_rep1.bam
samtools view -bS conA_rep2.sam > conA_rep2.bam
samtools view -bS conB_rep1.sam > conB_rep1.bam
samtools view -bS conB_rep2.sam > conB_rep2.bam

# --- Step 4: Sort BAM Files ---
samtools sort conA_rep1.bam -o conA_rep1_sorted.bam
samtools sort conA_rep2.bam -o conA_rep2_sorted.bam
samtools sort conB_rep1.bam -o conB_rep1_sorted.bam
samtools sort conB_rep2.bam -o conB_rep2_sorted.bam

# --- Step 5: Index Sorted BAM Files ---
samtools index conA_rep1_sorted.bam
samtools index conA_rep2_sorted.bam
samtools index conB_rep1_sorted.bam
samtools index conB_rep2_sorted.bam

# --- Step 6: Generate Gene Counts with featureCounts ---
featureCounts -a 'Copy of GCF_000146045.2_R64_genomic.gtf' -o hisat_gene_counts.txt \
  conA_rep1_sorted.bam conA_rep2_sorted.bam \
  conB_rep1_sorted.bam conB_rep2_sorted.bam
```

```r
# --- Step 1: Load Required Libraries ---
library(SGSeq)
library(pheatmap)

# --- Step 2: Define BAM Files and Sample Information ---
file_bam <- c("conA_rep1_sorted.bam", "conA_rep2_sorted.bam",
              "conB_rep1_sorted.bam", "conB_rep2_sorted.bam")
sample_name <- c("conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2")
bam_info <- data.frame(sample_name, file_bam)

# Retrieve BAM metadata
Baminfo <- getBamInfo(bam_info)
Baminfo

# --- Step 3: Import and Convert Transcript Annotations ---
tx <- importTranscripts("Copy of GCF_000146045.2_R64_genomic.gtf")
TxFeat <- convertToTxFeatures(tx)

# Analyze splicing features from BAM alignments
sgfc <- analyzeFeatures(Baminfo, features = TxFeat)
sgfc_annot <- annotate(sgfc, TxFeat)

# --- Step 4: Inspect Annotated Features ---
colData(sgfc_annot)
rowRanges(sgfc_annot)
counts(sgfc_annot)[1:5, ]
head(FPKM(sgfc_annot))

# Save annotated splicing feature counts
write.csv(counts(sgfc_annot), file = "sgfc_annot_counts.csv", row.names = TRUE)

# --- Step 5: Visualize Splicing for Selected Gene ---
plotFeatures(sgfc_annot,
             geneID = unique(rowRanges(sgfc_annot)$geneID)[1],
             color_novel = "red")

plotFeatures(sgfc_annot,
             geneID = unique(rowRanges(sgfc_annot)$geneID)[1],
             color_novel = "green")

# --- Step 6: Generate Heatmap of Top Splicing Features ---
expr_mat <- assay(sgfc_annot)
top_feats <- head(order(rowMeans(expr_mat), decreasing = TRUE), 50)

pheatmap(expr_mat[top_feats, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Top Splicing Features Across Samples")
```

---

## GENOME ASSEMBLY AND ANNOTATION
```bash
# --- Step 1: Genome Assembly with Velvet ---

# Initialize Velvet assembly directory with k-mer size 31 using paired-end reads
velveth Velvet_output_Contig 31 -shortPaired -separate -fastq SR1.fastq SR2.fastq

# Build the assembly graph; discard contigs with coverage <5
velvetg Velvet_output_Contig -cov_cutoff 5


# --- Step 2: Assembly Quality Assessment with QUAST ---

# Evaluate assembly quality (N50, total length, GC content, etc.)
quast.py contigs.fa -o quast_out


# --- Step 3: Gene Prediction with AUGUSTUS ---

# Predict genes on assembled contigs using E. coli K12 model; output in GFF format
augustus --species=E_coli_K12 contigs.fa > Augustus_out/contig_1.gff

# Extract predicted gene sequences (nucleotide FASTA) from GFF
gffread -g Velvet_output_Contig/contigs.fa -x Augustus_out/output_gene.fa Augustus_out/contig_1.gff

# Extract predicted protein sequences from GFF
gffread -g Velvet_output_Contig/contigs.fa -y Augustus_out/output_protein.fa Augustus_out/contig_1.gff


# --- Step 4: Assembly Completeness Check with BUSCO ---

# Assess completeness using bacterial ortholog dataset
busco -i Velvet_output_Contig/contigs.fa -m genome -l bacteria_odb10 -o busco_bacteria_results --cpu 4

# Example BUSCO output:
# C:95.1%[S:91.9%,D:3.2%],F:4.8%,M:0.1%,n:124


# --- Step 5: Reference-Guided Scaffolding with RagTag ---

# Correct contigs using a reference genome
ragtag.py correct GCA_014131755.1_ASM1413175v1_genomic.fna contigs.fa

# Scaffold corrected contigs against the reference genome
ragtag.py scaffold -u GCA_014131755.1_ASM1413175v1_genomic.fna ragtag_output/ragtag.correct.fasta

# Assess quality of scaffolded assembly
quast.py ragtag_output/ragtag.scaffold.fasta -o quast_ragtag_out


# --- Step 6: BUSCO on Corrected Assembly ---

# Evaluate completeness of corrected (ragtag) assembly
busco -i ragtag_output/ragtag.correct.fasta -o busco_ragtag_correct -m genome --cpu 4

# Example interpretation:
# C:95.1%[S:91.9%,D:3.2%],F:4.8%,M:0.1%,n:124
```
---

## VARIANT CALLING AND ANNOTATION
```bash
# --- Step 1: Prepare Reference Genome ---
bwa index GCF_000146045.2_R64_genomic.fna
samtools faidx GCF_000146045.2_R64_genomic.fna
gatk CreateSequenceDictionary -R GCF_000146045.2_R64_genomic.fna


# --- Step 2: Align Reads with BWA ---
bwa mem GCF_000146045.2_R64_genomic.fna conA_rep1.fq > conA_rep1.sam
bwa mem GCF_000146045.2_R64_genomic.fna conA_rep2.fq > conA_rep2.sam
bwa mem GCF_000146045.2_R64_genomic.fna conB_rep1.fq > conB_rep1.sam
bwa mem GCF_000146045.2_R64_genomic.fna conB_rep2.fq > conB_rep2.sam


# --- Step 3: Convert and Sort SAM Files ---
samtools view -bS conA_rep1.sam | samtools sort -o conA_rep1_sorted.bam
samtools view -bS conA_rep2.sam | samtools sort -o conA_rep2_sorted.bam
samtools view -bS conB_rep1.sam | samtools sort -o conB_rep1_sorted.bam
samtools view -bS conB_rep2.sam | samtools sort -o conB_rep2_sorted.bam


# --- Step 4: Merge Replicates ---
samtools merge conA_merged.bam conA_rep1_sorted.bam conA_rep2_sorted.bam
samtools merge conB_merged.bam conB_rep1_sorted.bam conB_rep2_sorted.bam


# --- Step 5: Mark Duplicates ---
java -jar picard.jar MarkDuplicates \
  I=conA_merged.bam O=conA_merged_markdup.bam M=conA_markdup_metrics.txt
samtools index conA_merged_markdup.bam

java -jar picard.jar MarkDuplicates \
  I=conB_merged.bam O=conB_merged_markdup.bam M=conB_markdup_metrics.txt
samtools index conB_merged_markdup.bam


# --- Step 6: Add Read Groups ---
java -jar picard.jar AddOrReplaceReadGroups \
  I=conA_merged_markdup.bam O=conA_merged_markdup_rg.bam \
  RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=conA
samtools index conA_merged_markdup_rg.bam

java -jar picard.jar AddOrReplaceReadGroups \
  I=conB_merged_markdup.bam O=conB_merged_markdup_rg.bam \
  RGID=2 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=conB
samtools index conB_merged_markdup_rg.bam


# --- Step 7: Variant Calling with GATK HaplotypeCaller ---
gatk HaplotypeCaller -R GCF_000146045.2_R64_genomic.fna \
  -I conA_merged_markdup_rg.bam -O conA_variants.vcf
gatk HaplotypeCaller -R GCF_000146045.2_R64_genomic.fna \
  -I conB_merged_markdup_rg.bam -O conB_variants.vcf


# --- Step 8: Base Quality Score Recalibration (BQSR) ---
gatk BaseRecalibrator -R GCF_000146045.2_R64_genomic.fna \
  -I conA_merged_markdup_rg.bam --known-sites conA_variants.vcf -O conA_BQSR.table
gatk ApplyBQSR -R GCF_000146045.2_R64_genomic.fna \
  -I conA_merged_markdup_rg.bam --bqsr-recal-file conA_BQSR.table -O conA_merged_markdup_rg_bqsr.bam

gatk BaseRecalibrator -R GCF_000146045.2_R64_genomic.fna \
  -I conB_merged_markdup_rg.bam --known-sites conB_variants.vcf -O conB_BQSR.table
gatk ApplyBQSR -R GCF_000146045.2_R64_genomic.fna \
  -I conB_merged_markdup_rg.bam --bqsr-recal-file conB_BQSR.table -O conB_merged_markdup_rg_bqsr.bam


# --- Step 9: Final Variant Calling ---
gatk HaplotypeCaller -R GCF_000146045.2_R64_genomic.fna \
  -I conA_merged_markdup_rg_bqsr.bam -O conA_final_variants.vcf
gatk HaplotypeCaller -R GCF_000146045.2_R64_genomic.fna \
  -I conB_merged_markdup_rg_bqsr.bam -O conB_final_variants.vcf


# --- Step 10: Extract and Filter SNPs ---
gatk SelectVariants -R GCF_000146045.2_R64_genomic.fna \
  -V conA_final_variants.vcf --select-type-to-include SNP -O conA_snps.vcf
gatk VariantFiltration -R GCF_000146045.2_R64_genomic.fna \
  -V conA_snps.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
  -O conA_snps_filtered.vcf

gatk SelectVariants -R GCF_000146045.2_R64_genomic.fna \
  -V conB_final_variants.vcf --select-type-to-include SNP -O conB_snps.vcf
gatk VariantFiltration -R GCF_000146045.2_R64_genomic.fna \
  -V conB_snps.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
  -O conB_snps_filtered.vcf


# --- Step 11: Extract and Filter INDELs ---
gatk SelectVariants -R GCF_000146045.2_R64_genomic.fna \
  -V conA_final_variants.vcf --select-type-to-include INDEL -O conA_indels.vcf
gatk VariantFiltration -R GCF_000146045.2_R64_genomic.fna \
  -V conA_indels.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 200.0" \
  -O conA_indels_filtered.vcf

gatk SelectVariants -R GCF_000146045.2_R64_genomic.fna \
  -V conB_final_variants.vcf --select-type-to-include INDEL -O conB_indels.vcf
gatk VariantFiltration -R GCF_000146045.2_R64_genomic.fna \
  -V conB_indels.vcf --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 200.0" \
  -O conB_indels_filtered.vcf


# --- Step 12: Rename Chromosomes for Ensembl VEP ---
cat > rename.txt << 'EOF'
NC_001133.9 chrI
NC_001134.8 chrII
NC_001135.5 chrIII
NC_001136.10 chrIV
NC_001137.3 chrV
NC_001138.5 chrVI
NC_001139.9 chrVII
NC_001140.6 chrVIII
NC_001141.2 chrIX
NC_001142.9 chrX
NC_001143.9 chrXI
NC_001144.5 chrXII
NC_001145.3 chrXIII
NC_001146.8 chrXIV
NC_001147.6 chrXV
NC_001148.4 chrXVI
NC_001224.1 chrM
EOF

bcftools annotate --rename-chrs rename.txt conA_variants_final.filtered.PASS.vcf -Ov -o conA_variants_renamed.vcf

# Use Ensembl VEP for functional annotation (select Saccharomyces cerevisiae)
```

---

## QIIME2 METAGENOMICS WORKFLOW
```bash
QIIME2 MICROBIOME ANALYSIS PIPELINE
# --- Step 1: Fix Manifest File Format ---
# Replace commas with tabs to ensure compatibility with QIIME2 manifest format
sed 's/,/\t/' manifest.tsv > manifest_fixed.tsv


# --- Step 2: Import FASTQ Data into QIIME2 ---
# Import single-end FASTQ files (Phred33 encoding) using manifest file
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_fixed.tsv \
  --output-path demux-single.qza \
  --input-format SingleEndFastqManifestPhred33V2


# --- Step 3: Summarize Demultiplexed Data ---
# Generate summary visualization of sequence quality and per-sample read counts
qiime demux summarize \
  --i-data demux-single.qza \
  --o-visualization demux-single.qzv


# --- Step 4: Denoising and Chimera Removal with DADA2 ---
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-single.qza \
  --p-trim-left 0 \
  --p-trunc-len 150 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza


# --- Step 5: Summarize Feature Table ---
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv


# --- Step 6: Taxonomic Classification ---
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


# --- Step 7: Phylogenetic Tree Construction ---
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


# --- Step 8: Core Diversity Metrics (Alpha and Beta Diversity) ---
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 100 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results


# --- Step 9: Alpha Diversity Significance Tests ---

# 1. Faith’s Phylogenetic Diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# 2. Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv


# --- Step 10: Beta Diversity Significance Tests ---

# 1. Unweighted UniFrac by body site
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column body-site \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

# 2. Unweighted UniFrac by subject
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column subject \
  --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv \
  --p-pairwise


# --- Step 11: PCoA Visualization (Emperor Plots) ---

# 1. Unweighted UniFrac PCoA
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

# 2. Bray-Curtis PCoA
qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv


# --- Step 12: Alpha Rarefaction Curves ---

# (First attempt failed because max_depth > max sample frequency)
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 4000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

# Corrected version with adjusted max depth (1800)
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1800 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv


# --- Step 13: Taxonomic Visualization ---

# View taxonomy metadata
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Create taxonomic bar plots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv


# --- Step 14: Filter Samples by Metadata ---

# Example: Keep only samples where body-site = "gut"
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[body-site]='gut'" \
  --o-filtered-table gut-table.qza
```

---

## BISMARK METHYLATION ANALYSIS
```bash
# --- Step 1: Prepare Reference Genome ---

# Create a directory for the reference genome
mkdir reference_genome1

# Move the reference genome FASTA file into the genome directory
mv ref.fa reference_genome1

# Prepare the reference genome for Bismark (builds bisulfite-converted indexes)
bismark_genome_preparation reference_genome1/


# --- Step 2: Align Paired-End Reads ---

# Align paired-end reads to the bisulfite-prepared reference genome
# -1 and -2 specify the paired FASTQ files
# -o specifies the output directory for alignment results
bismark --genome reference_genome1/ -1 R1.fastq.gz -2 R2.fastq.gz -o bis_output


# --- Step 3: Extract Methylation Information ---

# Extract methylation calls from the aligned BAM file
# --paired-end: indicates paired-end input reads
# --no_overlap: prevents double-counting overlapping read pairs
# --bedGraph: generates bedGraph files for genome browser visualization
# --cytosine_report: generates genome-wide cytosine methylation report
# --genome_folder: specifies reference genome location
bismark_methylation_extractor \
  --paired-end \
  --no_overlap \
  --bedGraph \
  --cytosine_report \
  --genome_folder reference_genome1/ \
  bis_output/R1_bismark_bt2_pe.bam
```
---

## CHIP-SEQ PEAK CALLING AND ANALYSIS
```bash
# --- Step 1: Quality Control and Trimming with fastp ---

# Trim low-quality bases, remove adapters, and generate QC reports
fastp \
  -i SO_4933_D_CHR1_R1.fastq.gz \
  -o SO_4933_D_CHR1_R1_clean.fastq.gz \
  -h SO_4933_D_CHR1_R1_fastp.html \
  -j SO_4933_D_CHR1_R1_fastp.json \
  -q 20 \
  -l 50 \
  --detect_adapter_for_pe false \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --disable_duplication_filter false \
  --thread 4

fastp \
  -i SO_4933_D_INR1_R1.fastq.gz \
  -o SO_4933_D_INR1_R1_clean.fastq.gz \
  -h SO_4933_D_INR1_R1_fastp.html \
  -j SO_4933_D_INR1_R1_fastp.json \
  -q 20 \
  -l 50 \
  --detect_adapter_for_pe false \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --disable_duplication_filter false \
  --thread 4


# --- Step 2: Index Reference Genome for Alignment ---
bwa index E.coli_BW25113.fasta


# --- Step 3: Align Reads to Reference Genome ---
bwa mem E.coli_BW25113.fasta SO_4933_D_CHR1_R1_clean.fastq.gz > SO_4933_D_CHR1_R1.sam
bwa mem E.coli_BW25113.fasta SO_4933_D_INR1_R1_clean.fastq.gz > SO_4933_D_INR1_R1.sam


# --- Step 4: Convert and Sort SAM Files ---
samtools view -Sb SO_4933_D_CHR1_R1.sam | samtools sort -o SO_4933_D_CHR1_R1_sorted.bam
samtools view -Sb SO_4933_D_INR1_R1.sam | samtools sort -o SO_4933_D_INR1_R1_sorted.bam


# --- Step 5: Index Sorted BAM Files ---
samtools index SO_4933_D_CHR1_R1_sorted.bam
samtools index SO_4933_D_INR1_R1_sorted.bam


# --- Step 6: Peak Calling using MACS3 ---

# For the first dataset
macs3 callpeak -t SO_4933_D_CHR1_R1_sorted.bam \
  -f BAM \
  -g 4.6e6 \
  -n D_CHR1_peaks \
  --outdir macs3_output \
  --keep-dup all

# For the second dataset
macs3 callpeak -t SO_4933_D_INR1_R1_sorted.bam \
  -f BAM \
  -g 4.6e6 \
  -n D_INR1_peaks \
  --outdir macs3_output \
  --keep-dup all \
  --qvalue 0.01 \
  --bdg \
  --nomodel \
  --extsize 150


# --- Step 7: Copy NarrowPeak Files for Merging ---
cp D_CHR1_peaks_peaks.narrowPeak D_CHR1.narrowPeak
cp D_INR1_peaks_peaks.narrowPeak D_INR1.narrowPeak


# --- Step 8: Merge Peaks using HOMER ---
mergePeaks D_CHR1.narrowPeak D_INR1.narrowPeak -d 100 > D_merged_HOMER.txt


# --- Step 9: Merge and Sort Peaks Manually ---
cat D_CHR1.narrowPeak D_INR1.narrowPeak | sort -k1,1 -k2,2n > D_merged_peaks.bed


# --- Step 10: Merge Nearby Peaks and Average Signal ---
bedtools merge \
  -i <(cat D_CHR1.narrowPeak D_INR1.narrowPeak | sort -k1,1 -k2,2n) \
  -c 5,7,8 -o mean,mean,mean > D_merged_peaks_with_signal.bed


# --- Step 11: Find Closest Genes to Peaks ---
bedtools closest \
  -a D_merged_peaks_with_signal.bed \
  -b E.coli_BW25113_genes.bed \
  -d > D_merged_peaks_closest_gene.bed


# --- Step 12: Filter Peaks Within 500 bp of Genes ---
awk '$15 <= 500' D_merged_peaks_closest_gene.bed > D_merged_peaks_near_genes.bed


# --- Step 13: Count Peaks per Gene ---
cut -f14 D_merged_peaks_closest_gene.bed | sort | uniq -c | sort -nr > peak_counts_per_gene.txt


# --- Step 14: Simplify BED File for Visualization ---
cut -f1-3,4,14,15 D_merged_peaks_closest_gene.bed > D_merged_peaks_simple.bed


# --- Step 15: Extract Sequences for Motif Discovery ---
bedtools getfasta \
  -fi E.coli_BW25113.fasta \
  -bed D_merged_peaks_simple.bed \
  -fo D_merged_peaks.fa


# --- Step 16: Motif Discovery using MEME ---
meme D_merged_peaks.fa \
  -dna \
  -oc meme_out \
  -mod zoops \
  -nmotifs 5 \
  -minw 6 \
  -maxw 20 \
  -revcomp
```
### The End of Analysis !!
