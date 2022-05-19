# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.14")
BiocManager::install("BiocParallel")

# load packages
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)

# load datasets
SRR6040092 <- read.delim("htseq_CDS_SRR6040092_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6040093 <- read.delim("htseq_CDS_SRR6040093_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6040094 <- read.delim("htseq_CDS_SRR6040094_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6040096 <- read.delim("htseq_CDS_SRR6040096_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6040097 <- read.delim("htseq_CDS_SRR6040097_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6156066 <- read.delim("htseq_CDS_SRR6156066_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6156067 <- read.delim("htseq_CDS_SRR6156067_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
SRR6156069 <- read.delim("htseq_CDS_SRR6156069_scaffold_06_Aligned.sortedByCoord.out.bam.txt", header = FALSE)
colnames(SRR6040092)=c('gene', 'SRR6040092')
colnames(SRR6040093)=c('gene', 'SRR6040093')
colnames(SRR6040094)=c('gene', 'SRR6040094')
colnames(SRR6040096)=c('gene', 'SRR6040096')
colnames(SRR6040097)=c('gene', 'SRR6040097')
colnames(SRR6156066)=c('gene', 'SRR6156066')
colnames(SRR6156067)=c('gene', 'SRR6156067')
colnames(SRR6156069)=c('gene', 'SRR6156069')

df <- merge(SRR6040092, SRR6040093, by = 'gene')
df <- merge(df, SRR6040094, by = 'gene')
df <- merge(df, SRR6040096, by = 'gene')
df <- merge(df, SRR6040097, by = 'gene')
df <- merge(df, SRR6156066, by = 'gene')
df <- merge(df, SRR6156067, by = 'gene')
df <- merge(df, SRR6156069, by = 'gene')
df <- df[-c(1:5),]

# create dataframes
df <- data.frame(row.names = df[,1],
                 leaf = df[2], root = df[3], aril_1 = df[4], stem = df[5],
                 aril_2 = df[6], aril_1_mn = df[7], aril_2_mn = df[8], aril_3_mn = df[9])

samples <- colnames(df)[2:9]
tissue_list <- c('leaf', 'root', 'aril_mk', 'stem', 'aril_mk',
            'aril_mon', 'aril_mon', 'aril_mon')
tissue <- data.frame(tissue_list)
rownames(tissue) <- samples

# run DESeq
DESdf <- DESeqDataSetFromMatrix(countData = df, colData = tissue, design = ~ tissue_list)
dds <- DESeq(DESdf)

# run PCA
rlog_dds <- rlog(dds)
plotPCA(rlog_dds, intgroup=c("tissue_list"))

# create heatmap
topVarGenes <- head(order(rowVars(assay(rlog_dds)), decreasing = TRUE), 10)

heatmap.2(assay(rlog_dds)[topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", margins=c(8,10), 
          col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

