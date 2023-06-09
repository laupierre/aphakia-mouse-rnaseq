## 220120_A00558_0160_AH5V77DSX3

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub ("star.IIT_RNM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)


annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]


torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
colnames (a) <- gsub ("star.IIT_RNM_", "", colnames (a))
a <- a[ ,-1]


pheno <- read.delim ("Moratalla phenodata 220303.txt")
pheno <- pheno[pheno$area == "SN", ]
pheno <- pheno[pheno$experiment ==1, ]
idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]

stopifnot (colnames (a) == pheno$sample)


dds <- DESeqDataSetFromMatrix(countData = round (a), colData = pheno, design = ~ genotype)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds

dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "AK", "WT"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Pitx3", ] 


## Sample to sample correlation, see https://rockefelleruniversity.github.io/RU_RNAseq/presentations/singlepage/RU_RNAseq_p3.html

sampleDists <- as.dist(1 - cor(log2 (counts(dds,normalized=TRUE)+1), method="pearson"))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Distance between samples plot.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors, annotation_row=df)
dev.off()


## heatmap plot
# Th, DDC, NR4A2, SLC6A3 (known target), SNCA, Pitx3, SLC18A2 (known target, see doi: 10.1111/j.1471-4159.2009.06404.x), Tcf7l2
select <- c("ENSMUSG00000000214.12", "ENSMUSG00000020182.17", "ENSMUSG00000026826.14", "ENSMUSG00000021609.7", "ENSMUSG00000025889.14", "ENSMUSG00000025229.16", "ENSMUSG00000025094.9", "ENSMUSG00000024985.21")

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Heatmap plot.pdf")
pheatmap( log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select, ], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off ()


## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, label=sample)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		geom_text_repel()  + 
		  coord_fixed () 

ggsave ("PCA plot.pdf")


## Sample to sample correlation, see https://rockefelleruniversity.github.io/RU_RNAseq/presentations/singlepage/RU_RNAseq_p3.html

sampleDists <- as.dist(1 - cor(log2 (counts(dds,normalized=TRUE)+1), method="pearson"))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)

df <- as.data.frame(colData(dds)[,c("genotype","sample")])

pdf ("Distance between samples plot.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors, annotation_row=df)
dev.off()







