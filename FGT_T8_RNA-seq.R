#Load the required libraries & load the files for tutorial8
library(Rsubread)
library(edgeR)
library(limma)
library(gplots)
library(DESeq2)
library(affy)
library(QuasR)
load("FGT_T8_RNA-seq.Rdata")

#The input files were generated using this code- 
#out and provided for information only
#alignment_files <- list.files(path = "./aligned_merged",pattern="*.bam",full.names=T)

#Generate a simple QC file to confirm the QC results
#cl <- makeCluster(15)
#qQCReport(alignment_files, 
#	pdfFilename="results/qc_report.pdf",clObj=cl,useSampleNames=T)
#stopCluster(cl)

#use featureCounts to load the BAM files
# this object has slots "$" to access
# mytable_features <- featureCounts(files=alignment_files,
# 	annot.ext="./omes/annot/Mus_musc ulus.GRCm38.99.gtf",
# 	isGTFAnnotationFile = TRUE, 
# 	countMultiMappingReads=F, 
# 	minMQS=30, 
# 	nthreads=10,
# 	minOverlap=10)

# Access count data and column names
table_rnaseq <- (mytable_feaures)$counts
colnames_rnaseq <- colnames(table_rnaseq)

# Shorten sample names eg"X.home.guest101.new_rnaseq.aligned.SRR5054353.sorted.bam" 
adf <- read.table("mychromosome_samples4.txt",
	sep='\t', 
	row.names=1,
	fill=T,
	header=T) 
# Order the elements in the table 
idx <- match(colnames_rnaseq, rownames(adf)) 
colnames(table_rnaseq) <- adf$ShortName[idx] 
adf <- adf[idx,] 
rownames(adf) <- colnames(table_rnaseq) 


# remove XO samples
keep_cols <- adf$Sex[idx] != "O"
# Create the required object for DESeq analysis
# The design for differential analysis is fed into this object
# The adf annotates each sample 
dds_rnaseq <- DESeqDataSetFromMatrix(
	countData=table_rnaseq[,keep_cols],
	colData=adf[keep_cols,],
	design=~ Line + Sex)

# check the object dimensions 
dim(dds_rnaseq) 
head(rownames(dds_rnaseq))
colnames(dds_rnaseq) 
# check the object type
# DESeqDataSet class, type S4
class(dds_rnaseq) 
typeof(dds_rnaseq)
# check what this object contains 
slotNames(dds_rnaseq) 
colnames(colData(dds_rnaseq)) 
colnames(assay(dds_rnaseq)) 

# Filter remove genes below threshold
# first filter by counts 
keep <- rowSums(counts(dds_rnaseq)) > 1 
dds_rnaseq <- dds_rnaseq[keep,]
nrow(dds_rnaseq) 
# second filter - require >10 in >= 3 samples 
keep <- rowSums(counts(dds_rnaseq) >= 10) >= 3 
dds_rnaseq <- dds_rnaseq[keep,] 
nrow(dds_rnaseq)

# size factor estimation for normalisation
# these are used to scale by library size when data is exported. 
dds_rnaseq <- estimateSizeFactors(dds_rnaseq) 
sizeFactors(dds_rnaseq)

# DESeq2 methods for data normalization
# Normalized log coiunts
counts_rnaseq <- log2(counts(dds_rnaseq, normalized=TRUE))
# fpm is similar to cpm and discussed here 
# https://rdrr.io/bioc/DESeq2/man/fpm.html
# Here the data is log transformed
# Also normalize by the sizeFactors
fpm_rnaseq <- log(fpm(dds_rnaseq))
# rlog and vsd are discussed here 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
vsd_rnaseq <- vst(dds_rnaseq, blind=T)
head(vsd_rnaseq, 3)
rld_rnaseq <- rlog(dds_rnaseq, blind=T)
head(rld_rnaseq, 3)

# rlog is just used for visualization
# but DE uses counts

# Show normalised data using boxplots
boxplot(counts_rnaseq,
        main="Counts")
boxplot(fpm_rnaseq,
        main="Fragments per million")
boxplot(assay(vsd_rnaseq),
        main="Variance-stabilizing transformed")
boxplot(assay(rld_rnaseq),
        main="Regularized-logarithm transformed")
# but suppose we look at non-normalised
na.rm=T
counts_un <- counts(dds_rnaseq, normalized=FALSE)
boxplot(log2(counts_un),
	main="Counts, non-normalized")
counts_un_rlog <- rlog(counts_un)

# Make mva.pairs-this is slow “mva.pairs.png” is provided 
#png(filename="mva.pairs.png",width=2000, height=2000)
#mva.pairs(counts_un_rlog)
#dev.off()

# Make a PCA plot
library(scatterplot3d)
library(ggplot2)
pca <- prcomp(t(na.omit(assay(rld_rnaseq))),
	scale=T)
s3d <- scatterplot3d(pca$x[,1:3], 
	pch=19, 
	color=colData(dds_rnaseq)$Colour)
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, 
     labels=colnames(rld_rnaseq),
     pos=3, offset=0.5, cex=0.5)
# or a quick plot in 2d
qplot(pca$x[,1],pca$x[,2],
	xlab="PCA1",
    ylab="PCA2",
    color=colData(dds_rnaseq)$Colour)

# Output a distance matrix
sampleDists <- dist(t(assay(rld_rnaseq)))
sampleDists
# Generate heatmap of distance matrix
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
    main="Heatmap of samples distance matrix")

# Alternative PCA plot
plotPCA(rld_rnaseq, 
    intgroup = c("Sex", "Line"))

	#----------------DEG using DESeq2

# Stats similar to Limma but for count data
# Two-factor analysis (Sex, Line)
dds_rnaseq <- DESeq(dds_rnaseq)
# Specify contrasts
result_sex <- results(dds_rnaseq,
	contrast=c("Sex","M","F"))
result_cell_line <- results(dds_rnaseq,
    contrast=c("Line","BC","CB"))
# Plot the DEG, plotMA()
plotMA(result_sex, 
	main="DESeq2", 
	ylim=c(-2,2))
table(result_sex$padj < 0.01)
table(result_cell_line$padj < 0.01)
table(result_sex$padj < 0.05)
table(result_cell_line$padj < 0.05)
result_sex_selected <- subset(result_sex, padj < 0.05)
result_sex_selected <- result_sex_selected[order(
	abs(result_sex_selected$log2FoldChange), 
	decreasing=TRUE), ]
head(result_sex_selected)

# Make heatmaps from each output table
# Make a heatmap for the top 50
top50 <- rownames(result_sex_selected)[1:50]
pheatmap(assay(rld_rnaseq)[top50,],
    scale="row",
    show_rownames=T,
    main="Heatmap of Top 50 DEGs: M versus F")
# Repeat for all genes
pheatmap(assay(rld_rnaseq)[rownames(result_sex_selected),],
    scale="row",
    show_rownames=F,
    main="Heatmap of all DEGs: M versus F")

# Repeat for the cell lines
result_cell_selected <- subset(result_cell_line, 
	padj < 0.05)
result_cell_selected <- result_cell_selected[order(
	abs(result_cell_selected$log2FoldChange), 
	decreasing=TRUE), ]
head(result_cell_selected)
top50c <- rownames(result_cell_selected)[1:50]
# Heatmap for top 50 by cell lines
pheatmap(assay(rld_rnaseq)[top50c,],
    scale="row",
    show_rownames=T,
    main="Heatmap of Top 50 DEGs: Line BC vs CB")
# Repeat for all by cell lines
pheatmap(assay(rld_rnaseq)[rownames(result_cell_selected),],
    scale="row",
    show_rownames=F,
    main="Heatmap of all DEGs: Line BC vs CB")

	#---------------GENE ANNOTATION

# Download Ensembl annotation using BiomaRt and rename the samples
library(biomaRt)
# UK ensembl is being updated so we use a USA mirror, "useast.ensembl.org"
ensembl_host <- "uswest.ensembl.org"
head(biomaRt::listMarts(host = ensembl_host), 15)
head(
	biomaRt::listAttributes(
	biomaRt::useDataset(
		dataset = "mmusculus_gene_ensembl",
		mart = useMart("ENSEMBL_MART_ENSEMBL",
			host = ensembl_host)
		)), 40) 
mart <- biomaRt::useDataset(
	dataset = "mmusculus_gene_ensembl",
	mart = useMart("ENSEMBL_MART_ENSEMBL",
		host = ensembl_host))
# resultAnnot <- biomaRt::getBM(
# 	values=rownames(dds_rnaseq),
# 	attributes = c("ensembl_gene_id","external_gene_name",
# 		"chromosome_name","start_position","end_position",
# 		"description","strand"),
# 	filters="ensembl_gene_id",
# 	mart=mart)

# Merge annotation with input data
names <- resultAnnot[,1]
resultAnnot <- as.data.frame(resultAnnot)
rownames(resultAnnot) <- names
idx <- match(rownames(dds_rnaseq), rownames(resultAnnot))
# Make sure annotation is in same order
all(rownames(dds_rnaseq) == rownames(resultAnnot))
grr <- resultAnnot[match(
	rownames(dds_rnaseq), 
	resultAnnot$ensembl_gene_id),]
all(rownames(dds_rnaseq) == rownames(grr))
resultAnnot <- grr
all(rownames(dds_rnaseq) == rownames(resultAnnot))
# make the nice names
nice_names <- paste(resultAnnot$ensembl_gene_id,
	resultAnnot$external_gene_name,
	sep = '_')
resultAnnot$nice_names <- nice_names
head(resultAnnot)
all(rownames(dds_rnaseq) == rownames(resultAnnot))
# check names
rld_rnaseq <- rlog(dds_rnaseq, blind = TRUE)
idx2 <- match(
	rownames(result_sex_selected)[1:50],
	rownames(dds_rnaseq))
plotme <- (rld_rnaseq)[rownames(result_sex_selected)[1:50],]
rownames(plotme) <- resultAnnot$nice_names[idx2]

# Make heatmap with candidate genes
png(filename="heatmap_candidates.png",
	width=600, height=1000)
pheatmap(assay(plotme),
    scale="row",
    fontsize_row=10,
    cellheight=12, cellwidth=12,
    treeheight_row=40, treeheight_col=40,
    main="Heatmap of candidate genes")
dev.off()

# But what about the marker profiles from known markers?
# Load some developmental markers to assess the expression profile
mkr <- read.table("markers.txt",
     stringsAsFactors = F,
     fill=T, row.names = 1,
     header=T, sep="\t")
rownames(mkr)
#%in% This makes sure the marker gene names are rownames in the object
valid_names <- rownames(mkr)[rownames(mkr) %in% rownames(dds_rnaseq)]
plotmetoo <- rld_rnaseq[valid_names,]
idx3 <- match(rownames(plotmetoo),rownames(dds_rnaseq))
rownames(plotmetoo) <- resultAnnot$nice_names[idx3]
png(filename="heatmap_markers.png",width=800, height=1000)
pheatmap(assay(plotmetoo),
     scale="row",
     fontsize_row=10,
     cellheight=12, cellwidth=12,
     treeheight_row=40, treeheight_col=40,
     main="Heatmap of developmental marker genes")
dev.off()

# Finally add column annotation
mkr_bit <- mkr[valid_names,]
rownames(mkr_bit) <- resultAnnot$nice_names[idx3]
Sys.setlocale("LC_ALL", "C")
png(filename="heatmap_markers2.png", width=800, height=1000) 
pheatmap(assay(plotmetoo),
         scale="row",
         fontsize_row=10,
annotation_row=as.data.frame(mkr_bit[,c("LongState","ShortState")]),
         cellheight=12, cellwidth=12,
         treeheight_row=40, treeheight_col=40,
         main="Heatmap of developmental marker genes (with col annotation)")
dev.off()

save.image("FGT_T8_RNA-seq.Rdata")

