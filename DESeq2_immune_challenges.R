library("tximport")
library("readr")
library("DESeq2")
library("tibble")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("ggrepel")

######################################
# personalize PCA plot function
plotPCA.DESeqTransform.custom = function(object, intgroup="condition",
                                  ntop=500, returnData=FALSE, pcsToUse=1:2)
{
  message(paste0("using ntop=",ntop," top features by variance"))
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, name=colnames(object), colData(object))
  colnames(d)[1:2] <- pcs
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcsToUse]
    return(d)
  }
  
  ggplot(data=d, aes_string(x=pcs[1], y=pcs[2], color="group",shape = "group")) +
    geom_point(size=3) + 
    xlab(paste0(pcs[1],": ",round(percentVar[pcsToUse[1]] * 100),"% variance")) +
    ylab(paste0(pcs[2],": ",round(percentVar[pcsToUse[2]] * 100),"% variance")) +
    coord_fixed()
}


######################################

salmon.dir <- "../Analysis_Draco/salmon"

# control 001-004, exposed 005-008, juvs 009-013
# table with names and specs of samples
samples <- read.table(file.path(".", "samples.txt"), header = TRUE)
rownames(samples) <- paste0("sample", sprintf("%02d",1:13))
samples$date <- as.factor(samples$date)
samples$condition <- as.factor(samples$condition)

# only analyzing adults - immune challenged and controls
samples <- samples[samples$ontogenic_status == "adult",]
samples <- droplevels(samples)

# matrix count files
files <- file.path(salmon.dir,paste("quant_", sprintf("%03d",samples$num),sep=""), "quant.sf")
all(file.exists(files))
names(files) <- rownames(samples)

#### transcript to gene correspondence based on Trinity gene trans map ####
# build the transcript to gene correspondence data.frame
trinity_gt_map <- read.table("./Trinity_output.fasta.gene_trans_map")
names(trinity_gt_map) <- c("gene","transcript")
trinity_gt_map <- trinity_gt_map[c("transcript","gene")]
# call tximport
# countsFromAbundance = "lengthScaledTPM" --> scales counts to library size and transcript length
# trying to correct for
# same sequencing depth but different quantity of RNA sequenced (800ng for adults,400, 200, or 100 ng for juveniles) 
count.matrix <- tximport(files,type = "salmon", countsFromAbundance = "lengthScaledTPM", tx2gene = trinity_gt_map)
names(count.matrix)
head(count.matrix$abundance) #TPM from salmon quant.sf
head(count.matrix$counts) #NumReads from salmon quant.sf
head(count.matrix$length)
head(count.matrix$countsFromAbundance)

### Check that sample names match in both files
all(colnames(count.matrix$counts) %in% rownames(samples))
all(colnames(count.matrix$counts) == rownames(samples))

# build DESeq dataset
# offset matrix not built automatically (see https://github.com/nf-core/rnaseq/issues/499#issuecomment-743384515)
# design to reflect paired samples (ctrl and inf on the same date/immune challenge experiment)
# date = replicate
dds <- DESeqDataSetFromTximport(txi = count.matrix, colData = samples, design = ~ date + condition)

dds$condition <- relevel(dds$condition, ref = "ctrl")

# View(counts(dds))
dds.sf <- estimateSizeFactors(dds)
normalized_counts.dds <- counts(dds.sf, normalized=TRUE)

#### visual inspection of the input data for DESeq2 ####
data.tr <- data.frame(round(count.matrix$counts))

#### quality control on the samples and clustering ####
### Transform counts for data visualization
# blind=TRUE -> unbiased to sample condition information
vst.data <- vst(dds, blind=TRUE)

# PCA
# plot by condition and by date
# then manually keep shape from condition and add colour by date
plotPCA.DESeqTransform.custom(vst.data,intgroup="condition")+
	theme_bw()

ggsave("250513_PCA_shape_condition.pdf",path=".",device = "pdf", dpi=300, width = 300, height = 210,units="mm")

plotPCA(vst.data,intgroup=c("date"))+
  theme_bw()

ggsave("250513_PCA_color_replicate.pdf",path=".",device = "pdf", dpi=300, width = 300, height = 210,units="mm")


# hierarchical clustering
vst_mat <- assay(vst.data) #  from the "SummarizedExperiment" package
### Compute pairwise correlation values
vst_cor <- cor(vst_mat)
head(vst_cor)   ## check the output of cor(), make note of the rownames and colnames
# plot heatmap
pheatmap(vst_cor, annotation = samples)
# same results as PCA: high clustering with the same date/replicate
# and then two by two without anything in common

# run DESeq2 analysis on paired samples
dds <- DESeq(dds)
resultsNames(dds)


# Apply fold change shrinkage (more accurate log2 foldchange estimates)
# https://doi.org/10.1093/bioinformatics/bty895 used apeglm
dds.res.shrunk <- lfcShrink(dds, coef = "condition_inf_vs_ctrl", res=dds.res)


#### extract only significant differentially expressed genes ####

##### p 0.01
summary(dds.res.shrunk, alpha = 0.01)

# out of 499513 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 8, 0.0016%
# LFC < 0 (down)     : 21, 0.0042%
# outliers [1]       : 0, 0%
# low counts [2]     : 464515, 93%
# (mean count < 45)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

padj.cutoff.01 <- 0.01

# subset only significant ones
sig.res.01 <- dplyr::filter(tb.dds.res.shrunk, padj < padj.cutoff.01)
## Order results by padj values
sig.res.01 <- dplyr::arrange(sig.res.01, padj)
sig.res.01 # A tibble: 74 × 6

# print to file significantly DE genes
write.table(sig.res.01,file="Significantly_DE_genes_ctrl_inf_01.csv",sep =",",row.names =  FALSE )

################################################
#### visualization of significantly DEG ####
################################################

# genes
genes.to.plot <- dplyr::pull(sig.res.01,gene)

# retrieve counts for each sample
counts.to.plot.list <- lapply(genes.to.plot, function(gene) {
	plotCounts(dds, gene = gene, intgroup = c("condition","date","ID"), normalized = TRUE, returnData = TRUE)
})

# add names of the genes
for (i in seq_along(genes.to.plot)) {
	counts.to.plot.list[[i]]$gene <- genes.to.plot[i]
}


# Combine into one dataframe
counts.data <- do.call(rbind, counts.to.plot.list)
counts.data$gene <- factor(counts.data$gene)
counts.data$date <- factor(counts.data$date) # should already be factor

# 'gene' factor levels are in the order of padj value
counts.data$gene <- factor(counts.data$gene, levels = unique(counts.data$gene))


# p<0.01
ggplot(counts.data[counts.data$gene %in% sig.res.01$gene,]) +
  geom_point(aes(x = condition, y = count, color = date, shape = condition),size=3,position=position_jitter(w = 0.05,h = 0)) +
  scale_color_manual(values = c("red","orange","#008000","#000045"))+
  scale_y_log10() +
  scale_x_discrete(expand = c(1.5, 0))+
  facet_wrap(vars(gene),
             nrow=1,
             strip.position = "bottom",
             labeller = labeller(gene = function(x) gsub("^TRINITY_", "TRINITY\n", x)))+
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top Significant DE Genes, p < 0.01") +
  theme_minimal() +
  theme(legend.position = "right", 
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 11, angle = 90),
        axis.text.y = element_text(size = 11,angle = 90),
        axis.line.y = element_line(size = .5),
        plot.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(0,"lines"),
        panel.border = element_rect(colour = "grey90", fill=NA, linewidth=0.5),
        
        strip.placement = "outside",
        strip.text = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

ggsave("top_DE_gene_plessthan_01.pdf",path=".",device = "pdf", dpi=300, width = 600, height = 150,units="mm")



#########################
#### volcano plot ####
#######################

## Obtain logical vector where TRUE values denote padj values < 0.01 and fold change > 1 in either direction

res.logical <- as_tibble(rownames_to_column(data.frame(dds.res),var="gene")) %>% 
	dplyr::mutate(threshold_OE = padj < 0.01 & abs(log2FoldChange) >= 0.5)
res.logical$up_down_threshold <- "down"
res.logical[!is.na(res.logical$log2FoldChange) & res.logical$log2FoldChange < 0,]$up_down_threshold <- "up"


## Volcano plot
ggplot(res.logical) +
	geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = interaction(threshold_OE,up_down_threshold))) +
  scale_color_manual(values = c("grey40","red","grey40","darkblue")) +
	xlab("log2 fold change") + 
	ylab("-log10 adjusted p-value") +
  xlim(-2.4,2.4) +
	#scale_y_continuous(limits = c(0,50)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
	theme(legend.position = "none",
				plot.title = element_text(size = rel(1.5), hjust = 0.5),
				axis.title = element_text(size = rel(1.25)))

ggsave("volcanoplot_p0_01.pdf",path=".",device = "pdf", dpi=300, width = 150, height = 75,units="mm")



