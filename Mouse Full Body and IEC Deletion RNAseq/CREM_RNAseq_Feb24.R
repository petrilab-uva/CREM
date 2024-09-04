library(tximport)
library(AnnotationDbi)
library(GenomicFeatures)

#Reading in data from Salmon quant results
#####
#reading in quant files and adding in gene information
dir = "/Users/audreybrown/Documents/UVA/Petri Lab/Data/CREM_Mouse_Feb24"
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
files <- file.path(dir, "salmon", samples$folder, "quant.sf")
names(files) <- paste0(c("162", "169", "170", "176", "177", "178", "179", "198"))

#make TxDB object with gene information
TxDb <- makeTxDbFromGFF(file = "/Users/audreybrown/Documents/UVA/Petri Lab/Data/CREM_Mouse_RNAseq/gencode.vM33.annotation.gtf.gz", format = "gtf")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDb, k, "GENEID", "TXNAME")
head(tx2gene)

#import salmon files with count matrix scaled to TPM
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
names(txi)
head(txi$counts)

txicheck<-as.data.frame(txi$counts)
txicheck["ENSMUSG00000107068.2",]



#DESeq2 and Visualization
#####
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(biomaRt)
library(tidyverse)
library(ggplot2)

#create sample table with genotype information
sampleTable <- data.frame(Genotype = factor(c("Deletion","Wild_Type","Deletion", "Deletion","Wild_Type", "Deletion","Wild_Type", "Wild_Type")))
rownames(sampleTable) <- colnames(txi$counts)
sampleTable$Infection<-as.factor(c("N", "Y", "N", "Y", "Y", "Y", "N", "N"))

#create DESeq object
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~Genotype+Infection)

#prefiltering, keep only genes with at least 1 counts across samples
nrow(dds)
keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]
nrow(dds)
rm(keep)

#prefilter ERa gene and non protein coding genes
biomart<-read.csv("mart_export_mouse_gene_type.csv")
filtergenes<-biomart %>% 
  filter(Gene.type != "protein_coding")
filtergenes2<-biomart[grep("Esr1", biomart$Gene.name), ]
filtergenes<-rbind(filtergenes, filtergenes2)

#Obtain the indices of only desired genes
genesToRemove <- which(!rownames(dds) %in% filtergenes$Gene.stable.ID.version)

#Cut your desired genes in the DESeq object
dds <- dds[genesToRemove, ]
nrow(dds)

#Verify that undesired genes are removed from DESeq object
genesToRemove %in% rownames(dds)

#run DESeq
dds <- DESeq(dds)
res<-results(dds, contrast = c("Genotype", "Deletion", "Wild_Type"))
summary(res)

#Log fold change shrinkage for visualization
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Genotype_Wild_Type_vs_Deletion", type="apeglm")

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="Genotype_Wild_Type_vs_Deletion_filtered_Feb24.csv")

#read back in results file and add gene name info from biomart
resOrder<-read.csv("Genotype_Wild_Type_vs_Deletion_filtered_Feb24.csv")
resOrder<-resOrder %>% 
  rename("Gene.stable.ID.version" = "X")

biomart<-read.csv("mart_export_mouse.csv")
resOrder<-left_join(resOrder, biomart, by = "Gene.stable.ID.version")

#overwrite with df with gene name info
write.csv(as.data.frame(resOrder), 
          file="Genotype_Wild_Type_vs_Deletion_filtered_Feb24.csv")

#visualizations
plotMA(res, ylim=c(-6,6))
plotMA(resLFC, ylim=c(-5,5))

vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

#dispersion plot
plotDispEsts(dds)

#cooks distance boxplot
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#PCA plot
nudge <- position_nudge(y = 1.5)
plotPCA(vsd, intgroup="Infection")+geom_text(aes(label = c("162","169","170","176", "177", "178", "179", "198")), position = nudge)+ggtitle("RNAseq PCA")+theme(plot.title = element_text(hjust = 0.5))+labs(color='Infection')
plotPCA(vsd, intgroup="Genotype")+geom_text(aes(label = c("162","169","170","176", "177", "178", "179", "198")), position = nudge)+ggtitle("RNAseq PCA")+theme(plot.title = element_text(hjust = 0.5))+labs(color='CREM Group')

resOrderdf<-as.data.frame(resOrder)

#save as jpeg 700x700
EnhancedVolcano(resOrderdf,
                lab = resOrderdf$Gene.name,
                selectLab =  c("Il1b", "Cxcl2", "Nlrp3", "Cxcl1"),
                title = 'CREM Deletion vs. Wild-Type',
                legendPosition = "bottom",
                legendLabels = c("Non-Significant","", "Significant by P-value"),
                subtitle ="Differential Expression in Murine Samples",
                caption = "",
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = NA,
                vline = c(-1,1),
                vlineType = "longdash",
                vlineCol = "black",
                vlineWidth = 0.4,
                pCutoff = 0.05,
                xlim = c(-10,10),
                ylim = c(0,12),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                arrowheads = FALSE)

#CREM plot
barplot(assay(dds)["ENSMUSG00000063889.17",],las=2, main=rownames(dds)[ "ENSMUSG00000063889.17"  ], names.arg=c("Deletion","Wild_Type","Deletion", "Deletion","Wild_Type", "Deletion","Wild_Type", "Wild_Type"))

Crem<-plotCounts(dds, gene = "ENSMUSG00000063889.17" , intgroup = "Genotype", returnData = TRUE)

ggplot(Crem, aes(x=Genotype, y=count))+ ylab("Normalized Count") + xlab("Genotype") + ggtitle("Crem Expression")+geom_boxplot()+geom_point()

#IL17a plot
barplot(assay(dds)["ENSMUSG00000025929.5",],las=2, main=rownames(dds)[ "ENSMUSG00000025929.5"  ], names.arg=c("Deletion","Wild_Type","Deletion", "Deletion","Wild_Type", "Deletion","Wild_Type", "Wild_Type"))

il17a<-plotCounts(dds, gene = "ENSMUSG00000025929.5" , intgroup = "Genotype", returnData = TRUE)

ggplot(il17a, aes(x=Genotype, y=count))+ ylab("Normalized Count") + xlab("Genotype") + ggtitle("Crem Expression")+geom_boxplot()+geom_point()


##### 
#GSEA
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(ggbreak)
library(enrichplot)
library(DOSE)
library(stringr)


resgsea<-res[order(-res$stat),]
head(resgsea)

gsea_list<-resgsea$stat
names(gsea_list)<-rownames(resgsea)
names(gsea_list)<-gsub("\\..*","", names(gsea_list))
head(gsea_list)

#running GSEA on genesets with at least 5 genes, adjusted for multiple testing
gsea<-gseGO(gsea_list,
            ont = "BP",
            keyType = "ENSEMBL",
            OrgDb = "org.Mm.eg.db",
            minGSSize = 10,
            pAdjustMethod = "BH",
            eps = 1e-300,
            seed = 2024)

#save as jpeg 700x550
dotplot(gsea, x ="NES", showCategory= 8, split=".sign", title ="GO Term Enrichment: CREM Deletion vs. Wild-Type", font.size=11)+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-3, 2.5)+
  scale_x_break(breaks=c(-2, 1.5), scales =1)

#cnetplot
#order genes by (adjusted) p-value, then filter gsea results to only include these genes
#I filtered because there are way to many core genes per node to be interpretable in a graph
resgsea_pval<-res[order(res$padj),]
gsea_list_pval<-resgsea_pval$padj
names(gsea_list_pval)<-rownames(resgsea_pval)
head(gsea_list_pval)

#pull put core genes
core.genes <- str_split(as.data.frame(gsea)[,"core_enrichment"] , "/")

#filter for list of genes that are sig and change nomeclature
gsea_list_pval2<-as.data.frame(gsea_list_pval)
gsea_list_pval2<-gsea_list_pval2 %>% 
  filter(gsea_list_pval <=0.05)

my.selected.genes <- row.names( gsea_list_pval2)
my.selected.genes <- gsub("\\..*","", (my.selected.genes))

#filter core genes by sig genes
filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),paste, collapse="/")

#write filtered core genes back to gsea object
gsea@result$core_enrichment <- filtered.core.genes

y<-setReadable(gsea, 'org.Mm.eg.db')

#remove decimals from the end of pvalue gene list so it will match with gsea results
names(gsea_list_pval) <- gsub("\\..*","", names(gsea_list_pval))
head(gsea_list_pval)

#save as jpeg 700x600
cnetplot(y, foldChange = gsea_list_pval, showCategory=c("myeloid leukocyte migration"), cex.params = list(category_label =0)) +scale_color_gradient(low = "blue", high = "red") + labs(color="Adjusted P-value", title = "Significantly Differentially Expressed Leading Edge Genes") +theme(plot.title = element_text(hjust = 0.5, size =16)) +guides(size=FALSE)

ydf<-as.data.frame(y$ID)
ydf$Description<-y$Description
gseaplot(y, y$ID[220], title=y$Description[220])

gsea_df<-as.data.frame(gsea@result)
write.csv(gsea_df, "gsea_res_mouse_feb24.csv")
