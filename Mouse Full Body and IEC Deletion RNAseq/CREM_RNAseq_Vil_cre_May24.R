library(tximport)
library(AnnotationDbi)
library(GenomicFeatures)

#Reading in data from Salmon quant results
#####
#reading in quant files and adding in gene information
dir = "/Users/audreybrown/Documents/UVA/Petri Lab/Data/CREM_Mouse_Feb24"
list.files(dir)
samples <- read.table(file.path(dir, "samples_May24_Vil.txt"), header = TRUE)
samples
files <- file.path(dir, "salmon_May24_Vil", samples$folder, "quant.sf")
names(files) <- paste0(c("247", "248", "336", "337", "392", "393", "473", "474", "487", "488"))

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


sampleTable <- data.frame(Genotype = factor(c("Wild_Type", "Deletion", "Wild_Type", "Deletion", "Wild_Type", "Deletion", "Deletion", "Wild_Type","Wild_Type", "Deletion")))
rownames(sampleTable) <- colnames(txi$counts)
sampleTable$Infection<-as.factor(c("Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N"))

#create DESeq object
#infections status not included because there is only one positive sample per genotype
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~Genotype)

#prefiltering, keep only genes with at least 1 counts across samples
nrow(dds)
keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]
nrow(dds)
rm(keep)

#prefilter ERa gene and lncRNA and pseudogenes
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
summary(resLFC)

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="Genotype_Wild_Type_vs_Deletion_Vil_cre_filtered_May24.csv")

#read back in results file and add gene name info from biomart
resOrder<-read.csv("Genotype_Wild_Type_vs_Deletion_Vil_cre_filtered_May24.csv")
resOrder<-resOrder %>% 
  rename("Gene.stable.ID.version" = "X")

biomart<-read.csv("mart_export_mouse.csv")
resOrder<-left_join(resOrder, biomart, by = "Gene.stable.ID.version")

#overwrite with df with gene name info
write.csv(as.data.frame(resOrder), 
          file="Genotype_Wild_Type_vs_Deletion_Vil_cre_filtered_May24.csv")

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
plotPCA(vsd, intgroup="Infection")+geom_text(aes(label = c("247", "248", "336", "337", "392", "393", "473", "474", "487", "488")))+ggtitle("RNAseq PCA")+theme(plot.title = element_text(hjust = 0.5))+labs(color='Infection')
plotPCA(vsd, intgroup="Genotype")+geom_text(aes(label = c("247", "248", "336", "337", "392", "393", "473", "474", "487", "488")))+ggtitle("RNAseq PCA")+theme(plot.title = element_text(hjust = 0.5))+labs(color='CREM Group')

resOrderdf<-as.data.frame(resOrder)

#save as jpeg 700x700
EnhancedVolcano(resOrderdf,
                lab = resOrderdf$Gene.name,
                selectLab =  c("Erc2", "B3galt1", "Cx3cl1"),
                title = 'IEC CREM Deletion vs. Wild-Type',
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
                labSize = 6.0,
                xlim = c(-8,8),
                ylim = c(0,3.2),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                arrowheads = FALSE)

#CREM plot
barplot(assay(dds)["ENSMUSG00000063889.17",],las=2, main=rownames(dds)[ "ENSMUSG00000063889.17"  ], names.arg=c("Wild_Type", "Deletion", "Wild_Type", "Deletion", "Wild_Type", "Deletion", "Deletion", "Wild_Type","Wild_Type", "Deletion"))

Crem<-plotCounts(dds, gene = "ENSMUSG00000063889.17" , intgroup = "Genotype", returnData = TRUE)

ggplot(Crem, aes(x=Genotype, y=count))+ ylab("Normalized Count") + xlab("Genotype") + ggtitle("Crem Expression Vil-CRE")+geom_boxplot()+geom_point()

Crem<-plotCounts(dds, gene = "ENSMUSG00000063889.17" , intgroup = "Infection", returnData = TRUE)
ggplot(Crem, aes(x=Infection, y=count))+ ylab("Normalized Count") + xlab("Infection") + ggtitle("Crem Expression Vil-CRE")+geom_boxplot()+geom_point()

#IL17a plot
barplot(assay(dds)["ENSMUSG00000025929.5",],las=2, main=rownames(dds)[ "ENSMUSG00000025929.5"  ], names.arg=c("Wild_Type", "Deletion", "Wild_Type", "Deletion", "Wild_Type", "Deletion", "Deletion", "Wild_Type","Wild_Type", "Deletion"))

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
dotplot(gsea, x ="NES", showCategory= 8, split=".sign", title ="GO Term Enrichment: CREM IEC Deletion vs. Wild-Type", font.size=11)+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-2.5, 3.5)+
  scale_x_break(breaks=c(-1.5, 2.5), scales =1)

gsea_df<-as.data.frame(gsea@result)
write.csv(gsea_df, "gsea_res_mouse_Vil_feb24.csv")
