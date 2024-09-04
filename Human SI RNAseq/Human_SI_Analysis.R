library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggfortify)
library(ggpubr)

#annotation table
anno<-read.csv("biomart.csv")

#read in children counts
cts <- read.csv("readcount_geneID.csv", row.names=1,check.names = FALSE)

#read in adult counts
cts_BA <- read.csv("readcount adults.csv", row.names=1,check.names = FALSE)

#read in metadata for adults
coldata_BA<-read.csv("coldata adults.csv", row.names = 1)

coldata_BA<-coldata_BA %>% 
  dplyr::select(c(sex, agey, site))
head(cts)

#read in metadata for children
coldata<-read.csv("coldata EED score all.csv", row.names = 1)
coldata<- as.data.frame(coldata)

coldata_clean<-coldata %>% 
  dplyr::select(c(sex, agey, site))

#combine adult with child metadata
coldata_all<-rbind(coldata_clean,coldata_BA)
cts_all<-cbind(cts, cts_BA)
cts_all<-cts_all[,-c(77,94)]

#add in the genotype data for rs2148483
rs483<-read.csv("rs2148483 format for ggsashimi.csv")
rs483$Sample[rs483$Sample == 'BC11007 - 1'] <- 'BC11007'

rs483<-rs483[-124,]
row.names(rs483) <- rs483$Sample
coldata_all<-dplyr::left_join(coldata_all %>%
                                mutate(Sample = rownames(coldata_all)),rs483, by = 'Sample')

row.names(coldata_all) <- coldata_all$Sample
coldata_all <- coldata_all %>% 
  dplyr::select(-Sample)
coldata_all<-coldata_all[-c(77,94),]
coldata_all$Call<-as.factor(coldata_all$Call)

#check format for DDS
rownames(coldata_all) == colnames(cts_all)

#create DEseq dataset 
dds<- DESeqDataSetFromMatrix(countData = cts_all,
                             colData = coldata_all,
                             design = ~sex+site+Call)

#prefiltering, keep only genes with at least 1 counts across samples
nrow(dds)
keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]
nrow(dds)
rm(keep)

#prefilter for protein coding genes
biomart<-read.csv("mart_export_human_gene_type.csv")

filtergenes<-biomart %>% 
  filter(Gene.type != "protein_coding")

#Obtain the indices of only desired genes
genesToRemove <- which(!rownames(dds) %in% filtergenes$Gene.stable.ID)

#Cut your desired genes in the DESeq object
dds <- dds[genesToRemove, ]
nrow(dds)

#Verify that undesired genes are removed from DESeq object
genesToRemove %in% rownames(dds)


#Differential expression analysis

#all individuals
dds <- DESeq(dds)

#risk v het 483
res<-results(dds,contrast = c("Call", "Allele 1","Heterozygote"))
summary(res)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="DEG_CREM_rs483_riskvhet.csv")

#read back in results file and add gene name info from biomart
resOrder<-read.csv("DEG_CREM_rs483_riskvhet.csv")
resOrder<-resOrder %>% 
  rename("Gene.stable.ID" = "X")

resOrder<-left_join(resOrder, biomart, by = "Gene.stable.ID")

#overwrite with df with gene name info
write.csv(as.data.frame(resOrder), 
          file="DEG_CREM_rs483_riskvhet.csv")


library(EnhancedVolcano)
resOrdereddf<-as.data.frame(resOrder)

#save as jpeg 700x700
EnhancedVolcano(resOrdereddf,
                lab = resOrdereddf$Gene.name,
                selectLab = c("CD79B", "CCR7", "CD19", "CXCR4", "CCL19"),
                title = 'rs2148483 Homozygous Risk vs. Heterozygous',
                legendPosition = "bottom",
                legendLabels = c("Non-Significant","", "Significant by P-value"),
                subtitle ="Differential Expression in Human Samples",
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
                xlim = c(-3,3),
                ylim = c(0,8),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                arrowheads =FALSE)

#Gene set enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(ggbreak)
library(enrichplot)
library(DOSE)
library(stringr)

#results for genotype comparisons
resgsea<-res[order(-res$stat),]
head(resgsea)

gsea_list<-resgsea$stat
names(gsea_list)<-rownames(resgsea)
head(gsea_list)

#running GSEA on genesets with at least 10 genes, adjusted for multiple testing
gsea<-gseGO(gsea_list,
            ont = "BP",
            keyType = "ENSEMBL",
            OrgDb = "org.Hs.eg.db",
            minGSSize = 10,
            pAdjustMethod = "BH",
            eps = 1e-300,
            seed = 2023)

#rs483
#save as jpeg 700x550
dotplot(gsea, x ="NES", showCategory= 8, split=".sign", title ="GO Term Enrichment: rs2148483 Homozygous Risk vs. Heterozygous", font.size=14, label_format = 80)+
  theme(plot.title = element_text(hjust = 0.7, size =18),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(hjust=0.62,size=18))+
  xlim(-3, 2.5)+
  scale_x_break(breaks=c(-2, 1.5), scales =1)+
  scale_y_discrete(labels=c("adaptive immune response", "B cell activation", "somatic recombination of immune receptors", "leukocyte cell-cell adhesion", "regulation of leukocyte cell-cell adhesion", "leukocyte mediated immunity", "regulation of T cell activation", "regulation of lymphocyte activation", "precursor metabolite & energy generation", "ribose phosphate biosynthetic process", "purine ribonucleotide biosynthetic process", "aerobic respiration", "mitochondrial gene expression", "ribonucleotide triphosphate biosynthetic process", "mitochondrial translation", "respiratory chain complex"))
                            
write.csv(gsea, "gsea_crem_rs483.csv")

gsea_df<-as.data.frame(gsea)

#cnetplot
#order genes by (unadjusted) p-value, then filter gsea results to only include these genes
#I filtered because there are way to many core genes per node to be interpretable in a graph
resgsea_pval<-res[order(res$padj),]
gsea_list_pval<-resgsea_pval$padj
names(gsea_list_pval)<-rownames(resgsea_pval)
head(gsea_list_pval)
y<-setReadable(gsea, 'org.Hs.eg.db')

my.selected.genes <- names( gsea_list_pval[abs(gsea_list_pval) <= .05] )  
filtered.core.genes <- sapply(lapply(core.genes, function(x) x[x %in% my.selected.genes]),paste, collapse="/")
y@result$core_enrichment <- filtered.core.genes

yfilter<-setReadable(gsea, 'org.Hs.eg.db')

#rs483
#Negative NES filtered only for significant genes
cnetplot(y, foldChange = gsea_list_pval, showCategory=c("adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", "B cell activation", "leukocyte cell-cell adhesion"), cex.params = list(category_label =0))+ scale_color_gradient(low = "blue", high = "red") +labs(color="Adjusted P-value", title = "Significantly Differentially Expressed Leading Edge Genes") +guides(size=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size =18))