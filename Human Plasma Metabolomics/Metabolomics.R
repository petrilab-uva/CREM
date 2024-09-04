library(tidyverse)
library(broom)
library(viridis)
library(vegan)
library(mvnormtest)

#### Read in Genotype data ####

#reading DBC, PROVIDE, and Crypto genotype data by child ID
genotypes<-read.csv("lab_bv_crem_prkca_010719.csv")

#reading in weight/height/age etc. data on children
DBC<-read.csv("DBC_data.csv")
Provide<-read.csv("Provide_data.csv")
Crypto<-read.csv("Crypto_data.csv")

#removing unneeded columns
DBC<- DBC %>% 
  select(-c(VTYPE, DOC, DOB, AGEM, WHOBAZ, NCHSHAZ, NCHSWAZ, NCHSWHZ, WTKG, HTCM))
Provide<- Provide %>% 
  select(c(CHILDID,AgeDay,HAZ,WAZ,WHZ,GENDER))
Crypto<- Crypto %>% 
  select(c(CHILDID,AgeDay,HAZ,WAZ,WHZ,GENDER, SEX))

#Wrangling Gender/Sex in M/F for all data sets
DBC<-DBC %>% 
  rename("SEX" = "GENDER")

DBC<-DBC %>%
  mutate(GENDER = case_when(SEX == 1 ~ "Male", SEX == 2 ~ "Female"))

DBC<-DBC %>% 
  select(-SEX)

Crypto<-Crypto %>% 
  select(-SEX)

#adding cohort before merging into one dataset
DBC<-DBC %>% 
  mutate(Cohort = "DBC")
Crypto<-Crypto %>% 
  mutate(Cohort = "Crypto")
Provide<-Provide %>% 
  mutate(Cohort = "Provide")

#Combining cohorts into one dataset
combodf<-rbind(DBC, Crypto)
metadf<-rbind(combodf, Provide)

meta_geno<-full_join(metadf, genotypes, by = "CHILDID")
meta_geno<-meta_geno %>% 
  select(-c(Genotype_ID,Study))

#read in BEED genotype data
rs483<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/RNAseq_Marie/rs2148483 format for ggsashimi.csv")

rs483$Call<-as.factor(rs483$Call)
rs483<-rs483 %>% 
  filter(Call != "No Call")

rs483<-rs483 %>% 
  mutate(CREM_rs2148483_A = case_when(Call == "Allele 1" ~ 0, Call == "Allele 2" ~ 2, Call == "Heterozygote" ~ 1))

rs483<-rs483 %>% 
  rename("CHILDID" = "Sample")

#Now that I have BEED and DBC/PROVIDE/Crypto data with genotypes. Filter birth cohorts down to one row per child. Then for both dfs eep only relevant columns.

meta_geno<-meta_geno %>% 
  group_by(CHILDID) %>% 
  slice_sample(n=1) 

meta_geno<-meta_geno %>% 
  select(c("CHILDID", "CREM_rs2148483_A"))

meta_geno$CHILDID<-as.character(meta_geno$CHILDID)

rs483<-rs483 %>% 
  select(c("CHILDID", "CREM_rs2148483_A"))

geno<-rbind(meta_geno, rs483)

#### Read in sample data ####

#Data file from Brett Moreau: Kept all features with at least 80% of samples >LOD. Remaining features were imputed (using Random Forest imputation), log2 transformed,  centered, scaled, and normalized to sample volume (different samples had different volumes extracted before analysis).
metab<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/Crem_genotype_waz_haz/Metabolomics/Metabolon untargeted metabolomics imputed_all samples.csv")

#Child AC27176 was identified by metabolon to be an outlier. This ID corresponds to parent sample 6716212. Remove this sample from the data. Also remove "BEED 101 sample"
metab<-metab %>% 
  filter(PARENT.SAMPLE.ID != 6716212) %>% 
  filter(PARENT.SAMPLE.ID != "BEED 101 plasma")

#read in metadata
meta<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/Crem_genotype_waz_haz/Metabolomics/Metabolomics Sample Metadata.csv")

#read in biochemical legend
biochem<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/Crem_genotype_waz_haz/Metabolomics/Biochemical_legend.csv")

#### Metab Result Data Wrangling ####

#Change CHILDID column to match ID numbers in meta_geno file
meta<-meta %>% 
  rename("CHILDID"= "SUBJECT.ID")

meta$CHILDID <- ifelse(meta$CHILDID == "", meta$CLIENT.IDENTIFIER, meta$CHILDID)

#add genotype data to metadata df
meta<-left_join(meta, geno, by = "CHILDID")

#filter children without CREM genotype
meta<-meta %>% 
  drop_na(CREM_rs2148483_A)

n_distinct(meta$CHILDID) #140 unique children have observations 

#keep only necessary columns
meta<-meta %>% 
  select(c("PARENT.SAMPLE.ID", "CHILDID", "CREM_rs2148483_A"))

#merge into metab df
metab<-inner_join(meta, metab, by = "PARENT.SAMPLE.ID")

#only one Bangladeshi BEED child is in df. remove this child

metab<-metab %>% 
  filter(CHILDID != "BC21073")

#save this for later when needed for heat map
metabheat<-metab

n_distinct(metab$CHILDID) #236 samples from 139 children in data set

#filter out Crypto Group A so all metab samples are independent from each other
metab<-metab %>% 
  filter(Group != "GR_CRYPTO_A")
n_distinct(metab$CHILDID) #139 samples from 139 children in data set

#### Exploratory data analysis ####

#Data has already had missing values imputed and been transformed, centered, and scaled.

#PCA plot
metab<-metab %>% 
  filter(CREM_rs2148483_A != "2")

forpca<-metab[,5:723]
pcagroup<-metab[,1:4]

pca<-prcomp(forpca)
plot(pca, type="l")

pcaload<-(pca$rotation)
pcaload<-as.data.frame(pcaload)

pcai<-data.frame(pca$x, Group=pcagroup$Group, Genotype = pcagroup$CREM_rs2148483_A)
pcai$Genotype<-as.factor(pcai$Genotype)
pcai$Group<-as.factor(pcai$Group)

MetabPCAplot <- ggplot(pcai,aes(x=PC1,y=PC2, color=Genotype)) + 
  geom_point(size=4, aes(shape = Group)) + 
  ggtitle("Comparative Metabolomics") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_manual(labels = c("Homozygous Risk", "Heterozygous"), values = c("darkslategray", "darkslategray3"))+
  scale_shape_manual(labels = c("American", "Bangladeshi"), values =c(17, 16)) +
  guides(color=guide_legend("rs2148483 Genotype"), shape=guide_legend("Study Site"))+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22), plot.title=element_text(size=26), legend.text=element_text(size=18), legend.title=element_text(size=18))
#save at jpeg 750x500
plot(MetabPCAplot)

#Check allele frequency by group
metab$CREM_rs2148483_A<-as.numeric(metab$CREM_rs2148483_A)
metab$Group<-as.factor(metab$Group)

metab %>% 
  group_by(Group) %>% 
  summarise(NonRisk=sum(CREM_rs2148483_A), Total=n()*2, MAF=NonRisk/Total)
#17.9% MAF in UVA kids
#24.7% MAF in Bangladeshi kids

metab$CREM_rs2148483_A<-as.factor(metab$CREM_rs2148483_A)

#T.test on homozygous risk vs heterozygous all groups
metabtest<-metab %>% 
  filter(CREM_rs2148483_A != "2") %>% 
  pivot_longer(5:723, names_to= "Metabolite", values_to = "value") %>% 
  group_by(Metabolite) %>% 
  do(tidy(t.test(value ~ CREM_rs2148483_A, data = .)))
metabtestadj<-p.adjust(metabtest$p.value, method = "BH", n =719)

metabtest<-as.data.frame(metabtest)
metabtest<-metabtest %>% 
  mutate(p.adjusted = metabtestadj)

#fold changes all groups
##jump to line 367 here for data wrangling into metaboanalyst format for fc

metabfc<-metab %>% 
  filter(CREM_rs2148483_A != "2") %>% 
  pivot_longer(5:723, names_to= "Metabolite", values_to = "value") %>% 
  group_by(Metabolite, CREM_rs2148483_A) %>%
  summarize(metmeanval = mean(value)) %>% 
  arrange(CREM_rs2148483_A)

View(metabfc)
riskfc<-metabfc$metmeanval[1:719]
hetfc<-metabfc$metmeanval[720:1438]

#data is already Log2 transformed so subtract mean(risk)-mean(nonrisk)
fc<-as_tibble(riskfc-hetfc)
fc<-fc %>% 
  mutate(Metabolite = metabfc$Metabolite[1:719]) %>% 
  rename("Log2FC"= "value")

metabtest$Metabolite <- str_remove(metabtest$Metabolite, "^X+")
fc$Metabolite <- str_remove(fc$Metabolite, "^X+")

metabtest<-left_join(metabtest, fc, by = "Metabolite")

#merge biochem data with stats results

biochem$COMP_ID <- str_remove(biochem$COMP_ID, "^X+")
biochem<-biochem %>% 
  rename("Metabolite" = "COMP_ID")

metabres<-left_join(metabtest, biochem, by = "Metabolite")

write.csv(metabres, "metabolomics_results.csv")                    

#### PerMANOVA Analysis ####
foradonis<-metab %>% 
  filter(CREM_rs2148483_A != "2")

mshapiro.test(as.matrix(foradonis[,5:273]))
#data is not multivariate normal (use non-parametric MANOVA aka PerMANOVA) 

foradonis2<-data.frame((foradonis[,5:273]))
set.seed(2024)
perMAN<-adonis2(foradonis2 ~ CREM_rs2148483_A + Group, data = foradonis, method = "eu")
print(perMAN)

#### Metabo analyst online version ####
#rework metab df to work with Metaboanalyst input

#centered scaled input
key<-as.data.frame(colnames(metab)) %>% 
  rename("Metabolite"= "colnames(metab)")
key$Metabolite <- str_remove(key$Metabolite, "^X+")

key<-left_join(key, biochem, by = "Metabolite")
key$BIOCHEMICAL <- ifelse(is.na(key$BIOCHEMICAL), key$Metabolite, key$BIOCHEMICAL)

#read in another key from metabolon that has info linking biochemical to HMDB and KEGG and CAS
key2<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/Crem_genotype_waz_haz/Metabolomics/Biochem_legend_expanded.csv")

key<-left_join(key, key2, by = "BIOCHEMICAL")

key<-key %>% 
  rename("HMDB" = "Group")

key<-key %>% 
  mutate(Name = ifelse(key$HMDB != "", key$HMDB, ifelse(key$KEGG != "", key$KEGG, ifelse(key$CAS != "", key$CAS, key$PUBCHEM))))

key<-key %>% 
  mutate(Name = ifelse(key$Name == "", key$BIOCHEMICAL, key$Name))
                    
key<-key %>% 
  mutate(Name = ifelse(is.na(key$Name), key$Metabolite, key$Name))

#this HMDB number occurs twice so switch it to CAS to differentiate
key<-key %>% 
  mutate(Name = ifelse(key$Name == "HMDB07219", key$CAS, key$Name))

write.csv(key, "names for conversion.csv")

keyvec<-as.vector(key$Name)
colnames(metab)<-keyvec

##keep only necessary columns for pathway analysis in metaboanalyst

metaboanalyst<-metab %>% 
  select(-c(CHILDID, Group)) %>% 
  filter(CREM_rs2148483_A != 2) %>% 
  mutate(CREM_rs2148483_A = ifelse(CREM_rs2148483_A == 0, "Risk", "Het"))

write.csv(metaboanalyst, "for_metaboanalyst.csv", row.names = FALSE)


#keep only necessary columns for  analysis in metaboanalyst

metaboanalystfc<-forfc %>% 
  select(-c(CHILDID, Group)) %>% 
  filter(CREM_rs2148483_A != 2) %>% 
  mutate(CREM_rs2148483_A = ifelse(CREM_rs2148483_A == 0, "Risk", "Het"))

write.csv(metaboanalystfc, "for_metaboanalyst_fc.csv", row.names = FALSE)

#read in results from statistical analysis [one-factor] on metaboanalyst
#run on "for_metaboanalyst_fc.csv" document using log10 transformation and pareto sccaling
fcres<-read.csv("/Users/audreybrown/Documents/UVA/Petri Lab/Data/Crem_genotype_waz_haz/fold_change_sig_frommetaboanalyst.csv")

fcres<-left_join(fcres, keyfc, by = "Name")


#subset df to only include significantly different metabolites
metabsig<-metabres %>% 
  filter(p.adjusted <= 0.1)
sigmets<-as.vector(metabsig$BIOCHEMICAL)

#filter stratified df to only include genes significantly different in meta analysis
metabresstrat_filt<-metabresstrat %>% 
  subset(BIOCHEMICAL %in% sigmets)

#check calculations for fold change all groups
ggplot(metabresstrat_filt, aes(Group, BIOCHEMICAL, fill= Log2FC)) + 
  geom_tile()+
  scale_fill_gradientn(colours = c("darkslategray","darkslategray4", "darkslategray1", "white", "papayawhip"), values = c(0,.3,.6,.9,1))+
  labs(fill='Fold Change')+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    plot.title = element_text(hjust = 0.5))+
  ggtitle("Metabolite Differential Abundance")+
  xlab("")+
  ylab("")+
  scale_x_discrete(labels = c("American", "Bangladeshi"))

