### Normalization ####
library(vegan)
Norm_ITS<-read.table("Norm_SIP_cinetique_ITS.txt",row.names = 1, header=TRUE)
df.normalized_ITS<-rrarefy(Norm_ITS,sample=9708)
write.csv2(df.normalized_ITS,file ="df_normalized_SIP_cinetique_ITS.csv")

## Negative control subtraction
tableau.final=array()
t<-read.table("Pour_temoin_neg.txt", h=T)
trans=t(t)
trans=as.data.frame(trans)
nom_col<-as.vector(trans[1,])
nom_col2<-as.data.frame(as.matrix(nom_col),stringsAsFactors=F)
trans2=trans[-1,] 
colnames(trans2)=nom_col2 
dimension=dim(trans2)[2]-1
for ( i in 1:dimension) {
  x<-as.numeric(as.character(trans2$TEXT))
    ifelse(as.numeric(as.character(trans2[,i]))<(10*x),0, as.numeric(as.character(trans2[,i])))->w
  tableau.final<-cbind(tableau.final,w)
}
tableau.final<-cbind(tableau.final,as.numeric(as.character(trans2$TEXT)))
tableau.final[,-1]->tab.tmp
colnames(tab.tmp)=colnames(trans2)
rownames(tab.tmp)=rownames(trans2)
tableau.final2=t(tab.tmp)
tab.shared <- data.frame(Group = colnames(trans2), numOtus = rep(dim(tableau.final2)[2], dim(tableau.final2)[1]), tableau.final2[,1:dim(tableau.final2)[2]])
write.csv2(tab.shared, file = "table_OTU.csv", row.names = F, quote = F)

## Normalization after negative control subtraction
Norm_ITS2<-read.table("Norm_ITS2.txt",row.names = 1, header=TRUE)
View(Norm_ITS2)
df.normalized_ITS<-rrarefy(Norm_ITS2,sample=3248)
View(df.normalized_ITS)
rowSums(df.normalized_ITS)
write.csv2(df.normalized_ITS,file ="df_normalized_SIP_cinetique_ITS2.csv")


### Phyloseq analysis ####
library(phyloseq)
library(ggplot2)
library("dplyr")
library("plyr")
library("ape")
library(vegan)
theme_set(theme_bw())

# Normalization of OTU_table with the mean of fractions
OTU_table<-read.table("OTU_table5.txt", h=T, row.names = 1)
df.normalized_OTU_table<-rrarefy(OTU_table,sample=3252)
rowSums(df.normalized_OTU_table)
write.table(df.normalized_OTU_table,file ="OTU_table6.txt", sep = "\t")

# Create phyloseq object
OTU_table<-read.table("OTU_table6.txt", h=T, row.names = 1)
sample_data<-read.table("env2.txt", h=T, row.names = 1)
sample_data$Fraction <- factor(sample_data$Fraction, levels = c("Light", "Heavy"))
taxa_table<-read.table("tax_table.txt", h=T, row.names = 1)
OTU_table <-as.matrix(OTU_table) 
taxa_table <-as.matrix(taxa_table)
class(OTU_table)
class(taxa_table)
OTU = otu_table(OTU_table, taxa_are_rows = FALSE)
TAX = tax_table(taxa_table)
OTU
TAX
physeq = phyloseq(OTU, TAX)
physeq
sampledata = sample_data(data.frame(sample_data, stringsAsFactors=FALSE))
sampledata
class(sampledata)
random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
plot(random_tree)
physeq1 = merge_phyloseq(physeq, sampledata, random_tree) #Merge new data xith current phyloseq object
physeq1
theme_set(theme_bw())

### Histogramme relative abundance ####
OTU_family<-read.table("relative_abundance_family.txt", h=T)
OTU_family$Family <- factor(OTU_family$Family, levels = c("Cordycipitaceae", "Davidiellaceae",  "Diatrypaceae","Dothioraceae","Erysiphaceae","Leotiaceae", "Phaeosphaeriaceae", "Pleosporaceae", "Saccharomycetaceae", "Corticiaceae", "Cystofilobasidiaceae", "Filobasidiaceae", "Hymenochaetaceae", "Malasseziaceae", "Meruliaceae", "Peniophoraceae", "Polyporaceae", "Trichocomaceae", "Wallemiaceae", "Other_fungi"))

subTi13CF <- subset(OTU_family, Comb3 == "Ti13CF")
subTi13CF$Family <- factor(subTi13CF$Family, levels = c("Cordycipitaceae", "Davidiellaceae",  "Diatrypaceae","Dothioraceae","Erysiphaceae","Leotiaceae", "Phaeosphaeriaceae", "Pleosporaceae", "Saccharomycetaceae", "Corticiaceae", "Cystofilobasidiaceae", "Filobasidiaceae", "Hymenochaetaceae", "Malasseziaceae", "Meruliaceae", "Peniophoraceae", "Polyporaceae", "Trichocomaceae", "Wallemiaceae", "Other_fungi"))
subTi13CF$Comb <- factor(subTi13CF$Comb, levels = c("Ti_Light", "Ti_Heavy", "4h_Light", "4h_Heavy", "10h_Light", "10h_Heavy", "30h_Light", "30h_Heavy"))
p = ggplot(data = subTi13CF, aes(x = Comb, y = Ab, fill = Family)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  xlab("") +
  ylab("Relative abundance (%)") +
  scale_fill_manual(values=c("#FFCCFF","#440154FF", "#CC3333", "#CC66CC", "#66CC33", "#3E4A89FF", "#FFCC33", "#CCFFFF", "#3399FF", "#FF3366", "#1F9E89FF",  "#FF6666", "#00CCCC", "#FF6600", "#3300FF", "#6DCD59FF", "#CC6666", "#B4DE2CFF", "#663366", "#FDE725FF")) +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14), legend.text = element_text(size = 14, face = "italic"), legend.title = element_text(size = 14))
p

subTi13CM <- subset(OTU_family, Comb3 == "Ti13CM")
subTi13CM$Family <- factor(subTi13CM$Family, levels = c("Cordycipitaceae", "Davidiellaceae",  "Diatrypaceae","Dothioraceae","Erysiphaceae","Leotiaceae", "Phaeosphaeriaceae", "Pleosporaceae", "Saccharomycetaceae", "Corticiaceae", "Cystofilobasidiaceae", "Filobasidiaceae", "Hymenochaetaceae", "Malasseziaceae", "Meruliaceae", "Peniophoraceae", "Polyporaceae", "Trichocomaceae", "Wallemiaceae", "Other_fungi"))
subTi13CM$Comb <- factor(subTi13CM$Comb, levels = c("Ti_Light", "Ti_Heavy", "4h_Light", "4h_Heavy", "10h_Light", "10h_Heavy"))
p = ggplot(data = subTi13CM, aes(x = Comb, y = Ab, fill = Family)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  xlab("") +
  ylab("Relative abundance (%)") +
  scale_fill_manual(values=c("#FFCCFF","#440154FF", "#CC3333", "#CC66CC", "#66CC33", "#3E4A89FF", "#FFCC33", "#CCFFFF", "#3399FF", "#FF3366", "#1F9E89FF",  "#FF6666", "#00CCCC", "#FF6600", "#3300FF", "#6DCD59FF", "#CC6666", "#B4DE2CFF", "#663366", "#FDE725FF")) +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14), legend.text = element_text(size = 14, face = "italic"), legend.title = element_text(size = 14))
p

### Enrichment factor ####
library(ggplot2)
theme_set(theme_bw())
EFg<-read.table("EF_Ab_genus.txt", header=TRUE)
EFg$Time <- factor(EFg$Time, levels = c("4h", "10h", "30h"))
EFg$Genus <- factor(EFg$Genus, levels = c("unclassified_Fungi", "Fungi", "unclassified_Basidiomycota", "Malassezia", "Basidiomycota", "unclassified_Pleosporales", "unclassified_Pleosporaceae", "unclassified_Ophiostomataceae", "unclassified_Diatrypaceae", "Sarocladium", "Saccharomyces", "Pezoloma", "Penicillium", "Cyberlindnera", "Cladosporium", "Candida", "Aureobasidium", "Aspergillus", "Alternaria", "Ascomycota"))

p <- ggplot(EFg, aes(EF, Genus))
p + geom_point(aes(colour = Time, shape = Sex, size = Ab)) +
  xlim(-5,25)
str(EFg)
fill1 = c("steelblue", "yellowgreen", "violetred1")
p <- ggplot(EFg, aes(x = EF, y = Genus, size = Ab, colour = Time))
p + geom_point(aes(shape = Sex)) +
  scale_colour_manual(values = fill1) +
  scale_size_area(max_size = 15) +
  labs(size = "Abundance (%)", x = "Enrichment factor", y ="") +
  theme(axis.text.y = element_text(face = c("italic","bold", "italic", "italic", "bold", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "bold"), 
                                   size = c(14,15,14,14,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15)),
        legend.text = element_text(size = 14)) +
  xlim(25,995)

p + geom_point(aes(shape = Sex)) +
  scale_colour_manual(values = fill1) +
  scale_size_area(max_size = 15) +
  labs(size = "Abundance (%)", x = "Enrichment factor", y ="") +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  xlim(25,995)

p <- ggplot(EFg, aes(EF, Genus))
p + geom_point(aes(colour = Time, shape = Sex, size = Ab)) +
  xlim(-5,25)
str(EFg)
fill1 = c("aquamarine2", "aquamarine4", "gray7")
p <- ggplot(EFg, aes(x = EF, y = Genus, size = Ab, colour = Time))
p + geom_point(aes(shape = Sex)) +
  scale_colour_manual(values = fill1) +
  scale_size_area(max_size = 15) +
  labs(size = "Abundance (%)", x = "Enrichment factor", y ="") +
  theme(axis.text.y = element_text(face = c("italic","bold", "italic", "italic", "bold", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "italic", "bold"), 
                                   size = c(14,15,14,14,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15)),
        legend.text = element_text(size = 14)) +
  xlim(-15,25)


### Correlation analyses ####
library(corrplot)

### Correlation Females ####
Corr1 <- read.table("Corr_F4h.txt", header = T, row.names = 1)
M1 <- cor(Corr1)
M1.1 <- M1[20:29,1:19]
p.mat1 <- cor.mtest(Corr1)
p.mat1.1 <- p.mat1$p[20:29,1:19]

Corr2 <- read.table("Corr_F10h.txt", header = T, row.names = 1)
M2 <- cor(Corr2)
M2.2 <- M2[20:29,1:19]
p.mat2 <- cor.mtest(Corr2)
p.mat2.2 <- p.mat2$p[20:29,1:19]

Corr3 <- read.table("Corr_F30h.txt", header = T, row.names = 1)
M3 <- cor(Corr3)
M3.3 <- M3[20:29,1:19]
p.mat3 <- cor.mtest(Corr3)
p.mat3.3 <- p.mat3$p[20:29,1:19]

MG <- rbind(M1.1, M2.2, M3.3)
p.matG <- rbind(p.mat1.1, p.mat2.2, p.mat3.3)
write.table(MG, file = "Resultat_corrplot_female.txt", sep = "\t")
write.table(p.matG, file = "Resultat_corrtest_female.txt", sep = "\t")

corrplot(MG, method="circle", na.label = " ", p.mat = p.matG, tl.cex = 1, tl.col = "black", sig.level = 0.05, insig = "label_sig", pch.cex = 1.6, pch.col = "white") 
#Corrplot with fungi by names
MG2 <- read.table("Resultat_corrplot_female2.txt", h=T, row.names = 1)
MG2 <- as.matrix(MG2)
p.matG2 <- read.table("Resultat_corrtest_female2.txt", h=T, row.names = 1)
p.matG2 <- as.matrix(p.matG2)
corrplot(MG2, method="circle", na.label = " ", p.mat = p.matG2, tl.cex = 1, tl.col = "black", sig.level = 0.05, insig = "label_sig", pch.cex = 1.6, pch.col = "white") 

# Correlation males ####
Corr4 <- read.table("Corr_M4h.txt", header = T, row.names = 1)
M4 <- cor(Corr4)
M4.4 <- M4[18:26,1:17]
p.mat4 <- cor.mtest(Corr4)
p.mat4.4 <- p.mat4$p[18:26,1:17]

Corr5 <- read.table("Corr_M10h.txt", header = T, row.names = 1)
M5 <- cor(Corr5)
M5.5 <- M5[18:25,1:17]
p.mat5 <- cor.mtest(Corr5)
p.mat5.5 <- p.mat5$p[18:25,1:17]

MG3 <- rbind(M4.4, M5.5)
p.matG3 <- rbind(p.mat4.4, p.mat5.5)
write.table(MG3, file = "Resultat_corrplot_male.txt", sep = "\t")
write.table(p.matG3, file = "Resultat_corrtest_male.txt", sep = "\t")

corrplot(MG2, method="circle", na.label = " ", p.mat = p.matG2, tl.cex = 1, tl.col = "black", sig.level = 0.05, insig = "label_sig", pch.cex = 1.6, pch.col = "white") 

#Corrplot with fungi by names
MG4 <- read.table("Resultat_corrplot_male2.txt", h=T, row.names = 1)
MG4 <- as.matrix(MG4)
p.matG4 <- read.table("Resultat_corrtest_male2.txt", h=T, row.names = 1)
p.matG4 <- as.matrix(p.matG4)
corrplot(MG4, method="circle", na.label = " ", p.mat = p.matG4, tl.cex = 1, tl.col = "black", sig.level = 0.05, insig = "label_sig", pch.cex = 1.6, pch.col = "white") 

