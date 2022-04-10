if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")

library(limma)
library(Biobase)
library(GEOquery)
library(pheatmap)
library(gplots)
library(ggplot2)
library(plyr)
setwd("C:/Users/mirza/Desktop/R_Programming/class_microarray_codes/")
gset=getGEO("GSE27536",GSEMatrix = TRUE,AnnotGPL = TRUE)
gset=getGEO(filename="C:/Users/mirza/Desktop/R_Programming/class_microarray_codes/data/GSE27536_series_matrix.txt.gz",GSEMatrix = TRUE,AnnotGPL = TRUE,destdir = "data/")
class(gset)
dim(gset)
colnames(gset)
sml=c(rep("normalB",12),rep("normalA",12),rep("smnormB",9),rep("smnormA",9),rep("patB",6),rep("patA",6))
length(sml)
sml <- factor(sml)
levels(sml)
ex<-exprs(gset)
#how to know is normal
max(ex)
min(ex)
# if it is not normal run
#ex<-log2(ex+1)
#exprs(gset)<-ex
pdf("Results/boxplots.pdf",width=6)
boxplot(ex)
dev.off()
#if data was not normal
#ex<-normalizeQuantiles(ex)
#exprs(gset)<-ex
x<-normalizeQuantiles(ex)
pdf("Results/boxplot_qnormal.pdf")
boxplot(x)
dev.off()

#Plot correlation HeatMap
pdf("Results/corr_heatmap.pdf")
pheatmap(cor(ex))
dev.off()
#corr_ex<-cor(ex)
#corr_ex[1:4,1:4]
pdf("results/corr_heatmap1.pdf",width = 10,height = 10)
pheatmap(cor(ex),labels_row = sml,labels_col = sml)
dev.off()

# Principal component analysis
pca<-prcomp(ex)
pdf("results/PCA.pdf")
plot(pca)
dev.off()

names(pca)
pca$sdev
colnames(pca$x)
pdf("results/PCAxx.pdf")
plot(pca$x[,1:2])
dev.off()

ex_scale=t(scale(t(ex),scale=F))
mean(ex_scale[1,])
pca2<-prcomp(ex_scale)
pdf("results/PCA_scales.pdf")
plot(pca2)
plot(pca2$x[,1:2])
dev.off()

plot(pca2$rotation)
pc.sample<-data.frame(pca2$r[,1:3],Group=sml)
head(pc.sample)
pdf("results/pca_samples.pdf")
ggplot(pc.sample,aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_bw()
dev.off()

pca2<-prcomp(t(ex_scale))
pc.sample<-data.frame(pca2$x[,1:3],Group=sml)
head(pc.sample)
pdf("results/pca_genes.pdf")
ggplot(pc.sample,aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_bw()
dev.off()

sml <- factor(sml)
levels(sml)
sml
gset$description <- sml
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(sml)
head(design)


#sml
gset$description <- sml
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(sml)
head(design)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(patA-patB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
colnames(tT)
head(ex)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, "Results/ptA_ptB.txt", row.names=F, sep="\t",quote = F)

ms.up=subset(tT,select=c("ID","adj.P.Val","logFC","Gene.symbol"))
ms.up=subset(tT,tT$logFC>1 & tT$adj.P.Val<0.05,select=c("ID","adj.P.Val","logFC","Gene.symbol"))
write.table(ms.up,file= "Results/ptA_ptB.txt", row.names=F, sep="\t",quote = F)
ms.down=subset(tT,tT$logFC< -1 & tT$adj.P.Val<0.05,select=c("ID","adj.P.Val","logFC","Gene.symbol"))
####################################################################################
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
colnames(tT)
ms.up<- subset(tT,logFC> 1 & adj.P.Val< 0.05)
dim(ms.up)
ms.up.genenames<- unique(ms.up$Gene.symbol)
write.table(ms.up.genenames,file = "Results/up_ptA_ptB.txt", quote = F)
write.table(ms.up.genenames,file = "Results/up_ptA_ptB.txt", quote = F,row.names = F,col.names = F)
#ms.up.genenames<- sub("///.*","",ms.up.genenames)
ms.up.genenames <- ms.up.genenames[ms.up.genenames != ""]
ms.up.genenames
ms.up.genenames<- strsplit2(ms.up.genenames,"///")
ms.up.genenames<- unique(ms.up.genenames)
ms.up.genenames<- as.character(ms.up.genenames)
ms.up.genenames
write.table(ms.up.genenames,file = "Results/up_ptA_ptB.txt", quote = F,row.names = F,col.names = F)


ms.down<- subset(tT,logFC< -1 & adj.P.Val< 0.05)
dim(ms.down)
ms.down.genenames<- unique(ms.down$Gene.symbol)
#ms.down.genenames<- sub("///.*","",ms.down.genenames)
ms.down.genenames <- ms.down.genenames[ms.down.genenames != ""]
ms.up.genenames<- strsplit2(ms.up.genenames,"///")
ms.up.genenames<- unique(ms.up.genenames)
ms.up.genenames<- as.character(ms.up.genenames)
ms.up.genenames
write.table(ms.down.genenames,file = "Results/down_ptA_ptB.txt", quote = F,row.names = F,col.names = F)




