#------------------ starting here--------


#installing required libraries and edger package 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)
library(limma)
library(RColorBrewer)
#BiocManager::install('mixOmics')
library(mixOmics)
#BiocManager::install("HTSFilter")
library(HTSFilter)
# BiocManager::install("tidyverse")
library(tidyverse)


#setting my directory
setwd("~/Desktop/Disseration HMS x Cambridge ")
directory <- getwd()
#list the files included in the working directory:
dir(directory)


#importing gene count matrix found in Hansen-gene-counts.final.tsv and viewing the table in R
# using the %>% unique() %>% command to remove some of the artifacts created by excel, excel changed some names into dates
fullrawCountTable <- read.csv("Hansen-gene-counts.final.tsv", header=TRUE, sep="\t") %>% unique() %>% 
  filter(gene_name != '1-Mar') %>%  
  filter(gene_name != '2-Mar')
#using this to set the row numbers as the gene names 
rownames(fullrawCountTable) <- fullrawCountTable$gene_name
view(fullrawCountTable)

#removing the first column which contains the gene names since I already set them as row names
fullrawCountTable$gene_name <- NULL
view(fullrawCountTable)

# these are the sample IDs for COMPARISON 3
listOfIDscomparison2 <- c('BCG_1',
                          'BCG_2',
                          'BCG_3',
                          'BCG_4',
                          'BCG_5',
                          'BCG_6',
                          'RhCMV_1',
                          'RhCMV_2',
                          'RhCMV_3',
                          'RhCMV_4',
                          'RhCMV_5',
                          'RhCMV_6')

#sub-setting data needed for this comparison only from the main fullrawCountTable 
counts_rhcmv_comparison <- subset(fullrawCountTable, select = listOfIDscomparison2)
#creating a factor with the samples, BCG versus rhcmv  
# bcg = 1 and rhcmv =2 
DE_BCG_comparison <- factor(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))  
y <- DGEList(counts = fullrawCountTable, group = group)

#viewing data 
fullrawCountTable
head(fullrawCountTable)  
nrow(fullrawCountTable)

#creating sample information table 
sampleInfo <- read.table("design_comparison3.csv", header=TRUE, sep=",", row.names=1)
View(sampleInfo)


#creating a DGElist object
#Using the as.matrix.data.frame function to convert my count table from a data.frame to a matrix 
fullrawCountTable1 <- as.matrix.data.frame(fullrawCountTable, rownames.force = NA)

#viewing data table after matrix conversion 
#setting class type as numeric this is the needed format for later analysis 
view(fullrawCountTable1)
class(fullrawCountTable1) <- "numeric"

#creating dgefull factor to store sample information 
dgeFull <- DGEList(counts_rhcmv_comparison, group=sampleInfo$condition)
dgeFull


#---------------Data exploration and quality assessment
#Extract pseudo-counts
pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts) 

#no need to plot histograms of count data becuase it was done in comparison 1 +2

# MDS for pseudo-counts (using limma package)
dev.new(width=10, height=10, units="in")
plotMDS(pseudoCounts, col=c("green","red"), xlim=c(-6,6))
legend("topright", legend=levels(group), col=3:2, ncol=2, pch=20, cex=0.8)


#heatmap for pseudo-counts (using mixOmics package)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cimColor
dev.new(width=10, height=10, units="in")
cim(sampleDists, color=cimColor, symkey=FALSE)
#this is what we expect, BCG cluster with BCG 

-------------#differential gene expression analysis
#removing genes with zero counts for all samples
  
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                     group=dgeFull$samples$group)
head(dgeFull$counts)

------------#filtering out by expression - part 1
#We filter based on the jaccard index threshold, we first need to claculate the normalization factors  
#estimate the normalization factors using TMM normalizaton method
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
#this result can be considered a table 
head(dgeFull$counts)

#saving the normalized and filtered counts
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)
#boxplot(pseudoNormCounts, col="gray", las=3, ylab="Normalized PseudoCounts" )

#box plot to compare the distribution of counts before and after normalization
boxplot(pseudoCounts, col="gray", las=3, ylab="PseudoCounts", main="Before Normalization", ylim=c(0 , 20))
boxplot(pseudoNormCounts, col="gray", las=3, ylab="Normalized PseudoCounts" , main="After Normalization", ylim=c(0 , 20))


##MDS plot to compare the distribution of counts before and after normalization
plotMDS(pseudoCounts, xlim=c(-7,7), ylim=c(-3, 1.5) , main='Before Normalization',
        
        col=c(rep("green",1),rep("green",1),rep("green",1),rep("green",1), rep("green",1),rep("green",1),
              rep("red",2),rep("red",2),rep("red",2),rep("red",2), rep("red",2),rep("red",2) ))

plotMDS(pseudoNormCounts, xlim=c(-7,7),  ylim=c(-3, 1.5) , main='After Normalization',
        
        col=c(rep("green",1),rep("green",1),rep("green",1),rep("green",1), rep("green",1),rep("green",1),
              rep("red",2),rep("red",2),rep("red",2),rep("red",2), rep("red",2),rep("red",2) ))


------------#Differential expression analysis 
  
#The first major step in the analysis of DGE data using the NB model is to estimate the dispersion parameter
# estimate common and tagwise dispersion


dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull
#plotting the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(dgeFull)
#result: a single estimate for the coefficient of variation is a bad model since 
#tagwise dispersion does not follow the model but instead decreases as the counts per million (CPM) increases.


#perform an exact test for the difference in expression between the two conditions 
dgeTest <- exactTest(dgeFull)
dgeTest
#top 10 differentially expressed genes with their FDR 
topTags(dgeTest, n=10)
dgeFull<- estimateTagwiseDisp(dgeFull)
dgeTest<- exactTest(dgeFull)
de1 <- decideTestsDGE(dgeTest, adjust.method="BH", p.value=0.05)

#Here the entries for -1, 0 and 1 are for down-regulated, non-differentially, and upregulated genes respectivley
summary(de1)


------------#filtering out by expression - part 2
#Independent filtering
# removing low expressed genes
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
dgeTestFilt

#Diagnostic plot for multiple testing
# plot an histogram of unadjusted p-values before vs after filtering 

hist(dgeTest$table[,"PValue"], breaks=50, main="Before Filtering", xlab="PValue" , ylim=c(0,2000))
hist(dgeTestFilt$table[,"PValue"], breaks=50, main="After Filtering", xlab="PValue", ylim=c(0,2000) )

#extract a summary of the differential expression genes
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)


#TOP 5% differentially expressed genes 
# compare the number of differentially expressed genes with and without filtering (risk: 5%)
sum(resNoFilt$table$FDR < 0.05)
# after independent filtering
sum(resFilt$table$FDR < 0.05)



#extract and sort deferentially expressed genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.05,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)
sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)
# write the results in csv files
write.csv(sigDownReg, file="sigDownReg_comparison_3.csv")
write.csv(sigUpReg, file="sigUpReg_comparison_3.csv")


#TOP 1% deferentially expressed genes 
sum(resNoFilt$table$FDR < 0.01)
# after independent filtering
sum(resFilt$table$FDR < 0.01)


#plotting MA plots showing log fold-change versus mean expression
#this MA plot shows deferentially expressed genes and non-differentially expressed genes after filtering 
de1tags12 <- rownames(dgeFull)[as.logical(de1)] 
plotSmear(dgeTest, de.tags=de1tags12, main='Before Filtering')
abline(h = c(-2, 2), col = "blue")


de1tags12 <- rownames(dgeFull)[as.logical(de1)] 
plotSmear(dgeTestFilt, de.tags=de1tags12, main='After Filtering')
abline(h = c(-2, 2), col = "blue")


#creating a Volcano plot for deferentially expressed genes only
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)

# transform the normalized counts in log-counts-per-million
y <- cpm(dgeFull, log=TRUE, prior.count = 1)
head(y)
# select 5% deferentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.05 & 
                                    abs(resFilt$table$logFC)>1.5],]
head(selY)
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
#ERROR
finalHM <- cim(t(selY), color=cimColor, symkey=FALSE) 
table(selY)
dev.new(width=30, height=20, units="in")



# result of the gene clustering
plot(finalHM$ddc, leaflab="none")
abline(h=10, lwd=2, col="pink")


#studying the number of clusters 
# cluster membership for each gene.

geneClust <- cutree(as.hclust(finalHM$ddc), h=10)
head(geneClust) 
length(unique(geneClust))
#naming the genes in cluster 1
names(which(geneClust==1))

de12 <- decideTestsDGE(dgeTestFilt, adjust.method="BH", p.value=0.05)
#Here the entries for -1, 0 and 1 are for down-regulated, non-differentially
#expressed and up-regulated tags respectively.
summary(de12)

