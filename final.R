library(WGCNA)
library(tidyverse)
library(DESeq2)
library(GEOquery)
library(gridExtra)
library(ggplot2)
library(CorLevelPlot)


Psoriasis <- getGEO("GSE14905", GSEMatrix = TRUE)[[1]] # read the Psoriasis data set from geo

fPsoriasis <- fData(Psoriasis) # get the fdata of Psoriasis which contains information on genesymbols, description etc....
genesymbols <- fPsoriasis$`Gene Symbol`# assign the gene symbol in the fdata to a variable 54675 elements

exprs_Psoriasis <- data.frame(exprs(Psoriasis)) # get the gene expression data 54675 elements, 82 columns
probeid <- row.names(exprs_Psoriasis) # assign the probeids of the genes to a variable 54675 elements


collapserows <- collapseRows(exprs_Psoriasis,genesymbols,probeid)[[2]]# perform collapse rows function to remove duplicate rows of genes  23,520 entries
exprs1_Psoriasis <- exprs_Psoriasis[collapserows[,2],] # filter the expression data using the collapsed rows output 517440


#OUTLIERS
htree <- hclust(dist(t(exprs1_Psoriasis)), method = "average")
plot(htree)

pca <- prcomp(t(exprs1_Psoriasis))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

outlier <- c('GSM372350')

exprs1_Psoriasis <- exprs1_Psoriasis[,-65]

pdata_Psoriasis  <- pData(Psoriasis)

condition <- pdata_Psoriasis$`characteristics_ch1`

phenodata_Psoriasis <- data.frame(condition)
row.names(phenodata_Psoriasis) <- row.names(pdata_Psoriasis) 
colnames(phenodata_Psoriasis) <- 'condition'

phenodata_Psoriasis <- data.frame(phenodata_Psoriasis[-65,])
row.names(phenodata_Psoriasis) <- colnames(exprs1_Psoriasis)
colnames(phenodata_Psoriasis) <- 'condition'


# Filtering genes 

dds <- DESeqDataSetFromMatrix(countData = round(exprs1_Psoriasis),
                              colData = phenodata_Psoriasis,
                              design = ~ 1) # not spcifying model

dds75 <- dds[rowSums(counts(dds) >= 10) >= 62,]
nrow(dds75) # 13284 genes

# get normalized counts
exp_psoriasis_filt <- assay(dds75)


#SOFT POWER
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(exp_psoriasis_filt,
                         powerVector = power,
                         networkType = "signed")



sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Dermatomyositis<- getGEO("GSE46239", GSEMatrix = TRUE) #load dermatomyositis data

fDermatomyositis <- fData(Dermatomyositis[[1]]) # get the fdata of dermatomyositis which contains information on genesymbols, description etc....
genesymbol_Dermatomyositis <- fDermatomyositis[,"GB_ACC"] # assign the gene symbol in the fdata to a variable

exprs_Dermatomyositis <- data.frame(exprs(Dermatomyositis[[1]])) # get the gene expression data 
probeid_Dermatomyositis <- row.names(exprs_Dermatomyositis) # assign the probeids of the genes to a variable


collapserows_Dermatomyositis <- collapseRows(exprs_Dermatomyositis,genesymbol_Dermatomyositis,probeid_Dermatomyositis)[[2]]# perform collapse rows function to remove duplicate rows of genes
exprs1_Dermatomyositis <- exprs_Dermatomyositis[collapserows_Dermatomyositis[,2],] # filter the expression data using the collapsed rows output
# 23520 genes


#OUTLIERS
htree <- hclust(dist(t(exprs1_Dermatomyositis)), method = "average")
plot(htree)

#PHENO DATA

pheno_Dermatomyositis <- pData(phenoData(Dermatomyositis[[1]]))

pheno_Dermatomyositis <- pheno_Dermatomyositis[,c(2,10)]

condition_Dermatomyositis <- pheno_Dermatomyositis$`characteristics_ch1`
pheno_cond_Dermatomyositis <- data.frame(condition_Dermatomyositis)
row.names(pheno_cond_Dermatomyositis) <- row.names(pheno_Dermatomyositis)


# Filtering the Dermatomyositis genes

dds <- DESeqDataSetFromMatrix(countData = round(exprs1_Dermatomyositis),
                              colData = pheno_cond_Dermatomyositis,
                              design = ~ 1) # not spcifying model

dds75 <- dds[rowSums(counts(dds) >= 10) >= 39,]
nrow(dds75) # 13284 genes

# get normalized counts
exp_dermatomyositis_filt <- assay(dds75)


# SOFT POWER
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(exp_dermatomyositis_filt,
                         powerVector = power,
                         networkType = "signed")



sft.data <- sft$fitIndices
library(ggplot2)
# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


# COMMON PROBES
commonProbesA <- intersect(row.names(exp_psoriasis_filt),row.names(exp_dermatomyositis_filt))



exp_Psoriasis <- exprs1_Psoriasis[which(rownames(exprs1_Psoriasis) %in% commonProbesA),]
exp_Dermatomyositis <- exprs1_Dermatomyositis[which(rownames(exprs1_Dermatomyositis) %in% commonProbesA),]



#adjacency matrix Psoriasis

softpower1 = 12 # Psoriasis
adjacencyPsoriasis = adjacency(t(exp_Psoriasis),power = 12, type= "signed")
dissTOMAPsoriasis = 1-TOMsimilarity(adjacencyPsoriasis, TOMType = "signed")
library(flashClust)
geneTreePsoriasis = flashClust(as.dist(dissTOMAPsoriasis), method = 'average')

plot(geneTreePsoriasis, main = "Gene clustering on TOM-based dissimilarity Psoriasis", labels = FALSE, hang = 0.04)



mcolorh = NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreePsoriasis, pamStage = FALSE, minClusterSize = 30, cutHeight = 0.99, deepSplit = ds, distM = dissTOMAPsoriasis)
  mcolorh = cbind(mcolorh, labels2colors(tree$labels))
}

plotDendroAndColors(geneTreePsoriasis, mcolorh, paste("dsplt =",0:3), main = "",dendroLabels = FALSE)

modulePsoriasis = mcolorh[,1]
plotDendroAndColors(geneTreePsoriasis, modulePsoriasis, "Modules", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendogram and module colors (Psoriasis)")






#adjacency derm

adjacencyderm = adjacency(t(exp_Dermatomyositis),power = 18, type = "signed")
dissTOMAderm = 1-TOMsimilarity(adjacencyderm, TOMType = "signed")
geneTreeDerm = flashClust(as.dist(dissTOMAderm),method = 'average')

plot(geneTreeDerm, main = "Gene clustering on TOM-based dissimilarity Dermatomyositis", labels = FALSE, hang = 0.04)

plotDendroAndColors(geneTreeDerm, modulePsoriasis, "Modules", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendogram and module colors (Dermatomyositis)")


# Module Eigegene

PCsPsoriasis = moduleEigengenes(t(exp_Psoriasis), colors = modulePsoriasis)
ME_psoriasis = PCsPsoriasis$eigengenes
distPCpsoriasis = 1-abs(cor(ME_psoriasis, use = 'p'))
distPCpsoriasis = ifelse(is.na(distPCpsoriasis),0, distPCpsoriasis)
pcTreePsoriasis = hclust(as.dist(distPCpsoriasis),method = 'a')
MDS_Psoriasis = cmdscale(as.dist(distPCpsoriasis),2)
colorsPsoriasis = names(table(modulePsoriasis))

PCsDermatomyositis = moduleEigengenes(t(exp_Dermatomyositis), colors=modulePsoriasis)
ME_Dermatomyositis = PCsDermatomyositis$eigengenes


# Gene Significance


pheno_Psoriasis <- data.frame(ifelse(grepl('normal', phenodata_Psoriasis$condition), 0, 1), row.names = row.names(phenodata_Psoriasis))
colnames(pheno_Psoriasis) <- 'condition'

nsample <- ncol(exp_Psoriasis)
ngenes <- nrow(exp_Psoriasis)

gene.sig.corr <- cor(ME_psoriasis, pheno_Psoriasis, use ='p')
gene.signf.corr.pvals <- corPvalueStudent(gene.sig.corr, nsample)

heatmap.data <- merge(ME_psoriasis, pheno_Psoriasis, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = "condition",
             y = names(heatmap.data)[1:7],
             col = c("blue1", "skyblue", "white", "pink", "red"))

nsample_derm <- ncol(exp_Dermatomyositis)
pheno_Dermatomyositis_binary <- data.frame(ifelse(grepl('DM', pheno_Dermatomyositis$characteristics_ch1), 1, 0), row.names = row.names(pheno_Dermatomyositis))
colnames(pheno_Dermatomyositis_binary) <- 'condition'


gene.sig.corr.derm <- cor(t(exp_Dermatomyositis), pheno_Dermatomyositis_binary, use ='p')
gene.signf.corr.pvals.derm <- corPvalueStudent(gene.sig.corr.derm, nsample_derm)

heatmap.data.derm <- merge(ME_Dermatomyositis, pheno_Dermatomyositis_binary, by = 'row.names')

head(heatmap.data.derm)

heatmap.data.derm <- heatmap.data.derm %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data.derm,
             x = "condition",
             y = names(heatmap.data.derm)[1:7],
             col = c("blue1", "skyblue", "white", "pink", "red"))


# Module Membership

geneModuleMembership1 = signedKME(t(exp_Psoriasis), ME_psoriasis)
colnames(geneModuleMembership1)=paste("PC",colorsPsoriasis,".cor",sep="");
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(exp_Psoriasis)[[2]]);
colnames(MMPvalue1)=c("blue","brown","green","grey","red","turquoise","yellow");
Gene = rownames(exp_Psoriasis)
kMEtable1 = cbind(Gene,Gene,modulePsoriasis)
for (i in 1:length(colorsPsoriasis))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i]) 
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), 
                                                    colnames(MMPvalue1))))


geneModuleMembership2 = signedKME(t(exp_Dermatomyositis), ME_Dermatomyositis) 
colnames(geneModuleMembership1)=paste("PC",colorsPsoriasis,".cor",sep="");
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(exp_Dermatomyositis)[[2]]); 
colnames(MMPvalue2)=c("blue","brown","green","grey","red","turquoise","yellow");
Gene1 <- row.names(exp_Dermatomyositis)
kMEtable2 = cbind(Gene1,Gene1,modulePsoriasis)
for (i in 1:length(colorsPsoriasis))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) 
colnames(kMEtable2)=colnames(kMEtable1)


#ANNOVA SINGLE FACTOR

"The significance of each module to the clinical condition of the datasets 
is calculated to find out which module is more related to the disease conditions using ANOVA test."

multiExprP = list(A1=list(data=t(exp_Psoriasis)),A2=list(data=t(exp_Dermatomyositis)))
multiColorP = list(A1 = modulePsoriasis) 
mp_p=modulePreservation(multiExprP,multiColorP,referenceNetworks=1,verbose=3,networkType="signed", nPermutations = 30)
stats = mp_p$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
write.csv(stats[order(-stats[,2]),c(1:2)],"preservation score all modules.csv")
stats[order(-stats[,2]),c(1:2)]
stats1 <- stats[which(stats$Zsummary.pres>10),]
stats1[,c(1:2)]

grouped =(cbind(ME_psoriasis,pheno_Psoriasis))

grouped %>%
  group_by(condition)

blue_psoriasis <- aov(MEblue ~ condition, grouped)
yellow_psoriasis <- aov(MEyellow ~ condition, grouped)
brown_psoriasis <- aov(MEbrown ~ condition, grouped)
green_psoriasis <- aov(MEgreen ~ condition, grouped)
grey_psoriasis <- aov(MEgrey ~ condition, grouped)
red_psoriasis <- aov(MEred ~ condition, grouped)
turquoise_psoriasis <- aov(MEturquoise ~ condition, grouped)

summary(blue_psoriasis)
summary(brown_psoriasis)
summary(green_psoriasis)
summary(grey_psoriasis)
summary(red_psoriasis)
summary(yellow_psoriasis)
summary(turquoise_psoriasis)

t.test(MEblue ~ condition, grouped)
t.test(MEyellow ~ condition, grouped)
t.test(MEbrown ~ condition, grouped)
t.test(MEgreen ~ condition, grouped)
t.test(MEgrey ~ condition, grouped)
t.test(MEred ~ condition, grouped)
t.test(MEturquoise ~ condition, grouped)


grouped_derm =(cbind(ME_Dermatomyositis,pheno_Dermatomyositis_binary))

grouped_derm %>%
  group_by(condition)

blue_derm <- aov(MEblue ~ condition, grouped_derm)
yellow_derm <- aov(MEyellow ~ condition, grouped_derm)
brown_derm <- aov(MEbrown ~ condition, grouped_derm)
green_derm <- aov(MEgreen ~ condition, grouped_derm)
grey_derm <- aov(MEgrey ~ condition, grouped_derm)
red_derm <- aov(MEred ~ condition, grouped_derm)
turquoise_derm <- aov(MEturquoise ~ condition, grouped_derm)


t.test(MEblue ~ condition, grouped_derm)
t.test(MEbrown ~ condition, grouped_derm)
t.test(MEgreen ~ condition, grouped_derm)
t.test(MEgrey ~ condition, grouped_derm)
t.test(MEred ~ condition, grouped_derm)
t.test(MEturquoise ~ condition, grouped_derm)
t.test(MEyellow ~ condition, grouped_derm)

summary(blue_derm)
summary(brown_derm)
summary(green_derm)
summary(grey_derm)
summary(red_derm)
summary(turquoise_derm)
summary(yellow_derm)

"Based on the Z summary score, ANOVA and t-test results in both Psoriasis and Dermatomyositis
the module Blue is seen to be associated with the disease condition"


"Performing furthur analysis on the genes in module blue to identify the driver genes for the diseases"


# IntraModular connectivity 

"calculated to identify genes which are highly correlated for the disease condition in the module of interest"


module_Blue <-rownames(data.frame(kMEtable1)[which(data.frame(kMEtable1)$Module == 'blue'),])


intramod_connectivity_psoriasis <- t(cor(ME_psoriasis$MEblue, t(exp_Psoriasis), use ='p'))
intramod_connectivity_psoriasis.pval <- data.frame(t(corPvalueStudent(intramod_connectivity_psoriasis, nsample)))
intramod_psoriasis <- cbind(intramod_connectivity_psoriasis, intramod_connectivity_psoriasis.pval)
colnames(intramod_psoriasis) <- c("Intramodular Connectivity","pval_blue")


blue_intramod_connectivity <- intramod_psoriasis[which(rownames(intramod_psoriasis) %in% module_Blue),]


gene_significance_psoriasis <- cor(t(exp_Psoriasis), pheno_Psoriasis, use ='p')
gene_significance_psoriasis.pval <- data.frame(corPvalueStudent(gene_significance_psoriasis, nsample))


module_Blue_derm <-rownames(data.frame(kMEtable2)[which(data.frame(kMEtable2)$Module == 'blue'),])


intramod_connectivity_dermatomyositis <- t(cor(ME_Dermatomyositis$MEblue, t(exp_Dermatomyositis), use ='p'))
intramod_connectivity_dermatomyositis.pval <- data.frame(t(corPvalueStudent(intramod_connectivity_dermatomyositis, nsample_derm)))
intramod_dermatomyositis <- cbind(intramod_connectivity_dermatomyositis, intramod_connectivity_dermatomyositis.pval)
colnames(intramod_dermatomyositis) <- c("Intramodular Connectivity","pval_blue")


blue_intramod_connectivity_dermatomyositis <- intramod_dermatomyositis[which(rownames(intramod_dermatomyositis) %in% module_Blue),]


gene_significance_dermatomyositis <- cor(t(exp_Dermatomyositis), pheno_Dermatomyositis_binary, use ='p')
gene_significance_dermatomyositis.pval <- data.frame(corPvalueStudent(gene_significance_dermatomyositis, nsample_derm))

