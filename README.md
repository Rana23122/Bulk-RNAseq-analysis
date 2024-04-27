# Bulk-RNAseq-analysis
Bulk RNAseq analysis using Edger package
#loading packages:
library(limma)
library(edgeR)
#loading data
library(readxl)
assesment_data <- read_excel("assesment data.xlsx")
assesment_data.full <- read_excel("assesment data.xlsx")
str(assesment_data)
# Remove the "Gene ID" column
library(dplyr)
assessment_data1 <- assesment_data %>%select(-"Gene ID")
#checking the null and NA reads
is.null(assessment_data1)
is.na(assessment_data1)
sum(is.na(assessment_data1))
#to see the distribution of the data :
hist(as.matrix(assessment_data1))
hist(log2(as.matrix(assessment_data1)))
#help(assessment_data1)
assessment_data1Groups <- c(
  "thy"   = c("thy_1", "thy_2", "thy_3"),
  "ckit"  = c("ckit_1", "ckit_2", "ckit_3"),
  "alpha" = c("alpha_1", "alpha_2", "alpha_3"))

# Create a rep grouping factor
rep_groups <- c(rep("thy", 3), 
                rep("ckit", 3), 
                rep("alpha", 3))

# Create a new DGEList object with the new grouping factor
d <- DGEList(counts = assessment_data1, 
             group = factor(rep_groups))
#filtering
dim(d)
d.full <- d # keep the old one in case we mess up
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
#Normalizing
d <- calcNormFactors(d, method="TMM")
#Data exploration
plotMDS(d, col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=23)
#GLM estimates of dispersion => For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25]. The CR method is derived to overcome
#the limitations of the qCML method as mentioned above. It takes care of multiple factors by
#fitting generalized linear models (GLM) with a design matrix
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d1 <- estimateGLMCommonDisp(d,design.mat)
d1<- estimateGLMTrendedDisp(d1,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d1 <- estimateGLMTagwiseDisp(d1,design.mat)
plotBCV(d1)
#GLM testing for differential expression
design.mat
fit <- glmQLFit(d1, design.mat)
# compare (group 1 - group 2) to 0:
lrt12 <- glmQLFTest(fit, contrast=c(1,-1,0))
lrt13 <- glmQLFTest(fit, contrast=c(1,0,-1))
lrt23 <- glmQLFTest(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)
topTags(lrt13, n=10)
topTags(lrt23, n=10)
#grouplrt12
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05 )
de2tags12 <- rownames(d1)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")
summary(de2)
View(de2)
# Save as CSV
write.csv(de2, file = "de2.csv", row.names = TRUE)
#grouplrt13
de3 <- decideTestsDGE(lrt13, adjust.method="BH", p.value = 0.05 )
de3tags13 <- rownames(d1)[as.logical(de3)]
plotSmear(lrt13, de.tags=de3tags13)
abline(h = c(-2, 2), col = "blue")
summary(de3)
View(de3)
# Save as CSV
write.csv(de3, file = "de3.csv", row.names = TRUE)
#grouplrt23
de4 <- decideTestsDGE(lrt23, adjust.method="BH", p.value = 0.05 )
de4tags23 <- rownames(d1)[as.logical(de4)]
plotSmear(lrt23, de.tags=de4tags23)
abline(h = c(-2, 2), col = "blue")
summary(de4)
View(de4)
# Save as CSV
write.csv(de4, file = "de4.csv", row.names = TRUE)
