## Supplement to Meng et al., 
## Dimension reduction techniques for the integrative analysis of multi-omics data
## ----setup---------------------------------
library(omicade4)
library(genefilter)
library(beeswarm)
library(SciViews)
library(knitr)

## ----Download package from git-----------------------------------------------------
library(devtools)
install_github("aedin/NCI60Example")

## ----loadData------------------------------------------------------------
library(NCI60Example)
data(nci60)

## ----parseData-----------------------------------------------------------
sapply(nci60, dim)
cancerType <- substr(colnames(nci60$mrna), 1, 2)

col <- as.factor(cancerType)
levels(col) <- c("blue", "brown", "gray75")
col <- as.character(col)

## ----fig1_pca------------------------------------------------------------
pc <- pcomp(t(nci60$mrna[subsetGenes, ]), center = TRUE, scale = TRUE)
par(mfrow = c(1,3))
plot(pc, col = "grey40", main = "", ylim = c(0,9))
biplot(pc, font = 2, cex = 0.6, xlim = c(-0.5,0.5), ylim = c(-0.4, 0.4), col = c("grey40", "darkred"))
vectorplot(correlation(pc), labels = sapply(strsplit(subsetGenes, "_"), "[", 1), font = 2, 
           col = "grey40", circle.col = "darkred", main = "", cex = 1, lwd = 2)

## ----fig2_mcia-----------------------------------------------------------
mcoin <- mcia(nci60, cia.nf = 2)
mcoin$mcoa$RV

data.col <- c(rgb(255, 165, 0, alpha = 100, maxColorValue = 255), 
              rgb(0, 255, 0, alpha = 100, maxColorValue = 255), 
              rgb(0, 0, 255, alpha = 100, maxColorValue = 255))

plot.mcia(mcoin, sample.lab = TRUE, sample.color = col, 
          phenovec = cancerType, 
          df.color = data.col, 
          gene.nlab = 0, sample.legend = FALSE)

## ----suppl_Cor_SR_otherCells---------------------------------------------
r <- cor(nci60$prot[, "LE.SR"], nci60$prot[, !colnames(nci60$prot) %in% "LE.SR"])
ct <- sapply(strsplit(colnames(r), "\\."), "[", 1)
bxplot(c(r) ~ ct, ylab="Correlation coefficient")
beeswarm(c(r) ~ ct, add=TRUE, col=c("green", "cyan", "brown"), pch=17, cex=1.2)

