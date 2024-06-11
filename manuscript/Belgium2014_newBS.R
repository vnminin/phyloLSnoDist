library(phyloLSnoDist)
library(seqinr)
library(here)
library(phangorn)

set.seed(842023)

# see bottom of document for scratch work determining which sequences to pick
# and where the gaps are

# import and drop gapped sites
env_fasta <- read.fasta(here("analysis/BelgiumData/set7_env_Vrancken.fasta"),forceDNAtolower=F)
env_fasta$A96cl7 <- env_fasta$A96cl7[-c(304:306,505:507,1486:1506)]
env_fasta$B90cl22 <- env_fasta$B90cl22[-c(304:306,505:507,1486:1506)]
env_fasta$F02cl15 <- env_fasta$F02cl15[-c(304:306,505:507,1486:1506)]
env_fasta$H96cl02 <- env_fasta$H96cl02[-c(304:306,505:507,1486:1506)]
env_fasta$I99cl29 <- env_fasta$I99cl29[-c(304:306,505:507,1486:1506)]

env_fasta <- env_fasta[which(is.element(names(env_fasta), c("A96cl7", "B90cl22", "F02cl15", "H96cl02", "I99cl29")))]
env_ABFHI <- as.phyDat(ape::as.alignment(env_fasta))

tree_nodist <- phylo.ls.nodist(env_ABFHI, search.all=TRUE)
tree_ols <- phylo.ls(env_ABFHI, search.all=TRUE)
tree_ml <- phylo.ML(env_ABFHI, search.all=TRUE)

# bootstrapping using phangorn
nodist_BS_trees <- bootstrap.phyDat(env_ABFHI, phylo.ls.nodist, bs=100)
ols_BS_trees <- bootstrap.phyDat(env_ABFHI, phylo.ls, bs = 100)
ml_BS_trees <- bootstrap.phyDat(env_ABFHI, phylo.ML, bs=100)

save.image("C:/Users/pchi01/Dropbox/Villanova/Research/LeastSquares/phyloLSnoDist/analysis/Belgium_2014_newBS.RData")



load("C:/Users/pchi01/Dropbox/Villanova/Research/LeastSquares/phyloLSnoDist/analysis/Belgium_2014_newBS.RData")

# Make plot
#tree_nodist$tip.label <- c("A","B","F","H","I")
#tree_ols$tip.label <- c("A","B","F","H","I")
#tree_ml$tip.label <- c("A","B","F","H","I")
# this seems to fuck up the branch labels

# the below results in warnings but it seems to run fine...
pdf(here("analysis", "Fig_env_2014.pdf"), width=9, height=4)
par(mfrow=c(1,3), mar=c(1,2,3,2), oma=c(0,0,0,0))
plotBS(tree_ols, ols_BS_trees, type='u',edge.width=1.5,show.tip.label=F, cex=1.3, p=50)
title("OLS tree", cex.main=2, line=1)
tiplabels("A",tip=1,frame="none",bg=NULL,cex=1.75,adj=c(1,1))
tiplabels("B",tip=2,frame="none",bg=NULL,cex=1.75,adj=c(0.45,0.1))
tiplabels("F",tip=3,frame="none",bg=NULL,cex=1.75,adj=c(0.1,1.05))
tiplabels("H",tip=4,frame="none",bg=NULL,cex=1.75,adj=c(-0.1,0.5))
tiplabels("I",tip=5,frame="none",bg=NULL,cex=1.75,adj=c(1.3,0.1))

plotBS(tree_nodist, nodist_BS_trees, type='u',edge.width=1.5,show.tip.label=F, cex=1.3)
title("New LS tree", cex.main=2, line=1)
tiplabels("A",tip=1,frame="none",bg=NULL,cex=1.75,adj=c(1,1))
tiplabels("B",tip=2,frame="none",bg=NULL,cex=1.75,adj=c(-0.1,0.5))
tiplabels("F",tip=3,frame="none",bg=NULL,cex=1.75,adj=c(0.1,1.05))
tiplabels("H",tip=4,frame="none",bg=NULL,cex=1.75,adj=c(0.5,0.1))
tiplabels("I",tip=5,frame="none",bg=NULL,cex=1.75,adj=c(1.3,0.1))

plotBS(tree_ml, ml_BS_trees, type='u',edge.width=1.5,show.tip.label=F, cex=1.3)
title("ML tree", cex.main=2, line=1)
tiplabels("A",tip=1,frame="none",bg=NULL,cex=1.75,adj=c(1,1))
tiplabels("B",tip=2,frame="none",bg=NULL,cex=1.75,adj=c(-0.1,0.5))
tiplabels("F",tip=3,frame="none",bg=NULL,cex=1.75,adj=c(0.1,1.05))
tiplabels("H",tip=4,frame="none",bg=NULL,cex=1.75,adj=c(0.5,0.1))
tiplabels("I",tip=5,frame="none",bg=NULL,cex=1.75,adj=c(1.3,0.1))

dev.off()



B <- 100
br_table <- matrix(NA, nrow=7, ncol=3)

# First row: OLS
toI <- NULL
for(i in 1:B){
  toI[i] <- ols_bsfixed[[i]]$edge.length[2]
}
toI <- toI[order(toI)]
ols_toI <- paste(signif(tree_ols$edge.length[2], 2), " (", signif(toI[3], 2), ", ", signif(toI[97], 2), ")", sep="")


toH <- NULL
for(i in 1:B){
  toH[i] <- ols_bsfixed[[i]]$edge.length[1]
}
toH <- toH[order(toH)]
ols_toH <- paste(signif(tree_ols$edge.length[1], 2), " (", signif(toH[3], 2), ", ", signif(toH[97], 2), ")", sep="")


int1 <- NULL
for(i in 1:B){
  int1[i] <- ols_bsfixed[[i]]$edge.length[4]
}
int1 <- int1[order(int1)]
ols_int1 <- paste(signif(tree_ols$edge.length[4], 2), " (", signif(int1[3], 2), ", ", signif(int1[97], 2), ")", sep="")


toB <- NULL
for(i in 1:B){
  toB[i] <- ols_bsfixed[[i]]$edge.length[3]
}
toB <- toB[order(toB)]
ols_toB <- paste(signif(tree_ols$edge.length[3], 2), " (", signif(toB[3], 2), ", ", signif(toB[97], 2), ")", sep="")


int2 <- NULL
for(i in 1:B){
  int2[i] <- ols_bsfixed[[i]]$edge.length[7]
}
int2 <- int2[order(int2)]
ols_int2 <- paste(signif(tree_ols$edge.length[7], 2), " (", signif(int2[3], 2), ", ", signif(int2[97], 2), ")", sep="")


toA <- NULL
for(i in 1:B){
  toA[i] <- ols_bsfixed[[i]]$edge.length[5]
}
toA <- toA[order(toA)]
ols_toA <- paste(signif(tree_ols$edge.length[5], 2), " (", signif(toA[3], 2), ", ", signif(toA[97], 2), ")", sep="")


toF <- NULL
for(i in 1:B){
  toF[i] <- ols_bsfixed[[i]]$edge.length[6]
}
toF <- toF[order(toF)]
ols_toF <- paste(signif(tree_ols$edge.length[6], 2), " (", signif(toF[3], 2), ", ", signif(toF[97], 2), ")", sep="")


br_table[,1] <- c(ols_toI, ols_toH, ols_int1, ols_toB, ols_int2, ols_toA, ols_toF)

# Second row: nodist
toI <- NULL
for(i in 1:B){
  toI[i] <- nodist_bsfixed[[i]]$edge.length[2]
}
toI <- toI[order(toI)]
nodist_toI <- paste(signif(tree_nodist$edge.length[2], 2), " (", signif(toI[3], 2), ", ", signif(toI[97], 2), ")", sep="")


toH <- NULL
for(i in 1:B){
  toH[i] <- nodist_bsfixed[[i]]$edge.length[1]
}
toH <- toH[order(toH)]
nodist_toH <- paste(signif(tree_nodist$edge.length[1], 2), " (", signif(toH[3], 2), ", ", signif(toH[97], 2), ")", sep="")


int1 <- NULL
for(i in 1:B){
  int1[i] <- nodist_bsfixed[[i]]$edge.length[4]
}
int1 <- int1[order(int1)]
nodist_int1 <- paste(signif(tree_nodist$edge.length[4], 2), " (", signif(int1[3], 2), ", ", signif(int1[97], 2), ")", sep="")


toB <- NULL
for(i in 1:B){
  toB[i] <- nodist_bsfixed[[i]]$edge.length[3]
}
toB <- toB[order(toB)]
nodist_toB <- paste(signif(tree_nodist$edge.length[3], 2), " (", signif(toB[3], 2), ", ", signif(toB[97], 2), ")", sep="")


int2 <- NULL
for(i in 1:B){
  int2[i] <- nodist_bsfixed[[i]]$edge.length[7]
}
int2 <- int2[order(int2)]
nodist_int2 <- paste(signif(tree_nodist$edge.length[7], 2), " (", signif(int2[3], 2), ", ", signif(int2[97], 2), ")", sep="")


toA <- NULL
for(i in 1:B){
  toA[i] <- nodist_bsfixed[[i]]$edge.length[5]
}
toA <- toA[order(toA)]
nodist_toA <- paste(signif(tree_nodist$edge.length[5], 2), " (", signif(toA[3], 2), ", ", signif(toA[97], 2), ")", sep="")


toF <- NULL
for(i in 1:B){
  toF[i] <- nodist_bsfixed[[i]]$edge.length[6]
}
toF <- toF[order(toF)]
nodist_toF <- paste(signif(tree_nodist$edge.length[6], 2), " (", signif(toF[3], 2), ", ", signif(toF[97], 2), ")", sep="")


br_table[,2] <- c(nodist_toI, nodist_toH, nodist_int1, nodist_toB, nodist_int2, nodist_toA, nodist_toF)


# third row: ML
toI <- NULL
for(i in 1:B){
  toI[i] <- ml_bsfixed[[i]]$edge.length[2]
}
toI <- toI[order(toI)]
ml_toI <- paste(signif(tree_ml$edge.length[2], 2), " (", signif(toI[3], 2), ", ", signif(toI[97], 2), ")", sep="")


toH <- NULL
for(i in 1:B){
  toH[i] <- ml_bsfixed[[i]]$edge.length[1]
}
toH <- toH[order(toH)]
ml_toH <- paste(signif(tree_ml$edge.length[1], 2), " (", signif(toH[3], 2), ", ", signif(toH[97], 2), ")", sep="")


int1 <- NULL
for(i in 1:B){
  int1[i] <- ml_bsfixed[[i]]$edge.length[4]
}
int1 <- int1[order(int1)]
ml_int1 <- paste(signif(tree_ml$edge.length[4], 2), " (", signif(int1[3], 2), ", ", signif(int1[97], 2), ")", sep="")


toB <- NULL
for(i in 1:B){
  toB[i] <- ml_bsfixed[[i]]$edge.length[3]
}
toB <- toB[order(toB)]
ml_toB <- paste(signif(tree_ml$edge.length[3], 2), " (", signif(toB[3], 2), ", ", signif(toB[97], 2), ")", sep="")


int2 <- NULL
for(i in 1:B){
  int2[i] <- ml_bsfixed[[i]]$edge.length[7]
}
int2 <- int2[order(int2)]
ml_int2 <- paste(signif(tree_ml$edge.length[7], 2), " (", signif(int2[3], 2), ", ", signif(int2[97], 2), ")", sep="")


toA <- NULL
for(i in 1:B){
  toA[i] <- ml_bsfixed[[i]]$edge.length[5]
}
toA <- toA[order(toA)]
ml_toA <- paste(signif(tree_ml$edge.length[5], 2), " (", signif(toA[3], 2), ", ", signif(toA[97], 2), ")", sep="")


toF <- NULL
for(i in 1:B){
  toF[i] <- ml_bsfixed[[i]]$edge.length[6]
}
toF <- toF[order(toF)]
ml_toF <- paste(signif(tree_ml$edge.length[6], 2), " (", signif(toF[3], 2), ", ", signif(toF[97], 2), ")", sep="")


br_table[,3] <- c(ml_toI, ml_toH, ml_int1, ml_toB, ml_int2, ml_toA, ml_toF)



library(xtable)
br_conf <- xtable(br_table)
rownames(br_conf) <- c("branch to I", "branch to H", "interior branch 1", "branch to B", "interior branch 2", "branch to A", "branch to F")

br_conf

