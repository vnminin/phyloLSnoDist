# library(tools)
# library(xtable)

# setwd("/homes/peterchi/LeastSquares")
# source("/homes/peterchi/LeastSquares/new_ls.R")
# source("/homes/peterchi/Recombination/dss-func.R")
library(phyloLSnoDist)
library(robustDist)


set.seed(11312)

reps <- 1000
scen <- c(1,3,7,8,12)
n.scen <- length(scen)
n.tips <- 6
n.br <- 2*n.tips - 3

reg.errors <- array(NA,c(reps,n.br,n.scen))
new.errors <- array(NA,c(reps,n.br,n.scen))
counts <- matrix(NA,nrow=n.scen, ncol=reps)

for(l in scen){
   cat('scenario',l,'\n')

   old.err<-as.list(rep(NA,reps))
   new.err<-as.list(rep(NA,reps))

   tree.file<-paste(paste("sim_phylo6-",l,sep=""),".tree",sep="")
   my.tree<-read.tree(here("analysis", "tree_files", tree.file))
   new.brlen<-matrix(NA,nrow=reps,ncol=n.br)
   reg.brlen<-matrix(NA,nrow=reps,ncol=n.br)

   for(i in 1:reps){
     my.align <- as.character(simSeq(my.tree, 1000))

     data.bin<-as.DNAbin(as.alignment(my.align))
     ols.tree <- nnls.tree(as.matrix(dist.dna(data.bin,model="JC69")),unroot(my.tree),trace=0)

     # re-order branches to match input
     for(j in 1:n.br){
       edge<-unroot(my.tree)$edge[j,]
       ind<-which(apply(apply(ols.tree$edge,1,`==`,edge),2,sum)==2)
       reg.brlen[i,j]<-ols.tree$edge.length[ind]
     }


     my.align <- read.phylosim.nuc(my.align)
     optim.out <- new.ls.fit.optimx(rep(0.1, n.br), unroot(my.tree), my.align)
     new.brlen[i,] <- optim.out$par.est
     counts[l,i] <- optim.out$count
   }

   true<-unroot(my.tree)$edge.length
   reg.errors[, , l] <- (reg.brlen - true)/sqrt(true)
   new.errors[, , l] <- (new.brlen - true)/sqrt(true)

}

save.image(here("analysis", "Fig_min3_6taxa.RData"))



