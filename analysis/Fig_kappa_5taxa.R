library(phyloLSnoDist)
library(here)

set.seed(4142020)
reps <- 100
kappa <- c(0.25, 0.5, 2, 4)
scen <- c(1,3,9) # tree scenarios
n.br <- 7


new.brlen <- array(NA,c(reps, n.br, length(kappa), length(scen)))
reg.brlen <- array(NA,c(reps, n.br, length(kappa), length(scen)))
reg.errors <- array(NA,c(reps,n.br, length(kappa), length(scen)))
new.errors <- array(NA,c(reps,n.br, length(kappa), length(scen)))
kappa.est <- array(NA,c(reps, length(kappa), length(scen)))

for(l in 1:length(scen)){
  for(k in 1:length(kappa)){
    for(i in 1:reps){
      tree.file<-paste(paste("sim_phylo5-",scen[l],sep=""),".tree",sep="")
      my.tree<-read.tree(here("analysis", "tree_files", tree.file))
      my.align <- simSeq(x = my.tree, l = 1000, Q = c(1, kappa[k], 1, 1, kappa[k], 1))

      D <- as.matrix(dist.dna(as.DNAbin(my.align),model="K80"))
      ols.tree <- nnls.tree(D,unroot(my.tree),trace=0)

      # re-order branches to match input
      for(j in 1:n.br){
        edge<-unroot(my.tree)$edge[j,]
        ind<-which(apply(apply(ols.tree$edge,1,`==`,edge),2,sum)==2)
        reg.brlen[i,j,k,l]<-ols.tree$edge.length[ind]
      }

      new.K80.out <- new.ls.fit.K80(unroot(my.tree), my.align) # 4/8/20 maybe this helps, but check this...
      new.brlen[i, ,k,l] <- new.K80.out$par.est[1:n.br]
      kappa.est[i,k,l] <- new.K80.out$par.est[n.br+1] # 4/8/20 got some crazy shit here. Even just a few times is problematic...
    }
  }
}

save.image(here("analysis", "Fig_kappa_5taxa.RData"))
