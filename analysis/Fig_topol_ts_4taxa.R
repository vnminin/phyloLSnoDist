library(phyloLSnoDist)
library(tictoc)

set.seed(11312)

reps<-100
n.sites<-c(200,500,2500,10000)
n.scen <- length(n.sites)
correct.tp.ls<-rep(0,n.scen)
correct.tp.nodist<-rep(0,n.scen)


my.tree1<-unroot(read.tree(here("analysis", "tree_files", "sim_phylo4-1.tree")))
my.tree2<-my.tree1
my.tree2$tip.label <- my.tree2$tip.label[c(1,4,3,2)]


for(s in 1:length(n.sites)){
  for(i in 1:reps){
    my.align1 <- simSeq(x = my.tree1, l = (n.sites[s]/2), Q=c(0,1,0,0,1,0))
    my.align2 <- simSeq(x = my.tree2, l = (n.sites[s]/2), Q=c(1,0,1,1,0,1))

    my.align <- phyDat(cbind(as.character(my.align1), as.character(my.align2))[,sample(n.sites[s])])

    ols.tree <- phylo.ls(my.align, search.all=TRUE, model="TS") # add ts stuff here (4/20/20 not done at all yet)
    nodist.tree <- phylo.ls.nodist(my.align, search.all=TRUE, ts=TRUE)

    correct.tp.ls[s] <- ifelse(all.equal(my.tree1, ols.tree, use.edge.length = F), correct.tp.ls[s]+1, correct.tp.ls[s])
    correct.tp.nodist[s] <- ifelse(all.equal(my.tree1, nodist.tree, use.edge.length = F), correct.tp.nodist[s]+1, correct.tp.nodist[s])

    print(i)

  }
}

save.image(here("analysis", "Fig_topol_ts_4taxa.RData"))


# 42 sec elapsed
# for 2 reps. So 10000 should take 58 hours
# 100 should take about half an hour. 4/6/20 try that now
