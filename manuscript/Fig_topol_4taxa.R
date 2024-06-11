library(phyloLSnoDist)
library(tictoc)

set.seed(11312)

reps<-10000
n.sites<-c(200,500,2500,10000)
n.scen <- length(n.sites)
correct.tp.ls<-rep(0,n.scen)
correct.tp.nodist<-rep(0,n.scen)
correct.tp.ml<-rep(0,n.scen)


my.tree<-unroot(read.tree(here("analysis", "tree_files", "sim_phylo4-6.tree")))


for(s in 1:length(n.sites)){
  for(i in 1:reps){
    my.align <- simSeq(x = my.tree, l = n.sites[s])

    ols.tree <- phylo.ls(my.align, search.all=TRUE)
    nodist.tree <- phylo.ls.nodist(my.align, search.all=TRUE)
    ml.tree <- phylo.ML(my.align, search.all=TRUE)

    correct.tp.ls[s] <- ifelse(all.equal(my.tree, ols.tree, use.edge.length = F), correct.tp.ls[s]+1, correct.tp.ls[s])
    correct.tp.nodist[s] <- ifelse(all.equal(my.tree, nodist.tree, use.edge.length = F), correct.tp.nodist[s]+1, correct.tp.nodist[s])
    correct.tp.ml[s] <- ifelse(all.equal(my.tree, ml.tree, use.edge.length = F), correct.tp.ml[s]+1, correct.tp.ml[s])

    print(i)

  }
}

save.image(here("analysis", "Fig_topol_4taxa_10k.RData"))


# 42 sec elapsed
# for 2 reps. So 10000 should take 58 hours
# 100 should take about half an hour. 4/6/20 try that now
