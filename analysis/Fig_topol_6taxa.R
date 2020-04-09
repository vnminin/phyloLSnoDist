library(phyloLSnoDist)
library(tictoc)

set.seed(11312)

reps<-100
n.sites<-c(200,1000,10000,100000)
n.scen <- length(n.sites)
correct.tp.ls<-rep(0,n.scen)
correct.tp.nodist<-rep(0,n.scen)
correct.tp.ml<-rep(0,n.scen)


my.tree<-unroot(read.tree(here("analysis", "tree_files", "sim_phylo6-14.tree")))

tic()
for(s in 1:length(n.sites)){
  for(i in 1:reps){
    my.align <- simSeq(x = my.tree, l = n.sites[s])

    ols.tree <- phylo.ls(my.align)
    nodist.tree <- phylo.ls.nodist(my.align)
    ml.tree <- phylo.ML(my.align)

    correct.tp.ls[s] <- ifelse(all.equal(my.tree, ols.tree, use.edge.length = F), correct.tp.ls[s]+1, correct.tp.ls[s])
    correct.tp.nodist[s] <- ifelse(all.equal(my.tree, nodist.tree, use.edge.length = F), correct.tp.nodist[s]+1, correct.tp.nodist[s])
    correct.tp.ml[s] <- ifelse(all.equal(my.tree, ml.tree, use.edge.length = F), correct.tp.ml[s]+1, correct.tp.ml[s])

    print(i)

  }
}
toc()

save.image(here("analysis", "Fig_topol_6taxa.RData"))


# 7312.33 sec elapsed
# for 2 reps. Redoing with faster code 4/7/20
# 3621.26 sec elapsted for 2 reps after speedup. Now would take 50 hours for 100 reps
# 4/7/20 Took out three scenarios, so maybe 30 hours?
