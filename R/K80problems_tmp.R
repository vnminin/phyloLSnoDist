library(phyloLSnoDist)

set.seed(1)
my.tree<-read.tree(here("analysis", "tree_files", "sim_phylo5-1.tree"))
my.align <- simSeq(x = my.tree, l = 1000, Q = c(1, 4, 1, 1, 4, 1))

ts.counts <- dist.dna(as.DNAbin(my.align), model="TS")
tv.counts <- dist.dna(as.DNAbin(my.align), model="TV")

ts.counts / tv.counts   # These will be around 2, NOT 4 as they should be...
#t2       t4       t3       t5
#t4 2.172414
#t3 2.079365 2.031250
#t5 2.030769 2.015152 9.000000
#t1 1.915254 1.982759 2.365385 2.365385

new.k80 <- new.ls.fit.K80(my.tree, my.align)
new.k80$par.est[8] # This will be close to 4 as it should be
#p8
#4.819845
