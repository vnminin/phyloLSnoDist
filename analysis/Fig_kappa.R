library(ape)
library(here)


# Load 4 taxa data
load(here("analysis", "Fig_kappa_4taxa.RData"))
tmp4 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,2], kappa.est[,,3]), plot=F)

k.norm.err4 <- array(NA,c(reps, length(kappa), length(scen)))

for(l in 1:length(scen)){
  for(k in 1:length(kappa)){
    k.norm.err4[,k,l] <- (kappa.est[,k,l] - kappa[k])/kappa[k]
  }
}

tmp4.err <- boxplot(cbind(k.norm.err4[,,1], k.norm.err4[,,2], k.norm.err4[,,3]), plot=F)

# branch length errors
all.sum.err4 <- matrix(NA, nrow=reps, ncol=(2*length(kappa)*length(scen)))

for(i in 1:reps){
  m <- 1
  for(l in 1:length(scen)){
    tree.file<-paste(paste("sim_phylo4-",scen[l],sep=""),".tree",sep="")
    my.tree<-unroot(read.tree(here("analysis", "tree_files", tree.file)))
    for(k in 1:length(kappa)){
      all.sum.err4[i,m] <- sum((reg.brlen[i, , k, l] - my.tree$edge.length)/my.tree$edge.length)
      m <- m+1
      all.sum.err4[i,m] <- sum((new.brlen[i, , k, l] - my.tree$edge.length)/my.tree$edge.length)
      m <- m+1
    }
  }
}

tmp4.br <- boxplot(all.sum.err4, plot=F)

# Load 5 taxa data
load(here("analysis", "Fig_kappa_5taxa.RData"))
tmp5 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,2], kappa.est[,,3]), plot=F)

k.norm.err5 <- array(NA,c(reps, length(kappa), length(scen)))

for(l in 1:length(scen)){
  for(k in 1:length(kappa)){
    k.norm.err5[,k,l] <- (kappa.est[,k,l] - kappa[k])/kappa[k]
  }
}

tmp5.err <- boxplot(cbind(k.norm.err5[,,1], k.norm.err5[,,2], k.norm.err5[,,3]), plot=F)

# branch length errors
all.sum.err5 <- matrix(NA, nrow=reps, ncol=(2*length(kappa)*length(scen)))

for(i in 1:reps){
  m <- 1
  for(l in 1:length(scen)){
    tree.file<-paste(paste("sim_phylo5-",scen[l],sep=""),".tree",sep="")
    my.tree<-unroot(read.tree(here("analysis", "tree_files", tree.file)))
    for(k in 1:length(kappa)){
      all.sum.err5[i,m] <- sum((reg.brlen[i, , k, l] - my.tree$edge.length)/my.tree$edge.length)
      m <- m+1
      all.sum.err5[i,m] <- sum((new.brlen[i, , k, l] - my.tree$edge.length)/my.tree$edge.length)
      m <- m+1
    }
  }
}

tmp5.br <- boxplot(all.sum.err5, plot=F)

# start figure
pdf(here("analysis", "Fig_kappa.pdf"), width=14, height=18)
par(mfrow=c(3,2))


# kappa straight up
plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,6),xaxs="i",yaxs="i")
rect(4.5,-1,8.5,6,border="gray90",col="gray90")
bxp(tmp4, add=T, xaxt="n")
text(2.5, 5.75, "Balanced")
text(6.5, 5.75, "ULE1")
text(10.5, 5.75, "ULE2")
axis(1,line=0,at=1:12,rep(kappa, 3))
mtext("Four taxa",side=3,line=0.5)
mtext("True kappa", side=1, line=2.5)
mtext("Estimated kappa", side=2, line=2.5)
box()

plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,6),xaxs="i",yaxs="i")
rect(4.5,-1,8.5,6,border="gray90",col="gray90")
bxp(tmp5, add=T, xaxt="n")
text(2.5, 5.75, "Balanced")
text(6.5, 5.75, "ULE1")
text(10.5, 5.75, "ULE2")
axis(1,line=0,at=1:12,rep(kappa, 3))
mtext("Five taxa",side=3,line=0.5)
mtext("True kappa", side=1, line=2.5)
box()


# kappa normalized
plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(4.5,-3,8.5,3,border="gray90",col="gray90")
bxp(tmp4.err, add=T, xaxt="n")
text(2.5, 1.75, "Balanced")
text(6.5, 1.75, "ULE1")
text(10.5, 1.75, "ULE2")
axis(1,line=0,at=1:12,rep(kappa, 3))
mtext("True kappa", side=1, line=2.5)
mtext("Normalized Errors, kappa", side=2, line=2.5)
box()

plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(4.5,-3,8.5,3,border="gray90",col="gray90")
bxp(tmp5.err, add=T, xaxt="n")
text(2.5, 1.75, "Balanced")
text(6.5, 1.75, "ULE1")
text(10.5, 1.75, "ULE2")
axis(1,line=0,at=1:12,rep(kappa, 3))
mtext("True kappa", side=1, line=2.5)
box()


# sum of normalized branch lengths
plot.new()
plot.window(xlim=c(0.5,24.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(8.5,-3,16.5,3,border="gray90",col="gray90")
boxplot(all.sum.err4,xaxt="n",col=rep(c("white","gray60"),12),add=T)
text(4.5, 1.75, "Balanced")
text(12.5, 1.75, "ULE1")
text(20.5, 1.75, "ULE2")
axis(1,line=0,at=seq(1.5, 23.5, by=2),rep(kappa, 3))
mtext("True kappa", side=1, line=2.5)
mtext("Total sum of Normalized Errors, branch lengths", side=2, line=2.5)
legend("topleft",fill=c("white","gray60"),c(expression(L[1]),expression(L[2])))
box()

plot.new()
plot.window(xlim=c(0.5,24.5),ylim=c(-3.5,3.5),xaxs="i",yaxs="i")
rect(8.5,-4,16.5,4,border="gray90",col="gray90")
boxplot(all.sum.err5,xaxt="n",col=rep(c("white","gray60"),12),add=T)
text(4.5, 3.15, "Balanced")
text(12.5, 3.15, "ULE1")
text(20.5, 3.15, "ULE2")
axis(1,line=0,at=seq(1.5, 23.5, by=2),rep(kappa, 3))
mtext("True kappa", side=1, line=2.5)
box()


dev.off()


