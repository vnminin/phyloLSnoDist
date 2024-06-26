library(ape)
library(here)


# Load 4 taxa data
load(here("analysis", "Fig_kappa_4taxa_6272022.RData"))
tmp4 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,4], kappa.est[,,2], kappa.est[,,3]), plot=F)

k.norm.err4 <- array(NA,c(reps, length(kappa), length(scen)))

for(l in 1:length(scen)){
  for(k in 1:length(kappa)){
    k.norm.err4[,k,l] <- (kappa.est[,k,l] - kappa[k])/kappa[k]
  }
}

tmp4.err <- boxplot(cbind(k.norm.err4[,,1], k.norm.err4[,,4], k.norm.err4[,,2], k.norm.err4[,,3]), plot=F)

# branch length errors
all.sum.err4 <- matrix(NA, nrow=reps, ncol=(2*length(kappa)*length(scen)))

for(i in 1:reps){
  m <- 1
  for(l in c(1,4,2,3)){  # due to ULI coming last
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
load(here("analysis", "Fig_kappa_5taxa_6272022.RData"))
tmp5 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,2], kappa.est[,,3], kappa.est[,,4], kappa.est[,,5]), plot=F)

k.norm.err5 <- array(NA,c(reps, length(kappa), length(scen)))

for(l in 1:length(scen)){
  for(k in 1:length(kappa)){
    k.norm.err5[,k,l] <- (kappa.est[,k,l] - kappa[k])/kappa[k]
  }
}

tmp5.err <- boxplot(cbind(k.norm.err5[,,1], k.norm.err5[,,2], k.norm.err5[,,3], k.norm.err5[,,4], k.norm.err5[,,5]), plot=F)

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
pdf(here("analysis", "Fig_kappa.pdf"), width=14, height=16)
par(mfrow=c(3,2), mar=c(7.5,5.75,4,2.1))


# kappa straight up
plot.new()
plot.window(xlim=c(0.5,16.5),ylim=c(0,6.25),xaxs="i",yaxs="i")
rect(4.5,-1,8.5,6.25,border="gray90",col="gray90")
rect(12.5,-1,16.5,6.25,border="gray90",col="gray90")
bxp(tmp4, add=T, xaxt="n", cex.axis=2)
text(2.5, 6, "Balanced", cex=2)
text(6.5, 6, "ULI1", cex=2)
text(10.5, 6, "ULE1", cex=2)
text(14.5, 6, "ULE2", cex=2)
axis(1,line=0,at=1:16,rep(kappa, 4), cex.axis=2, las=2)
mtext("Four taxa",side=3,line=0.5, cex=3)
mtext("True kappa", side=1, line=6, cex=2)
mtext("Estimated kappa", side=2, line=3.5, cex=2)
box()

plot.new()
plot.window(xlim=c(0.5,20.5),ylim=c(0,6.25),xaxs="i",yaxs="i")
rect(4.5,-1,8.5,6.25,border="gray90",col="gray90")
rect(12.5,-1,16.5,6.25,border="gray90",col="gray90")
bxp(tmp5, add=T, xaxt="n", cex.axis=2)
text(2.5, 6, "Balanced", cex=2)
text(6.5, 6, "ULI1", cex=2)
text(10.5, 6, "ULI2", cex=2)
text(14.5, 6, "ULE1", cex=2)
text(18.5, 6, "ULE2", cex=2)
axis(1,line=0,at=1:20,rep(kappa, 5), cex.axis=2, las=2)
mtext("Five taxa",side=3,line=0.5, cex=3)
mtext("True kappa", side=1, line=6, cex=2)
box()


# kappa normalized
plot.new()
plot.window(xlim=c(0.5,16.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(4.5,-3,8.5,3,border="gray90",col="gray90")
rect(12.5,-3,16.5,3,border="gray90",col="gray90")
bxp(tmp4.err, add=T, xaxt="n", cex.axis=2)
abline(h=0, lty=2)
text(2.5, 1.75, "Balanced", cex=2)
text(6.5, 1.75, "ULI1", cex=2)
text(10.5, 1.75, "ULE1", cex=2)
text(14.5, 1.75, "ULE2", cex=2)
axis(1,line=0,at=1:16,rep(kappa, 4), cex.axis=2, las=2)
mtext("True kappa", side=1, line=6, cex=2)
mtext("Normalized Errors, kappa", side=2, line=3.5, cex=2)
box()

plot.new()
plot.window(xlim=c(0.5,20.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(4.5,-3,8.5,3,border="gray90",col="gray90")
rect(12.5,-3,16.5,3,border="gray90",col="gray90")
bxp(tmp5.err, add=T, xaxt="n", cex.axis=2)
abline(h=0, lty=2)
text(2.5, 1.75, "Balanced", cex=2)
text(6.5, 1.75, "ULI1", cex=2)
text(10.5, 1.75, "ULI2", cex=2)
text(14.5, 1.75, "ULE1", cex=2)
text(18.5, 1.75, "ULE2", cex=2)
axis(1,line=0,at=1:20,rep(kappa, 5), cex.axis=2, las=2)
mtext("True kappa", side=1, line=6, cex=2)
box()


# sum of normalized branch lengths
plot.new()
plot.window(xlim=c(0.5,32.5),ylim=c(-2,2),xaxs="i",yaxs="i")
rect(8.5,-3,16.5,3,border="gray90",col="gray90")
rect(24.5,-3,32.5,3,border="gray90",col="gray90")
boxplot(all.sum.err4,xaxt="n",col=rep(c("white","gray60"),16),add=T, cex.axis=2)
abline(h=0, lty=2)
text(4.5, 1.75, "Balanced", cex=2)
text(12.5, 1.75, "ULI1", cex=2)
text(20.5, 1.75, "ULE1", cex=2)
text(28.5, 1.75, "ULE2", cex=2)
axis(1,line=0,at=seq(1.5, 31.5, by=2),rep(kappa, 4), cex.axis=2, las=2)
mtext("True kappa", side=1, line=6, cex=2)
mtext("Normalized Errors, branch lengths", side=2, line=3.5, cex=2)
box()

plot.new()
plot.window(xlim=c(0.5,40.5),ylim=c(-5.25,5.25),xaxs="i",yaxs="i")
rect(8.5,-5.25,16.5,5.25,border="gray90",col="gray90")
rect(24.5,-5.25,32.5,5.25,border="gray90",col="gray90")
boxplot(all.sum.err5,xaxt="n",col=rep(c("white","gray60"),20),add=T, cex.axis=2)
abline(h=0, lty=2)
text(4.5, 4.6, "Balanced", cex=2)
text(12.5, 4.6, "ULI1", cex=2)
text(20.5, 4.6, "ULI2", cex=2)
text(28.5, 4.6, "ULE1", cex=2)
text(36.5, 4.6, "ULE2", cex=2)
axis(1,line=0,at=seq(1.5, 40.5, by=2),rep(kappa, 5), cex.axis=2, las=2)
mtext("True kappa", side=1, line=6, cex=2)
# legend("bottomright",fill=c("white","gray60"),c(expression(L[1]),expression(L[2])), cex=3)
legend("bottomright",fill=c("white","gray60"),c("OLS","New LS"), cex=3)
#text(38.5, -3.1, "L", vfont=c("script","bold"), cex=3)
#text(38.5, -4.2, "L", vfont=c("script","bold"), cex=3)
#text(39.5, -3.3, "1", cex=1.6)
#text(39.5, -4.4, "2", cex=1.6)
box()


dev.off()


