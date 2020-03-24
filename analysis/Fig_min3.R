library(phyloLSnoDist)
library(robustDist)

####################################################
# 1 row for each of 5, 6, 7 taxa
#

min.errors5<-matrix(NA,nrow=reps,ncol=30)

load(here("analysis","Fig_min3_5taxa.RData"))
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo5-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    min.errors5[k,((6*(l-1)) + c(1,3,5))]<-reg.errors[k, min3, l]
    min.errors5[k,((6*(l-1)) + c(2,4,6))]<-new.errors[k, min3, l]
  }
}



min.errors6<-matrix(NA,nrow=reps,ncol=30)

load(here("analysis","Fig_min3_6taxa.RData"))
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo6-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    min.errors6[k,((6*(l-1)) + c(1,3,5))]<-reg.errors[k, min3, l]
    min.errors6[k,((6*(l-1)) + c(2,4,6))]<-new.errors[k, min3, l]
  }
}




min.errors7<-matrix(NA,nrow=reps,ncol=30)

load(here("analysis","Fig_min3_7taxa.RData"))
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo7-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    min.errors7[k,((6*(l-1)) + c(1,3,5))]<-reg.errors[k, min3, l]
    min.errors7[k,((6*(l-1)) + c(2,4,6))]<-new.errors[k, min3, l]
  }
}


pdf(here("analysis", "Fig_min3.pdf"),height=8,width=8)
par(mfrow=c(3,1),oma=c(4,3,2,0),mar=c(2,3,0,1))

boxplot(min.errors5,xaxt="n")
rect(6.5,-10,12.5,10,border="gray95",col="gray95")
rect(18.5,-10,24.5,10,border="gray95",col="gray95")
boxplot(min.errors5,xaxt="n",col=rep(c("white","gray60"),15),add=T)
abline(h=0,col="gray40",lty=2)
legend("topleft",fill=c("white","gray60"),c(expression(L[1]),expression(L[2])))
mtext("Normalized Error", side=2,line=2.5)
mtext("Five taxa", side=2, line=4,cex=1.25)
box()


boxplot(min.errors6,xaxt="n")
rect(6.5,-10,12.5,10,border="gray95",col="gray95")
rect(18.5,-10,24.5,10,border="gray95",col="gray95")
boxplot(min.errors6,xaxt="n",col=rep(c("white","gray60"),15),add=T)
abline(h=0,col="gray40",lty=2)
mtext("Normalized Error", side=2,line=2.5)
mtext("Six taxa", side=2, line=4,cex=1.25)
box()

boxplot(min.errors7,xaxt="n")
rect(6.5,-10,12.5,10,border="gray95",col="gray95")
rect(18.5,-10,24.5,10,border="gray95",col="gray95")
boxplot(min.errors7,xaxt="n",col=rep(c("white","gray60"),15),add=T)
axis(1,seq(1.5,29.5,by=2),label=rep(c(1,2,3),5))
abline(h=0,col="gray40",lty=2)
mtext("Normalized Error", side=2,line=2.5)
mtext("Seven taxa", side=2, line=4,cex=1.25)
mtext(c("Balanced","ULI1","ULI2","ULE1","ULE2"), side=1,line=3,at=seq(3.5,29.5,by=6))
box()
dev.off()


