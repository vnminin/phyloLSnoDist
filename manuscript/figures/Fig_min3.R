library(phyloLSnoDist)

####################################################
# 1 row for each of 5, 6, 7 taxa
#


load(here("analysis","Fig_min3_5taxa.RData"))
min.errors5<-matrix(NA,nrow=reps,ncol=30)
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




load(here("analysis","Fig_min3_6taxa.RData"))
min.errors6<-matrix(NA,nrow=reps,ncol=30)
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





load(here("analysis","Fig_min3_7taxa.RData"))
min.errors7<-matrix(NA,nrow=reps,ncol=30)
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
rect(6.5,-10,12.5,10,border="gray90",col="gray90")
rect(18.5,-10,24.5,10,border="gray90",col="gray90")
boxplot(min.errors5,xaxt="n",col=rep(c("white","gray60"),15),add=T)
abline(h=0,col="gray40",lty=2)
# legend("topleft",fill=c("white","gray60"),c(expression(L[1]),expression(L[2])))
legend("topleft",fill=c("white","gray60"),c("OLS", "New LS"), cex=2)
#text(1.75, 0.325, "L", vfont=c("script","bold"), cex=2.25)
#text(1.75, 0.245, "L", vfont=c("script","bold"), cex=2.25)
#text(2.2, 0.305, "1", cex=1.2)
#text(2.2, 0.225, "2", cex=1.2)
mtext("Normalized Error", side=2,line=2.5)
mtext("Five taxa", side=2, line=4,cex=1.25)
box()


boxplot(min.errors6,xaxt="n")
rect(6.5,-10,12.5,10,border="gray90",col="gray90")
rect(18.5,-10,24.5,10,border="gray90",col="gray90")
boxplot(min.errors6,xaxt="n",col=rep(c("white","gray60"),15),add=T)
abline(h=0,col="gray40",lty=2)
mtext("Normalized Error", side=2,line=2.5)
mtext("Six taxa", side=2, line=4,cex=1.25)
box()

boxplot(min.errors7,xaxt="n")
rect(6.5,-10,12.5,10,border="gray90",col="gray90")
rect(18.5,-10,24.5,10,border="gray90",col="gray90")
boxplot(min.errors7,xaxt="n",col=rep(c("white","gray60"),15),add=T)
axis(1,seq(1.5,29.5,by=2),label=rep(c(1,2,3),5))
abline(h=0,col="gray40",lty=2)
mtext("Normalized Error", side=2,line=2.5)
mtext("Seven taxa", side=2, line=4,cex=1.25)
mtext(c("Balanced","ULI1","ULI2","ULE1","ULE2"), side=1,line=3,at=seq(3.5,29.5,by=6))
box()
dev.off()




#################################
# Efficiency

load(here("analysis","Fig_min3_5taxa.RData"))
tmp.reg <- matrix(NA, nrow=reps, ncol=3)
tmp.new <- matrix(NA, nrow=reps, ncol=3)

min3.eff5 <- rep(NA, 15)
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo5-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    tmp.reg[k, ] <- reg.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
    tmp.new[k, ] <- new.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
  }
  mse.reg <- (apply(tmp.reg, 2, mean) - true[min3])^2 + apply(tmp.reg, 2, var)
  mse.new <- (apply(tmp.new, 2, mean) - true[min3])^2 + apply(tmp.new, 2, var)
  min3.eff5[3*(l-1)+c(1,2,3)] <- mse.new/mse.reg
}



load(here("analysis","Fig_min3_6taxa.RData"))
tmp.reg <- matrix(NA, nrow=reps, ncol=3)
tmp.new <- matrix(NA, nrow=reps, ncol=3)

min3.eff6 <- rep(NA, 15)
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo6-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    tmp.reg[k, ] <- reg.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
    tmp.new[k, ] <- new.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
  }
  mse.reg <- (apply(tmp.reg, 2, mean) - true[min3])^2 + apply(tmp.reg, 2, var)
  mse.new <- (apply(tmp.new, 2, mean) - true[min3])^2 + apply(tmp.new, 2, var)
  min3.eff6[3*(l-1)+c(1,2,3)] <- mse.new/mse.reg
}



load(here("analysis","Fig_min3_7taxa.RData"))
tmp.reg <- matrix(NA, nrow=reps, ncol=3)
tmp.new <- matrix(NA, nrow=reps, ncol=3)

min3.eff7 <- rep(NA, 15)
for(l in 1:n.scen){
  tree.file<-paste(paste("sim_phylo7-",scen[l],sep=""),".tree",sep="")
  my.tree<-read.tree(here("analysis", "tree_files", tree.file))
  true<-unroot(my.tree)$edge.length
  min3<-order(unroot(my.tree)$edge.length)[1:3]

  for(k in 1:reps){
    tmp.reg[k, ] <- reg.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
    tmp.new[k, ] <- new.errors[k, min3, l] * sqrt(true[min3]) + true[min3]
  }
  mse.reg <- (apply(tmp.reg, 2, mean) - true[min3])^2 + apply(tmp.reg, 2, var)
  mse.new <- (apply(tmp.new, 2, mean) - true[min3])^2 + apply(tmp.new, 2, var)
  min3.eff7[3*(l-1)+c(1,2,3)] <- mse.new/mse.reg
}



pdf(here("analysis", "Fig_min3_eff1.pdf"),height=8,width=8)
par(mfrow=c(3,1),oma=c(4,3,2,0),mar=c(2,3,0,1))

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
points(min3.eff5~c(1:15), pch=19, ylim=c(0, 2), cex=2)
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Five taxa", side=2, line=4,cex=1.25)
box()

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
points(min3.eff6~c(1:15), pch=19, ylim=c(0, 2), cex=2)
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Six taxa", side=2, line=4,cex=1.25)
box()

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
axis(1,c(1:15),label=rep(c(1,2,3),5))
points(min3.eff7~c(1:15), pch=19, ylim=c(0, 2), cex=2)
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Seven taxa", side=2, line=4,cex=1.25)
mtext(c("Balanced","ULI1","ULI2","ULE1","ULE2"), side=1,line=3,at=c(2,5,8,11,14))
box()

dev.off()



pdf(here("analysis", "Fig_min3_eff2.pdf"),height=8,width=8)
par(mfrow=c(3,1),oma=c(4,3,2,0),mar=c(2,3,0,1))

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
points(min3.eff5~c(1:15), pch=19, ylim=c(0, 2), cex=2, col=ifelse(min3.eff5>1, "gray60", "black"))
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Five taxa", side=2, line=4,cex=1.25)
box()

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
points(min3.eff6~c(1:15), pch=19, ylim=c(0, 2), cex=2, col=ifelse(min3.eff6>1, "gray60", "black"))
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Six taxa", side=2, line=4,cex=1.25)
box()

plot(NULL, xlim=c(0.5,15.5), ylim=c(0,2), xlab="", ylab="", xaxt="n")
rect(3.5,-10,6.5,10,border="gray90",col="gray90")
rect(9.5,-10,12.5,10,border="gray90",col="gray90")
axis(1,c(1:15),label=rep(c(1,2,3),5))
points(min3.eff7~c(1:15), pch=19, ylim=c(0, 2), cex=2, col=ifelse(min3.eff7>1, "gray60", "black"))
abline(h=1,col="gray40",lty=2)
mtext(expression(MSE(L[2])/MSE(L[1])), side=2,line=2.25)
mtext("Seven taxa", side=2, line=4,cex=1.25)
mtext(c("Balanced","ULI1","ULI2","ULE1","ULE2"), side=1,line=3,at=c(2,5,8,11,14))
box()

dev.off()


