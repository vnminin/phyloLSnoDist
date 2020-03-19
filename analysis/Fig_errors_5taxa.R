library(robustDist)
setwd("C:/Users/pchi01/Dropbox/Villanova/Research/LeastSquares")


reps <- 1000
my.tree<-read.newick("sim_phylo5-9.tree")
n.br<-length(unroot(my.tree)$edge.length)

new.brlen<-matrix(NA,nrow=reps,ncol=n.br)
reg.brlen<-matrix(NA,nrow=reps,ncol=n.br)
new.errors<-matrix(NA,nrow=reps,ncol=n.br)
reg.errors<-matrix(NA,nrow=reps,ncol=n.br)
counts <- rep(NA, reps)



for(i in 1:reps){
  my.align <- as.character(simSeq(my.tree, 2000))

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
  counts[i] <- optim.out$count


  true<-unroot(my.tree)$edge.length

  reg.errors[i,] <- (reg.brlen[i,] - true)/sqrt(true)
  new.errors[i,] <- (new.brlen[i,] - true)/sqrt(true)

}

# re-order to match exam figure

br.names<-rep(NA,n.br)
n.tips<-length(my.tree$tip.label)

for(k in 1:n.br){			# to put node/tip labels on the MSE table, and fill it in
  nodes<-unroot(my.tree)$edge[k,]
  first<-ifelse(nodes[1]<=n.tips,my.tree$tip.label[nodes[1]],paste("n",nodes[1],sep=""))
  second<-ifelse(nodes[2]<=n.tips,my.tree$tip.label[nodes[2]],paste("n",nodes[2],sep=""))
  br.names[k]<-paste(first,second,sep="-")
}


all.errors<-matrix(NA,nrow=reps,ncol=14)

for(i in 1:reps){
  all.errors[i,c(1,3,5,7,9,11,13)]<-(reg.brlen[i,]-true)/true
  all.errors[i,c(2,4,6,8,10,12,14)]<-(new.brlen[i,]-true)/true
}




br.names<-br.names[c(2,3,1,7,4,6,5)]
errors.reord<-as.data.frame(all.errors[,c(3,4,5,6,1,2,13,14,7,8,11,12,9,10)])


pdf("Fig_errors_5taxa.pdf",height=6.5,width=6.5)
par(mfrow=c(2,1),oma=c(3,2.5,1.5,0),mar=c(1,1,0,1))
plot(unroot(my.tree),type="u",show.tip.label=F,y.lim=c(-0.05,0.5))
mtext("Five taxa, Scenario ULE2",side=3,line=0.5)
nodelabels(my.tree$tip.label[1],1,bg=NULL,frame="none",font=2,adj=c(1.25,1.25))
nodelabels(my.tree$tip.label[2],2,bg=NULL,frame="none",font=2,adj=c(0.25,1.25))
nodelabels(my.tree$tip.label[3],3,bg=NULL,frame="none",font=2,adj=c(-0.25,0.4))
nodelabels(my.tree$tip.label[4],4,bg=NULL,frame="none",font=2,adj=c(0.4,-0.5))
nodelabels(my.tree$tip.label[5],5,bg=NULL,frame="none",font=2,adj=c(1.25,0))
nodelabels("n6",6,bg="white")
nodelabels("n7",7,bg="white")
nodelabels("n8",8,bg="white",adj=c(-1,2.5))
arrows(x0=0.5,y0=0.35,x1=0.472,y1=0.39,length=0.1)
box()

boxplot(errors.reord, col=c("white", "gray"),xlab="",ylab="",xaxt="n")
legend("topleft",c(expression(L[1]),expression(L[2])),fill=c("white","gray"))
axis(1,line=0,at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5),br.names)
mtext("Normalized errors",side=2,line=2)
abline(h=0,lty=3)
mtext("branch",side=1,line=2.5)
dev.off()

