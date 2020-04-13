library(here)
source(here("analysis", "LCO_generator_all_n.R"))

load(here("analysis","Fig_topol_4taxa.RData"))
# 4/9/20 These indices are based on having run 6 scenarios and only needing 4. Redo if needed.
results.4<-c(correct.tp.ls[1],correct.tp.nodist[1],correct.tp.ml[1],
             correct.tp.ls[2],correct.tp.nodist[2],correct.tp.ml[2],
             correct.tp.ls[3],correct.tp.nodist[3],correct.tp.ml[3],
             correct.tp.ls[4],correct.tp.nodist[4],correct.tp.ml[4])


# Run the LCO generator for a sample size of reps (100 currently)
lco.bounds <- LCO.CI(reps,0.95,3)

res4.box<-boxplot(t(as.matrix(results.4/reps)),plot=F)
res4.box$stats <- sapply(res4.box$stats[1,],function(x) c(lco.bounds[(x*reps),2],x,x,x,lco.bounds[(x*reps),3]))
#res4.box$stats<-sapply(t(as.matrix(results.4/reps)),function(x) c(x-1.96*sqrt(x*(1-x))/sqrt(reps),x,x,x,x+1.96*sqrt(x*(1-x))/sqrt(reps)))
n.sites4 <- n.sites[c(1,2,4,6)] # 4/9/20 change indeces if re-run

load(here("analysis","Fig_topol_5taxa.RData"))
results.5<-c(correct.tp.ls[1],correct.tp.nodist[1],correct.tp.ml[1],
             correct.tp.ls[2],correct.tp.nodist[2],correct.tp.ml[2],
             correct.tp.ls[3],correct.tp.nodist[3],correct.tp.ml[3],
             correct.tp.ls[4],correct.tp.nodist[4],correct.tp.ml[4])

res5.box<-boxplot(t(as.matrix(results.5/reps)),plot=F)
res5.box$stats <- sapply(res5.box$stats[1,],function(x) c(lco.bounds[(x*reps),2],x,x,x,lco.bounds[(x*reps),3]))
# res5.box$stats<-sapply(t(as.matrix(results.5/reps)),function(x) c(x-1.96*sqrt(x*(1-x))/sqrt(reps),x,x,x,x+1.96*sqrt(x*(1-x))/sqrt(reps)))
n.sites5 <- n.sites


load(here("analysis","Fig_topol_6taxa.RData"))
results.6<-c(correct.tp.ls[1],correct.tp.nodist[1],correct.tp.ml[1],
             correct.tp.ls[2],correct.tp.nodist[2],correct.tp.ml[2],
             correct.tp.ls[3],correct.tp.nodist[3],correct.tp.ml[3],
             correct.tp.ls[4],correct.tp.nodist[4],correct.tp.ml[4])

res6.box<-boxplot(t(as.matrix(results.6/reps)),plot=F)
res6.box$stats <- sapply(res6.box$stats[1,],function(x) c(lco.bounds[(x*reps),2],x,x,x,lco.bounds[(x*reps),3]))
# res5.box$stats<-sapply(t(as.matrix(results.5/reps)),function(x) c(x-1.96*sqrt(x*(1-x))/sqrt(reps),x,x,x,x+1.96*sqrt(x*(1-x))/sqrt(reps)))
n.sites6 <- n.sites


# 4, 5 and 6 taxa together

pdf(here("analysis", "Fig_topol.pdf"),height=12,width=8)
par(mfrow=c(3,1),oma=c(4,2,5,0),mar=c(3,2,5,2))

plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,1.05),xaxs="i",yaxs="i")
rect(3.5,-1,6.5,1.05,border="gray90",col="gray90")
rect(9.5,-1,12.5,1.05,border="gray90",col="gray90")
bxp(res4.box,outline=F,medpch=rep(c(16,17,18),4),medcex=1.7,medlty="blank",whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.12,axes=F,add=T)
axis(2,at=seq(0,1,by=0.25),label=c(0,0.25,0.5,0.75,1))
axis(1,at=seq(2,13,by=3),label=n.sites4)
mtext("Percent correct",line=2.5,side=2)
mtext("Sequence length (nucleotides)",side=1,line=2.5)
mtext("Four taxa",cex=1.5,line=1.25)
legend("topleft",pch=c(16,17,18),c("LS","No Dist","MLE"))
box()



######################
plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,1.05),xaxs="i",yaxs="i")
rect(3.5,-1,6.5,1.05,border="gray90",col="gray90")
rect(9.5,-1,12.5,1.05,border="gray90",col="gray90")
bxp(res5.box,outline=F,medpch=rep(c(16,17,18),4),medcex=1.7,medlty="blank",whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.12,axes=F,add=T)
axis(2,at=seq(0,1,by=0.25),label=c(0,0.25,0.5,0.75,1))
axis(1,at=seq(2,13,by=3),label=n.sites5)
mtext("Percent correct",line=2.5,side=2)
mtext("Sequence length (nucleotides)",side=1,line=2.5)
mtext("Five taxa",cex=1.5,line=1.25)
box()



######################
plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,1.05),xaxs="i",yaxs="i")
rect(3.5,-1,6.5,1.05,border="gray90",col="gray90")
rect(9.5,-1,12.5,1.05,border="gray90",col="gray90")
bxp(res6.box,outline=F,medpch=rep(c(16,17,18),4),medcex=1.7,medlty="blank",whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.12,axes=F,add=T)
axis(2,at=seq(0,1,by=0.25),label=c(0,0.25,0.5,0.75,1))
axis(1,at=seq(2,13,by=3),label=n.sites6)
mtext("Percent correct",line=2.5,side=2)
mtext("Sequence length (nucleotides)",side=1,line=2.5)
mtext("Six taxa",cex=1.5,line=1.25)
box()


dev.off()






#######################3
# 4 taxa results only
pdf(here("analysis", "Fig_topol.pdf"),height=6,width=12)
#par(mfrow=c(1,2),oma=c(3,2,1,0),mar=c(1,2,2,2))

plot.new()
plot.window(xlim=c(0.5,18.5),ylim=c(0,1),xaxs="i",yaxs="i")
rect(3.5,-1,6.5,1,border="gray90",col="gray90")
rect(9.5,-1,12.5,1,border="gray90",col="gray90")
rect(15.5,-1,18.5,1,border="gray90",col="gray90")
bxp(res4.box,outline=F,medpch=rep(c(16,17,18),6),medcex=1.6,whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.12,axes=F,add=T)
axis(2,at=seq(0,1,by=0.25),label=c(0,0.25,0.5,0.75,1))
axis(1,at=seq(2,17,by=3),label=n.sites)
mtext("Percent correct",line=2.5,side=2)
mtext("Sequence length (nucleotides)",side=1,line=2.5)
mtext("Four taxa",cex=1.5,line=1.25)
legend("topleft",pch=c(16,17,18),c("LS","No Dist","MLE"))
box()
dev.off()

