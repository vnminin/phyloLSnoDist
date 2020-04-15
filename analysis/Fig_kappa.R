library(here)


# Load 4 taxa data
load(here("analysis", "Fig_kappa_4taxa.RData"))
tmp4 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,2], kappa.est[,,3]))

# Load 5 taxa data
load(here("analysis", "Fig_kappa_5taxa.RData"))
tmp5 <- boxplot(cbind(kappa.est[,,1], kappa.est[,,2], kappa.est[,,3]))


# start figure
pdf(here("analysis", "Fig_kappa.pdf"), width=14, height=6)
par(mfrow=c(1,2))

plot.new()
plot.window(xlim=c(0.5,12.5),ylim=c(0,6),xaxs="i",yaxs="i")
rect(4.5,-1,8.5,6,border="gray90",col="gray90")
bxp(tmp4, add=T, xaxt="n")
text(2.5, 5.75, "Balanced")
text(6.5, 5.75, "ULE1")
text(10.5, 5.75, "ULE2")
axis(1,line=0,at=1:12,rep(kappa, 3))
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
mtext("True kappa", side=1, line=2.5)
box()

dev.off()
