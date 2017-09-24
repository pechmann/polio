postscript("figures/Figure2/figure_dnds.ps", width=3, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mar=c(6,6,4,1))
par(mfrow=c(1,1))

wt <- array(wt, dim=c(7,2) )
ga <- array(ga, dim=c(7,2) )

plot(-1, -1, t='p',  xlim=c(0.6, 2.4), ylim=c(0, 0.002), axes=F, xlab="", ylab="")
arrows(0.75, median(wt[,1]), 0.95, median(wt[,1]), length=0, lwd=2)
arrows(1.15, median(ga[,1]), 1.35, median(ga[,1]), length=0, lwd=2)
arrows(1.75, median(wt[,2]), 1.95, median(wt[,2]), length=0, lwd=2)
arrows(2.15, median(ga[,2]), 2.35, median(ga[,2]), length=0, lwd=2)

points(rep(0.8,7)+runif(7,0, 0.1), wt[,1], col="#29ABE2", pch=16, cex=1.2)
points(rep(1.2,7)+runif(7,0, 0.1), ga[,1], col=2, pch=16, cex=1.2)

points(rep(1.8,7)+runif(7,0, 0.1), wt[,2], pch=16, cex=1.2, col="#29ABE2")
points(rep(2.2,7)+runif(7,0, 0.1), ga[,2], col=2, pch=16, cex=1.2)

axis(1, at=c(1,2), labels=c("dS", "dN"), cex.axis=1.1, line=0, tick=F)
#axis(2, at=c(0, 0.005, 0.001, 0.0015), labels=c(0, format(0.0005, scientific=F), 0.001, 0.0015), cex.axis=1.1)
axis(2, at=c(0, 0.001, 0.002), labels=c(0, 0.001, 0.002), cex.axis=1.1)


mtext("Evolutionary rate", side=2, line=3, cex=1.3)

legend("topright", legend=c("Hsp90N", "hsp90i"), col=c("#29ABE2", 2), pch=16, cex=1, xpd=T)

pval1 = wilcox.test( wt[,1], ga[,1] )$p.val
mtext(paste("p =", sprintf("%11.3e",pval1) , sep=' '), side=3, line=2, cex=0.8)

pval2 = wilcox.test( wt[,2], ga[,2] )$p.val
mtext(paste("p =", sprintf("%11.3e",pval2) , sep=' '), side=3, line=1, cex=0.8)

dev.off()

