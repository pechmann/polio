postscript("figures/Supplement/figure_supplement_coverage.ps", width=6, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,1))
par(mar=c(6,8,2,1))

wt <- array(wt, dim=c(2218,7) )
ga <- array(ga, dim=c(2218,7) )

plot( log10(rowMeans( wt )), t='l', col="#29ABE2", lwd=1, axes=F, xlab="", ylab="", ylim=c(1.8, 5.5))

lines( log10(rowMeans( ga )), lwd=1, col=2)

axis(1, cex.axis=1.2)
axis(2, at=c(2,3, 4, 5), labels=c("100", "1000", "10,000", "100,000"), las=2, cex.axis=1.2)

mtext("Number of reads", side=2, line=5, cex=1.3)
mtext("Poliovirus coding sequence", side=1, line=2.5, cex=1.2)

legend(10, 3.5, legend=c("Hsp90N", "hsp90i"), lty=1, lwd=2, col=c("#29ABE2", 2), bty='n', xpd=T, cex=1.2)

dev.off()

