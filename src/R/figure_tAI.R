postscript("figures/Figure4/figure_tAI.ps", width=5, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,1))
par(mar=c(6,6,2,1))

data <- array(data, dim=c(2200, 8, 7) )
s <- rowSums(data[,6,] > 10) > 6 & rowSums(data[,8,] > 10) > 6
# require that sites have >10 mutations in all passages 
# this is looking at a suble effect, so it is necessary to filter out sites
# without mutations as they won't hold any discriminative information
# this threshold is arbitrary but conservative and the results are robust


plot( -1, -1, xlim=c(-0.45, 0.45), ylim=c(0, 3), axes=F, xlab="", ylab="")
arrows(0,0,0,3, length=0, lty=2, lwd=1, col=1)

d.rand <- density(data[s,2,1]-data[s,1,1], bw=0.075, na.rm=T)		# rand is same for all passages
a.rand <- d.rand$x > -0.45 & d.rand$x < 0.45
lines(d.rand$x[a.rand], d.rand$y[a.rand], col="#777777", lwd=3)

d.wt <- density(data[s,3,7]-data[s,1,1], bw=0.075, na.rm=T)	# refseq tAI is same for all passages
a.wt <- d.wt$x > -0.45 & d.wt$x < 0.45
lines(d.wt$x[a.wt], d.wt$y[a.wt], col="#29ABE2", lwd=3)

d.ga <- density(data[s,4,7]-data[s,1,1], bw=0.075, na.rm=T)
a.ga <- d.ga$x > -0.45 & d.ga$x < 0.45
lines(d.ga$x[a.ga], d.ga$y[a.ga], col=2, lwd=3, lty=4)

delta_tAI_rand <- data[s,2,1]-data[s,1,1]
delta_tAI_WT   <- rowMeans(data[s,3,])-data[s,1,1]
delta_tAI_GA   <- rowMeans(data[s,4,])-data[s,1,1]

boxplot( delta_tAI_rand, delta_tAI_WT, delta_tAI_GA, horizontal=T, add=T, axes=F, at=c(0.2, 0.5, 0.8), boxwex=0.25, col=c("#777777", "#29ABE2", 2))

axis(1, at=c(-0.4, -0.2, 0, 0.2, 0.4), labels=c(-0.4, -0.2, 0, 0.2, 0.4), cex.axis=1.1)
axis(2, at=c(0, 1, 2, 3), cex.axis=1.1)

mtext("Normalized frequency", side=2, line=2.5)
mtext("Relative tAI of mutations", side=1, line=2.5)

legend(0.2, 2.8, legend=c("Random", "Hsp90N", "hsp90i"), pch=15, col=c("#777777", "#29ABE2", 2), bty='n', xpd=T, cex=1.2)

pval_wt <- wilcox.test( delta_tAI_WT, delta_tAI_rand  )$p.value
pval_ga <- wilcox.test( delta_tAI_GA, delta_tAI_rand  )$p.value

text(-0.2, 2.5, paste("Hsp90N: p =", sprintf("%11.3e",pval_wt) , sep=' '), cex=0.6)
text(-0.2, 2.2, paste("hsp90i: p =", sprintf("%11.3e",pval_ga) , sep=' '), cex=0.6)


dev.off()

