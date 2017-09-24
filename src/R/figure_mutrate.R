postscript("figures/Figure2/figure_mutrate.ps", width=2.2, height=3.5, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,1))
par(mar=c(6,6,3,1))

data <- as.matrix(data)

boxplot(data[,1], data[,2], col=c("#29ABE2", 2) , ylim=c(0.00005, 0.00009), axes=F, at=c(0.2, 0.7), xlim=c(0, 0.9), boxwex=0.4)
axis(1, at=c(0.2, 0.7), labels=c("Hsp90N", "hsp90i"), cex.axis=1, tick=F, las=2)
axis(2, cex.axis=1)
mtext("Mutation rate", side=2, line=2.5)

pval = wilcox.test( data[,1], data[,2] )$p.val
mtext(paste("p =", sprintf("%11.3e",pval) , sep=' '), side=3, line=2)

dev.off()



