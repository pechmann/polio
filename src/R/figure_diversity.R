postscript("figures/Figure2/figure_diversity.ps", width=2.2, height=3.5, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,1))
par(mar=c(6,6,3,1))

wt <- array(wt, dim=c(2200,7) )
ga <- array(ga, dim=c(2200,7) )

wt <- wt[rowSums(wt)>0,]
ga <- ga[rowSums(ga)>0,]

boxplot(1-rowMeans(wt), 1-rowMeans(ga), col=c("#29ABE2", 2) , ylim=c(0.2, 1), axes=F, at=c(0.2, 0.7), xlim=c(0, 0.9), boxwex=0.4)
axis(1, at=c(0.2, 0.7), labels=c("Hsp90N", "hsp90i"), cex.axis=1, tick=F, las=2)
axis(2, at=c(0.2, 0.6, 1), labels=c(0.2, 0.6, 1), cex.axis=1)
mtext("Mutational diversity (1-D)", side=2, line=2.5)

pval = wilcox.test( rowMeans(wt), rowMeans(ga) )$p.val
mtext(paste("p =", sprintf("%11.3e",pval) , sep=' '), side=3, line=2)

dev.off()



