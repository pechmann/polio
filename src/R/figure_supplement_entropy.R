postscript("figures/Supplement/figure_supplement_diversity.ps", width=4, height=3.5, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,2))
par(mar=c(6,6,3,1))

entwt <- array(entwt, dim=c(2200,7) )
entga <- array(entga, dim=c(2200,7) )
divwt <- array(divwt, dim=c(2200,7) )
divga <- array(divga, dim=c(2200,7) )

entwt <- entwt[rowSums(entwt)>0,]
entga <- entga[rowSums(entga)>0,]
divwt <- divwt[rowSums(divwt)>0,]
divga <- divga[rowSums(divga)>0,]

boxplot(rowMeans(entwt), rowMeans(entga), col=c("#29ABE2", 2) , ylim=c(0.2, 1.8), axes=F, at=c(0.2, 0.7), xlim=c(0, 0.9), boxwex=0.4)
axis(1, at=c(0.2, 0.7), labels=c("Hsp90N", "hsp90i"), cex.axis=1, tick=F, las=2)
axis(2, at=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5), cex.axis=1)
mtext("Mutational entropy", side=2, line=2.5)

pval = wilcox.test( rowMeans(entwt), rowMeans(entga) )$p.val
mtext(paste("p =", sprintf("%11.3e",pval) , sep=' '), side=3, line=2, cex=0.8)


boxplot(1-rowMeans(divwt), 1-rowMeans(divga), col=c("#29ABE2", 2) , ylim=c(0.2, 1), axes=F, at=c(0.2, 0.7), xlim=c(0, 0.9), boxwex=0.4)
axis(1, at=c(0.2, 0.7), labels=c("Hsp90N", "hsp90i"), cex.axis=1, tick=F, las=2)
axis(2, at=c(0.2, 0.6, 1), labels=c(0.2, 0.6, 1), cex.axis=1)
mtext("Mutational diversity (1-D)", side=2, line=2.5)

pval = wilcox.test( rowMeans(divwt), rowMeans(divga) )$p.val
mtext(paste("p =", sprintf("%11.3e",pval) , sep=' '), side=3, line=2, cex=0.8)



dev.off()



