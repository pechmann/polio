postscript("figures/Figure3/figure_sites_sel_properties.ps", width=4, height=4, paper="special", horizontal=F)

aggwt <- as.matrix(agg_wt)
aggga <- as.matrix(agg_ga)
stabwt <- as.matrix(stab_wt)
stabga <- as.matrix(stab_ga)
hydwt <- as.matrix(hyd_wt)
hydga <- as.matrix(hyd_ga)
rand <- as.matrix(rand)


sel.wt.agg <- aggwt[,2] < 0.05 & aggwt[,3] > 1 & aggwt[,6] > 2
sel.ga.agg <- aggga[,2] < 0.05 & aggga[,3] < 1 & aggga[,6] > 2
sel.wt.stab <- stabwt[,2] < 0.05 & stabwt[,3] > 1 & stabwt[,5] > 1
sel.ga.stab <- stabga[,2] < 0.05 & stabga[,3] < 1 & stabga[,5] > 1
sel.wt.hyd <- hydwt[,2] < 0.05 & hydwt[,3] > 1 & hydwt[,4] > 2
sel.ga.hyd <- hydga[,2] < 0.05 & hydga[,3] < 1 & hydga[,4] > 2


par(mfrow=c(1,2))
par(mar=c(6,6,2,1))



boxplot(rand[,3], stabwt[sel.wt.stab,6], stabga[sel.ga.stab,6], col=c("#999999", "#29ABE2","red"), ylim=c(-20, 60), xlab="", ylab="", axes=F,  boxwex=0.4, at=c(0.2, 0.8, 1.3), xlim=c(0, 1.5), range=0.5 )
axis(2, at=c(0, 30, 60), cex.axis=1)
mtext(expression(paste(Delta, "AP")), side=2, line=2.5)
pval_agg  <- wilcox.test(stabwt[sel.wt.stab,6], stabga[sel.ga.stab,6])$p.value
text(0.2, 50, paste("p =", sprintf("%11.2e",pval_agg) , sep=' '), cex=0.6)


boxplot(rand[,1], stabwt[sel.wt.stab,4], stabga[sel.ga.stab,4], col=c("#999999", "#29ABE2","red"), ylim=c(-10, 10), xlab="", ylab="", axes=F,  boxwex=0.4, at=c(0.2, 0.8, 1.3), xlim=c(0, 1.5), range=0.5 )
axis(2, at=c(-6, -3, 0, 3, 6) )
mtext(expression(paste(Delta, "H")), side=2, line=2.5)
pval_hyd  <- wilcox.test(stabwt[sel.wt.stab,4], stabga[sel.ga.stab,4])$p.value
text(0.2, 1, paste("p =", sprintf("%11.2e",pval_hyd) , sep=' '), cex=0.6)



dev.off()
