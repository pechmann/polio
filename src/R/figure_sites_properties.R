postscript("figures/Figure3/figure_sites_properties.ps", width=8, height=4, paper="special", horizontal=F)

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


par(mfrow=c(1,3))
par(mar=c(6,6,2,1))


boxplot(rand[,2], stabwt[sel.wt.stab,5], stabga[sel.ga.stab,5], aggwt[sel.wt.agg,5], aggga[sel.ga.agg,5], hydwt[sel.wt.hyd,5], hydga[sel.ga.hyd,5], col=c("#999999", "#29ABE2","red", "#29ABE2","red", "#29ABE2","red"), ylim=c(-2, 6.5), xlab="", ylab="", axes=F,  boxwex=0.4, at=c(0.2, 0.8, 1.3, 1.9, 2.4, 3, 3.5), xlim=c(0, 3.7), range=0.5 )
axis(2, at=c(-2, 0, 2, 4), cex.axis=1)
mtext(expression(paste(Delta, Delta, "G")), side=2, line=2.5)

pval_stab_stab <- wilcox.test( stabwt[sel.wt.stab,5], stabga[sel.ga.stab,5]  )$p.val 
pval_stab_agg  <- wilcox.test( aggwt[sel.wt.agg,5], aggga[sel.ga.agg,5]  )$p.val 
pval_stab_hyd <- wilcox.test(  hydwt[sel.wt.hyd,5], hydga[sel.ga.hyd,5]  )$p.val 

text(1, 2, paste("s/s p =", sprintf("%11.3e",pval_stab_stab) , sep=' '), cex=0.6)
text(1, 1, paste("s/a p =", sprintf("%11.3e",pval_stab_agg) , sep=' '), cex=0.6)
text(1, 0, paste("s/h p =", sprintf("%11.3e",pval_stab_hyd) , sep=' '), cex=0.6)

pval_stab_randagg  <- wilcox.test( rand[,2], aggwt[sel.wt.agg,5]   )$p.val 
pval_stab_randhyd <- wilcox.test(  rand[,2], hydwt[sel.wt.hyd,5]  )$p.val 
text(2, 1, paste("r/a p =", sprintf("%11.3e",pval_stab_randagg) , sep=' '), cex=0.6)
text(2, 0, paste("r/h p =", sprintf("%11.3e",pval_stab_randhyd) , sep=' '), cex=0.6)


boxplot(rand[,3], stabwt[sel.wt.stab,6], stabga[sel.ga.stab,6], aggwt[sel.wt.agg,6], aggga[sel.ga.agg,6], hydwt[sel.wt.hyd,6], hydga[sel.ga.hyd,6], col=c("#999999", "#29ABE2","red", "#29ABE2","red", "#29ABE2","red"), ylim=c(-20, 300), xlab="", ylab="", axes=F,  boxwex=0.4, at=c(0.2, 0.8, 1.3, 1.9, 2.4, 3, 3.5), xlim=c(0, 3.7), range=0.5 )
axis(2, at=c(0, 100, 200, 300), cex.axis=1)
mtext(expression(paste(Delta, "AP")), side=2, line=2.5)

pval_agg_stab <- wilcox.test( stabwt[sel.wt.stab,6], stabga[sel.ga.stab,6] )$p.val
pval_agg_agg <- wilcox.test( aggwt[sel.wt.agg,6], aggga[sel.ga.agg,6] )$p.val
pval_agg_hyd <- wilcox.test(  hydwt[sel.wt.hyd,6], hydga[sel.ga.hyd,6] )$p.val

text(1, 100, paste("a/s p =", sprintf("%11.3e",pval_agg_stab) , sep=' '), cex=0.6)
text(1, 90, paste("a/a p =", sprintf("%11.3e",pval_agg_agg) , sep=' '), cex=0.6)
text(1, 80, paste("a/h p =", sprintf("%11.3e",pval_agg_hyd) , sep=' '), cex=0.6)

pval_agg_randaggp  <- wilcox.test( rand[,3],  aggwt[sel.wt.agg,6]  )$p.val 
pval_agg_randaggn <-  wilcox.test( rand[,3],  aggga[sel.ga.agg,6]  )$p.val 
text(2.5, 100, paste("r/a p =", sprintf("%11.3e",pval_agg_randaggp) , sep=' '), cex=0.6)
text(2.5, 90, paste("r/h p =", sprintf("%11.3e",pval_agg_randaggn) , sep=' '), cex=0.6)


boxplot(rand[,1], stabwt[sel.wt.stab,4], stabga[sel.ga.stab,4], aggwt[sel.wt.agg,4], aggga[sel.ga.agg,4], hydwt[sel.wt.hyd,4], hydga[sel.ga.hyd,4], col=c("#999999", "#29ABE2","red", "#29ABE2","red", "#29ABE2","red"), ylim=c(-10, 10), xlab="", ylab="", axes=F,  boxwex=0.4, at=c(0.2, 0.8, 1.3, 1.9, 2.4, 3, 3.5), xlim=c(0, 3.7), range=0.5 )
axis(2, at=c(-6, -3, 0, 3, 6) )
mtext(expression(paste(Delta, "H")), side=2, line=2.5)

pval_hyd_stab <- wilcox.test( stabwt[sel.wt.stab,4], stabga[sel.ga.stab,4] )$p.val
pval_hyd_agg  <- wilcox.test( aggwt[sel.wt.agg,4], aggga[sel.ga.agg,4] )$p.val
pval_hyd_hyd  <- wilcox.test( hydwt[sel.wt.hyd,4], hydga[sel.ga.hyd,4] )$p.val

text(1, 2, paste("h/s p =", sprintf("%11.3e",pval_hyd_stab) , sep=' '), cex=0.6)
text(1, 1, paste("h/a p =", sprintf("%11.3e",pval_hyd_agg) , sep=' '), cex=0.6)
text(1, 0, paste("h/h p =", sprintf("%11.3e",pval_hyd_hyd) , sep=' '), cex=0.6)

pval_hyd_stabp <- wilcox.test( rand[,1], stabwt[sel.wt.stab,4] )$p.val
pval_hyd_stabn <- wilcox.test( rand[,1], stabga[sel.ga.stab,4] )$p.val
text(2.5, 2, paste("h/sp p =", sprintf("%11.3e",pval_hyd_stabp) , sep=' '), cex=0.6)
text(2.5, 1, paste("h/sn p =", sprintf("%11.3e",pval_hyd_stabn) , sep=' '), cex=0.6)


dev.off()
