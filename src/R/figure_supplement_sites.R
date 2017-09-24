postscript("figures/Supplement/figure_supplement_sites.ps", width=9, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")
par(mfrow=c(1,3))
par(mar=c(6,6,3,1))

stab <- array(stab, dim=c(2200, 2) )
stab[,1] <- p.adjust(stab[,1], method="BH")

agg <- array(agg, dim=c(2200, 2) )
agg[,1] <- p.adjust(agg[,1], method="BH")

hyd <- array(hyd, dim=c(2200, 2) )
hyd[,1] <- p.adjust(hyd[,1], method="BH")

sig_positions <- function(d){
  d <- as.matrix(d)
  s <- d[,1] < 0.05
  
  pos <- sum(d[s,2] > 1, na.rm=T)
  neg <- sum(d[s,2] < 1, na.rm=T)
  
  out <- c(pos, neg)
  out
}


b <-barplot(cbind(sig_positions(stab[1:880,]), sig_positions(stab[1566:1747,]), sig_positions(stab[1748:2200,]))+1, space=c(0.2, 0.1, 0.8, 0.1, 0.2, 0.1), beside=T, col=c("#29ABE2", 2), border=F, axes=F, ylim=c(0, 120))
axis(1, at=colMeans(b), labels=c("P1", "3C", "3D"),  cex.axis=1, tick=F)
axis(2, at=c(0, 40, 80, 120), cex.axis=1.2)
mtext("Sites", side=2, line=2.5, cex=1.2)
mtext("Destabilizing mutations", side=3, line=1, cex=0.75)
legend("topright", legend=c("Hsp90N", "hsp90i"), bty='n', xpd=T, pch=15, cex=1, col=c("#29ABE2", 2))

# get p-values
P1_pval  <- fisher.test( array(c(sig_positions(stab[1:880,]), 880-sig_positions(stab[1:880,])), dim=c(2,2)) )$p.value
P3C_pval <- fisher.test( array(c(sig_positions(stab[1566:1747,]), (1747-1566) - sig_positions(stab[1566:1747,])), dim=c(2,2)) )$p.value
P3D_pval <- fisher.test( array(c(sig_positions(stab[1748:2200,]), (2200-1748) - sig_positions(stab[1748:2200,])), dim=c(2,2)) )$p.value

text(colMeans(b)[2], 44, paste("P1 p =", sprintf("%11.3e",P1_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 38, paste("P3C p =", sprintf("%11.3e",P3C_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 32, paste("P3D p =", sprintf("%11.3e",P3D_pval) , sep=' '), cex=0.6)


b <-barplot(cbind(sig_positions(agg[1:880,]), sig_positions(agg[1566:1747,]), sig_positions(agg[1748:2200,]))+1, space=c(0.2, 0.1, 0.8, 0.1, 0.2, 0.1), beside=T, col=c("#29ABE2", 2), border=F, axes=F, ylim=c(0, 50))
axis(1, at=colMeans(b), labels=c("P1", "3C", "3D"),  cex.axis=1, tick=F)
axis(2, at=c(0, 25, 50), cex.axis=1.2)
mtext("Sites", side=2, line=2.5, cex=1.2)
mtext("Aggregation prone mutations", side=3, line=1, cex=0.75)
legend("topright", legend=c("Hsp90N", "hsp90i"), bty='n', xpd=T, pch=15, cex=1, col=c("#29ABE2", 2))

# get p-values
P1_pval  <- fisher.test( array(c(sig_positions(agg[1:880,]), 880-sig_positions(agg[1:880,])), dim=c(2,2)) )$p.value
P3C_pval <- fisher.test( array(c(sig_positions(agg[1566:1747,]), (1747-1566) - sig_positions(agg[1566:1747,])), dim=c(2,2)) )$p.value
P3D_pval <- fisher.test( array(c(sig_positions(agg[1748:2200,]), (2200-1748) - sig_positions(agg[1748:2200,])), dim=c(2,2)) )$p.value

text(colMeans(b)[2], 44, paste("P1 p =", sprintf("%11.3e",P1_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 38, paste("P3C p =", sprintf("%11.3e",P3C_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 32, paste("P3D p =", sprintf("%11.3e",P3D_pval) , sep=' '), cex=0.6)



b <-barplot(cbind(sig_positions(hyd[1:880,]), sig_positions(hyd[1566:1747,]), sig_positions(hyd[1748:2200,]))+1, space=c(0.2, 0.1, 0.8, 0.1, 0.2, 0.1), beside=T, col=c("#29ABE2", 2), border=F, axes=F, ylim=c(0, 100))
axis(1, at=colMeans(b), labels=c("P1", "3C", "3D"),  cex.axis=1, tick=F)
axis(2, at=c(0, 50, 100), cex.axis=1.2)
mtext("Sites", side=2, line=2.5, cex=1.2)
mtext("Hydrophobic mutations", side=3, line=1, cex=0.75)
legend("topright", legend=c("Hsp90N", "hsp90i"), bty='n', xpd=T, pch=15, cex=1, col=c("#29ABE2", 2))


# get p-values
P1_pval  <- fisher.test( array(c(sig_positions(hyd[1:880,]), 880-sig_positions(hyd[1:880,])), dim=c(2,2)) )$p.value
P3C_pval <- fisher.test( array(c(sig_positions(hyd[1566:1747,]), (1747-1566) - sig_positions(hyd[1566:1747,])), dim=c(2,2)) )$p.value
P3D_pval <- fisher.test( array(c(sig_positions(hyd[1748:2200,]), (2200-1748) - sig_positions(hyd[1748:2200,])), dim=c(2,2)) )$p.value

text(colMeans(b)[2], 44, paste("P1 p =", sprintf("%11.3e",P1_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 38, paste("P3C p =", sprintf("%11.3e",P3C_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 32, paste("P3D p =", sprintf("%11.3e",P3D_pval) , sep=' '), cex=0.6)


dev.off()
