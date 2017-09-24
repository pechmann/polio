postscript("figures/Figure3/figure_sites_agg.ps", width=2.8, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")
par(mfrow=c(1,1))
par(mar=c(6,6,3,1))

data <- array(data, dim=c(2200, 2) )
data[,1] <- p.adjust(data[,1], method="BH")

sig_positions <- function(d){
  d <- as.matrix(d)
  s <- d[,1] < 0.05
  
  pos <- sum(d[s,2] > 1, na.rm=T)
  neg <- sum(d[s,2] < 1, na.rm=T)
  
  out <- c(pos, neg)
  out
}

b <-barplot(cbind(sig_positions(data[1:880,]), sig_positions(data[1566:1747,]), sig_positions(data[1748:2200,]))+1, space=c(0.2, 0.1, 0.8, 0.1, 0.2, 0.1), beside=T, col=c("#29ABE2", 2), border=F, axes=F, ylim=c(0, 75)) # pseudo-count of 1 added here as there are no sites in P3D otherwise)
axis(1, at=colMeans(b), labels=c("P1", "3C", "3D"),  cex.axis=1, tick=F)
axis(2, at=c(0, 25, 50, 75), cex.axis=1.2)
mtext("Sites", side=2, line=2.5, cex=1.2)
mtext("Aggregation prone mutations", side=3, line=1, cex=0.75)
legend("topright", legend=c("Hsp90N", "hsp90i"), bty='n', xpd=T, pch=15, cex=1, col=c("#29ABE2", 2))

# get p-values
P1_pval  <- fisher.test( array(c(sig_positions(data[1:880,]), 880-sig_positions(data[1:880,])), dim=c(2,2)) )$p.value
P3C_pval <- fisher.test( array(c(sig_positions(data[1566:1747,]), (1747-1566) - sig_positions(data[1566:1747,])), dim=c(2,2)) )$p.value
P3D_pval <- fisher.test( array(c(sig_positions(data[1748:2200,]), (2200-1748) - sig_positions(data[1748:2200,])), dim=c(2,2)) )$p.value

text(colMeans(b)[2], 44, paste("P1 p =", sprintf("%11.3e",P1_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 38, paste("P3C p =", sprintf("%11.3e",P3C_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 32, paste("P3D p =", sprintf("%11.3e",P3D_pval) , sep=' '), cex=0.6)

dev.off()



