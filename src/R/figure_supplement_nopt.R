postscript("figures/Supplement/figure_supplement_nopt.ps", width=8, height=3.5, horizontal = FALSE,  onefile = FALSE, paper = "special")
par(mfrow=c(1,2))
par(mar=c(6,6,4,2))


format_data <- function(data){

	data <- array(data, dim=c(2200, 2) )
	data <- data[1:880,]                                    # only P1
	data[,1] <- p.adjust(data[,1], method="BH")             # correct for multiple testing
	s <- data[,1] < 0.05
	data[!s,2] <- 1
	data <- log2(data[,2])

 	}

r1hela <- format_data(r1hela)
r1human <- format_data(r1human)
r2hela <-format_data(r2hela)
r2human <- format_data(r2human)


plot(-10, -10, xlim=c(0, 880), ylim=c(-2.75, 3.5), axes=F, xlab="", ylab="")
polygon(x=c(370, 410, 410, 370), y=c(-2.75, -2.75, 2, 2), col="lightgreen", border=F)
polygon(x=c(580, 620, 620, 580), y=c(-2.75, -2.75, 2, 2), col="lightgreen", border=F)

lines(r1hela, lwd=2)
lines(r1human, col=2, lwd=2, lty=2)

axis(2, at=c(-2, 0, 2), cex.axis=1.1)
axis(1, cex.axis=1.1)

mtext("Mutations to nonoptimal codons", side=2, line=3.2, cex=1.2)
mtext("Enrichment Hsp90N vs. hsp90i: log2 OR", side=2, line=2.3, cex=0.9)
mtext("Codon sequence", side=1, line=2.5, cex=1.2)
mtext("Amino acid sequence", side=3, line=2, cex=0.9)

text(72, 3.1, "Chain 4", cex=0.6)
text(147, 3.1, "Chain 2", cex=0.6)
text(415, 3.1, "Chain 3", cex=0.6)
text(655, 3.1, "Chain 1", cex=0.6)

c4 <- c(1, 69)
c2 <- c(70, 341)
c3 <- c(342, 579)
c1 <- c(580, 840)

plotDom <- function(d,cc){
  d <- d+35
  polygon(x=c(d[1],d[2],d[2],d[1]), y=c(3.25,3.25,3.45,3.45), col=cc, border=F)
}

plotDom(c1,"#333333")
plotDom(c2,"red")
plotDom(c3,"blue")
plotDom(c4,"gold")

axis(3, at=c(40, 240, 440, 640, 840), labels=c(0, 200, 400, 600, 800), cex.axis=0.8)

legend(700, -2,legend=c("tAI HeLa", "tAI human"), bty='n', xpd=T, lwd=2, col=c(1,2), lty=c(1,2))

plot(-10, -10, xlim=c(0, 880), ylim=c(-2.75, 3.5), axes=F, xlab="", ylab="")
polygon(x=c(360, 380, 380, 360), y=c(-2.75, -2.75, 2, 2), col="lightgreen", border=F)
polygon(x=c(560, 570, 570, 560), y=c(-2.75, -2.75, 2, 2), col="lightgreen", border=F)

lines(r2hela, lwd=2)
lines(r2human, col=2, lwd=2, lty=2)

axis(2, at=c(-2, 0, 2), cex.axis=1.1)
axis(1, cex.axis=1.1)

mtext("Mutations to nonoptimal codons", side=2, line=3.2, cex=1.2)
mtext("Enrichment Hsp90N vs. hsp90i: log2 OR", side=2, line=2.3, cex=0.9)
mtext("Codon sequence", side=1, line=2.5, cex=1.2)
mtext("Amino acid sequence", side=3, line=2, cex=0.9)

text(72, 3.1, "Chain 4", cex=0.6)
text(147, 3.1, "Chain 2", cex=0.6)
text(415, 3.1, "Chain 3", cex=0.6)
text(655, 3.1, "Chain 1", cex=0.6)

c4 <- c(1, 69)
c2 <- c(70, 341)
c3 <- c(342, 579)
c1 <- c(580, 840)

plotDom <- function(d,cc){
  d <- d+35
  polygon(x=c(d[1],d[2],d[2],d[1]), y=c(3.25,3.25,3.45,3.45), col=cc, border=F)
}

plotDom(c1,"#333333")
plotDom(c2,"red")
plotDom(c3,"blue")
plotDom(c4,"gold")

axis(3, at=c(40, 240, 440, 640, 840), labels=c(0, 200, 400, 600, 800), cex.axis=0.8)

legend(700, -2,legend=c("tAI HeLa", "tAI human"), bty='n', xpd=T, lwd=2, col=c(1,2), lty=c(1,2))

dev.off()


