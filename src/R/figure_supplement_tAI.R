postscript("figures/Supplement/figure_supplement_tAI.ps", width=7, height=4, horizontal = FALSE,  onefile = FALSE, paper = "special")

par(mfrow=c(1,2))
par(mar=c(6,6,2,1))

r1_hela <- array(r1hela, dim=c(2200, 8, 7) )
r1_human <- array(r1human, dim=c(2200, 8, 7) )
r2_hela <- array(r2hela, dim=c(2200, 8, 7) )
r2_human <- array(r2human, dim=c(2200, 8, 7) )

# filter out sites with less than 10 mutations in each passage 
s_r1_he <- rowSums(r1_hela[,6,] > 10) > 6 & rowSums(r1_hela[,8,] > 10) > 6
s_r1_hu <- rowSums(r1_human[,6,] > 10) > 6 & rowSums(r1_human[,8,] > 10) > 6
s_r2_he <- rowSums(r1_hela[,6,] > 10) > 6 & rowSums(r2_hela[,8,] > 10) > 6
s_r2_hu <- rowSums(r1_human[,6,] > 10) > 6 & rowSums(r2_human[,8,] > 10) > 6

boxplot( r1_hela[s_r1_he,2,1]-r1_hela[s_r1_he,1,1], rowMeans(r1_hela[s_r1_he,3,])-r1_hela[s_r1_he,1,1], rowMeans(r1_hela[s_r1_he,4,])-r1_hela[s_r1_he,1,1], r1_human[s_r1_hu,2,1]-r1_human[s_r1_hu,1,1], rowMeans(r1_human[s_r1_hu,3,])-r1_human[s_r1_hu,1,1], rowMeans(r1_human[s_r1_hu,4,])-r1_human[s_r1_hu,1,1], axes=F, at=c(0.2, 0.5, 0.8, 1.2, 1.5, 1.8), boxwex=0.25, col=c("#777777", "#29ABE2", 2), main="Replica 1")
axis(1, at=c(0.2, 0.5, 0.8, 1.2, 1.5, 1.8), labels=c("Random", "Hsp90N", "hsp90i", "Random", "Hsp90N", "hsp90i"), cex.axis=0.7, las=2, tick=F, line=0)
axis(1, at=c(0.5, 1.5), labels=c("HeLa tAI", "Human tAI"), cex.axis=0.7, line=1.5, tick=F)
axis(2, cex.axis=1.1)
mtext("Relative tAI of mutations", side=2, line=2.5)

boxplot(r2_hela[s_r2_he,2,1]-r2_hela[s_r2_he,1,1], rowMeans(r2_hela[s_r2_he,3,])-r2_hela[s_r2_he,1,1], rowMeans(r2_hela[s_r2_he,4,])-r2_hela[s_r2_he,1,1], r2_human[s_r2_hu,2,1]-r2_human[s_r2_hu,1,1], rowMeans(r2_human[s_r2_hu,3,])-r2_human[s_r2_hu,1,1], rowMeans(r2_human[s_r2_hu,4,])-r2_human[s_r2_hu,1,1], axes=F, at=c(0.2, 0.5, 0.8, 1.2, 1.5, 1.8), boxwex=0.25, col=c("#777777", "#29ABE2", 2), main="Replica 2")
axis(1, at=c(0.2, 0.5, 0.8, 1.2, 1.5, 1.8), labels=c("Random", "Hsp90N", "hsp90i", "Random", "Hsp90N", "hsp90i"), cex.axis=0.7, las=2, tick=F, line=0)
axis(1, at=c(0.5, 1.5), labels=c("HeLa tAI", "Human tAI"), cex.axis=0.7, line=1.5, tick=F)
axis(2, cex.axis=1.1)
mtext("Relative tAI of mutations", side=2, line=2.5)

dev.off()

