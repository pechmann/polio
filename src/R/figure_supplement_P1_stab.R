#load data
struct <- read.table("data/reference/polio.asa", header=T)
struct <- as.matrix(struct)
data <- array(stab, dim=c(2200, 2)) 
data <- cbind(data[1:2200, 1:2], struct[1:2200,2:4])


#function to extract counts
getNumbs <- function(data, sel, theta){
  t <- sum(sel, na.rm=T)
  p <- sum( data[sel,1] < theta & data[sel,2] > 1, na.rm=T)
  n <- sum( data[sel,1] < theta & data[sel,2] < 1, na.rm=T)
  out <- (c(t, p, n))
  out
}


mutationsIntSurf <- function(data, norm=F){
  
  theta <- 0.01
  data <- as.matrix(data)
  
  out <- matrix(NA, ncol=3, nrow=8)
  
  P1 <- data[1:880,]
  P3C <- data[1566:1747,]
  P3D <- data[1748:2200,]
  
  sel.P1.buried <- P1[,5] < 30
  if (norm){
      out[1,] <- getNumbs(P1, sel.P1.buried, theta) / sum(sel.P1.buried, na.rm=T)
  } else { out[1,] <- getNumbs(P1, sel.P1.buried, theta)}
  
  sel.P1.monoint <- !sel.P1.buried & P1[,5] - P1[,4] > 30
  sel.P1.capsidint <- !sel.P1.buried & P1[,4] - P1[,3] > 30
  sel.P1.surf <- !sel.P1.capsidint & P1[,3] > 50
  
  out[2,] <- getNumbs(P1, sel.P1.monoint, theta) 
 
  if (norm){ 
      out[3,] <- getNumbs(P1, sel.P1.capsidint, theta) / sum(sel.P1.capsidint, na.rm=T) 
  } else { out[3,] <- getNumbs(P1, sel.P1.capsidint, theta) }

  if (norm){
  out[4,] <- getNumbs(P1, sel.P1.surf, theta) / sum(sel.P1.surf, na.rm=T)
  } else { out[4,] <- getNumbs(P1, sel.P1.surf, theta) }
  
  sel.P3C.buried <- P3C[,5] < 30
  out[5,] <- getNumbs(P3C, sel.P3C.buried, theta) 
  
  sel.P3C.exposed <- P3C[,5] > 30
  out[6,] <- getNumbs(P3C, sel.P3C.exposed, theta) 
  
  sel.P3D.buried <- P3D[,5] < 30
  out[7,] <- getNumbs(P3D, sel.P3D.buried, theta)
  
  sel.P3D.exposed <- P3D[,5] > 30
  out[8,] <- getNumbs(P3D, sel.P3D.exposed, theta)
  
  
  out
}



postscript("figures/Supplement/figure_supplement_P1_stab.ps", width=2, height=2, horizontal = FALSE,  onefile = FALSE, paper = "special") #_normalize
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))  

mutations_surface <- mutationsIntSurf(data, F)
mutations_surface_norm <- mutationsIntSurf(data, T)

b <- barplot(t(mutations_surface[c(1,3,4),2:3]), beside=T, border=F, ylim=c(0, 40), space=c(0.5, 0.1, 0.5, 0.1, 0.5, 0.1), axes=F, xlab="", ylab="", col=c("#29ABE2", 2))
axis(1, at=c(colMeans(b)), labels=c("b", "i", "s"), cex.axis=0.75, las=1, tick=F, line=0)
axis(2, at=c(0, 20, 40), cex.axis=0.75)
mtext("Sites", side=2, line=2.5, cex=1)
#legend("topright", legend=c("b: buried", "i: monomer interface", "s: capsid surface"), xpd=T, bty='n', cex=0.7)

# get p-values
buried_pval  <- fisher.test( array( c(mutations_surface[1,2], mutations_surface[1,1]-mutations_surface[1,2], mutations_surface[1,3], mutations_surface[1,1]-mutations_surface[1,3]), dim=c(2,2)) )$p.value
interface_pval <- fisher.test( array( c(mutations_surface[3,2], mutations_surface[3,1]-mutations_surface[3,2], mutations_surface[3,3], mutations_surface[3,1]-mutations_surface[3,3]), dim=c(2,2)) )$p.value
surface_pval <- fisher.test( array( c(mutations_surface[4,2], mutations_surface[4,1]-mutations_surface[4,2], mutations_surface[4,3], mutations_surface[4,1]-mutations_surface[4,3]), dim=c(2,2)) )$p.value

text(colMeans(b)[2], 25, paste("buried p =", sprintf("%11.2e",buried_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 21, paste("int p =", sprintf("%11.2e",interface_pval) , sep=' '), cex=0.6)
text(colMeans(b)[2], 17, paste("surf p =", sprintf("%11.2e",surface_pval) , sep=' '), cex=0.6)


dev.off()
