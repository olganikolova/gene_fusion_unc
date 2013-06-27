# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# makes sc's per gene and per fusion 
# =======================================

# set population size
# popsize <- c()

files <- as.list(c("BLCA.scat.fus","BRCA.scat.fus", "CESC.scat.fus","HNSC.scat.fus",
  			"KICH.scat.fus","KIRC.scat.fus", "KIRP.scat.fus","LGG.scat.fus", 
				"LIHC.scat.fus","LUAD.scat.fus","LUSC.scat.fus","PAAD.scat.fus",
				"PRAD.scat.fus","SKCM.scat.fus","THCA.scat.fus"))

myPlot <- function(file, size){

	d <- read.table(file, sep="\t", header=TRUE, as.is=TRUE);
	dd <- d[with(d, order(-Frequency)),];
	x <- 1:nrow(dd);
	tmp <- sub("^([^.]*).*","\\1",file);
	names <- as.character(d[,1])
	pdf(paste(file, ".pdf", sep=""));
	plot(x,dd[,2], xlab="Fusions", ylab="Frequency in Patients", main=paste(tmp, " (n = ", size, ")", sep=""), 
	xlim=c(0, nrow(dd)), ylim=c(0, 1))
	#axis(side=1, at=x, labels=names, las=2)
	dev.off()

}

lapply(as.list(1:length(files)), function(x){
		myPlot(files[[x]], popsize[x])
})

files <- as.list(c("BLCA.scat.gene","BRCA.scat.gene", "CESC.scat.gene","HNSC.scat.gene",
				"KICH.scat.gene","KIRC.scat.gene", "KIRP.scat.gene","LGG.scat.gene", 
				"LIHC.scat.gene","LUAD.scat.gene","LUSC.scat.gene","PAAD.scat.gene",
				"PRAD.scat.gene","SKCM.scat.gene","THCA.scat.gene"))

myPlot <- function(file, size){

	d <- read.table(file, sep="\t", header=TRUE);
	dd <- d[with(d, order(-Frequency)),];
	x <- 1:nrow(dd);
	tmp <- sub("^([^.]*).*","\\1",file);
	names <- as.character(d[,1])
	pdf(paste(file, ".pdf", sep=""));
	plot(x,dd[,2], xlab="Genes", ylab="Frequency in Patients", main=paste(tmp, " (n = ", size, ")", sep=""),
	xlim=c(0, nrow(dd)), ylim=c(0, 1))
	#axis(side=1, at=x, labels=names, las=2)
	dev.off()

}

lapply(as.list(1:length(files)), function(x){
		myPlot(files[[x]], popsize[x])
})





text(1:3, par("usr")[3] - 0.25, srt = 45, adj = 1, labels = names, xpd = TRUE)

foo <- plot(x,dd[,2], xlab="Fusions", ylab="% Observations", main=file, xaxt="n" )
	
text(1:3, par("usr")[3] - 0.2, labels = names, srt = 45, adj = 1, xpd = TRUE, cex=.9)

mtext(1, text = "X Axis Label", line = 6)
