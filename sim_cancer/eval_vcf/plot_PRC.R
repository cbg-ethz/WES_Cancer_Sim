#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")

args <- commandArgs(TRUE)
path2PR <- args[1]
myPRcurveFile <- args[2]
currTitle <- args[3]


print(myPRcurveFile)
pdfToPlot <- paste(dirname(myPRcurveFile),"/",currTitle,"_PRCurve.pdf",sep="")

if(!file.exists(pdfToPlot)){
	PRvalues <- read.table(myPRcurveFile)
	pdf(pdfToPlot)
	precision=as.vector(as.numeric(PRvalues[,1]))
	recall=as.vector(as.numeric(PRvalues[,2]))
	numSNVs=length(recall)
	if(basename(myPRcurveFile) == "TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.PRcurve"){
		plot_IDX=seq(1,numSNVs,1000)
	} else {
		plot_IDX=seq(1,numSNVs,300)
	}
	plot(recall[plot_IDX], precision[plot_IDX], pch=20,main=currTitle, xlab="Recall", ylab="Precision",ylim=c(0,1),xlim=c(0,1))
	dev.off()
}


