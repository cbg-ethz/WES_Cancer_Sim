
library("RColorBrewer")



##./plot_coverage_profiles.R $Rsummary $meanCov $sampleName
## $Rsummary $meanCov $sampleName $medianCov

args <- commandArgs(TRUE)
Rsummary <- args[1]
print(Rsummary)
dataToPlot <- read.table(Rsummary)
meanCov <- args[2]
sampleName <- args[3]

if (length(args)>3){
	medianCov <- args[4]
	percentCoverage=c(as.vector(as.numeric(dataToPlot[,1])))
	pdfToPlot <- paste(dirname(Rsummary),"/Coverage_Profile_",sampleName,".pdf",sep="")
		pdf(pdfToPlot)
		plot(percentCoverage,type="h",xlab="Coverage",ylab="Percentage covered with at least X or more reads",main=paste("Sample ",sampleName," has mean coverage ",meanCov,"X and median coverage ",medianCov,"X",sep=""),xaxt="n",ylim=c(0,100))
		axis(1,at=c(seq(1,500,50)),labels=paste(c(seq(1,500,50)),"x",sep=""),cex.axis=1)
		dev.off()
} else {

	percentCoverage=c(as.vector(as.numeric(dataToPlot[,2])))
	print(percentCoverage)

	pdfToPlot <- paste(dirname(Rsummary),"/Coverage_Profile_",sampleName,".pdf",sep="")
	#if(!file.exists(pdfToPlot)){
		pdf(pdfToPlot)
		plot(percentCoverage,type="h",xlab="Coverage",ylab="Percentage covered with at least X or more reads",main=paste("Sample ",sampleName," has mean coverage ",meanCov,sep=""),xaxt="n",ylim=c(0,100))
		axis(1,at=c(seq(1,51)),labels=ifelse(c(seq(1,51))%%5==0, paste(c(seq(1,51)),"x",sep=""),""),cex.axis=0.3)
		dev.off()
	#}
}

