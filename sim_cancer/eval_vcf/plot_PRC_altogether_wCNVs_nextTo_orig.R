#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")

args <- commandArgs(TRUE)
tools_PRC=list()
for(i in 1:18){
	tools_PRC[[i]] <- args[i]
}

print(tools_PRC)
pdfToPlot <- paste(dirname(tools_PRC[[1]]),"/Precision_Recall_Curves_all_wCNVs_and_orig.pdf",sep="")

myTitles <- c('deepSNV', 'GATK HP', 'GATK UG', 'JointSNVMix2', 'MuTect', 'SAMtools', 'SiNVICT', 'SomaticSniper', 'VarScan2')


if(!file.exists(pdfToPlot)){
	mySize=20
	pdf(pdfToPlot,width=mySize,height=mySize, pointsize=8)
	par(mfrow=c(3,3),mai=c(1,1,1,1))
	cnt=1
	for(i in c(1,3,5,7,9,11,13,15,17)){
		PRvalues <- read.table(tools_PRC[[i]])
		title <- sub("TU.wCont20.final.RG.50perc","",sub("_indel0.PRcurve","",basename(tools_PRC[[i]])))
		precision=as.vector(as.numeric(PRvalues[,1]))
		recall=as.vector(as.numeric(PRvalues[,2]))
		numSNVs=length(recall)
		if(basename(tools_PRC[[i]]) == "TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.PRcurve"){
			plot_IDX=seq(1,numSNVs,5000)
		} else {
			plot_IDX=seq(1,numSNVs,2000)
		}
		#plot(recall[plot_IDX], precision[plot_IDX], pch=20,main=title, xlab="Recall", ylab="Precision",ylim=c(0,1),xlim=c(0,1),cex.axis=3,cex.lab=3,cex.main=3)
		plot(recall[plot_IDX], precision[plot_IDX], pch=20,main=myTitles[cnt], xlab="Recall", ylab="Precision",ylim=c(0,1),xlim=c(0,1),cex.axis=5,cex.lab=5,cex.main=5, col="red")

		PRvalues=read.table(tools_PRC[[i+1]])
		precision=as.vector(as.numeric(PRvalues[,1]))
                recall=as.vector(as.numeric(PRvalues[,2]))
		numSNVs=length(recall)
		if(basename(tools_PRC[[i]]) == "TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.PRcurve"){
                        plot_IDX=seq(1,numSNVs,5000)
                } else {
                        plot_IDX=seq(1,numSNVs,2000)
                }
		points(recall[plot_IDX], precision[plot_IDX], col="black", pch=20)
		cnt=cnt+1
	}
	dev.off()
}


