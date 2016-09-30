#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")

## ./plot_auPRC_repeatedSubsampling.R ${eval_dir}/all_auPRC_5_10_repeatedSubsampling_FDR${fdr}.pdf $fdr $all_auPRC_all_tools 
args <- commandArgs(TRUE)

pdfToPlot <- args[1]
fdr <- as.numeric(args[2])
tools=list()
numTools=0
for(i in 3:length(args)){
	numTools=numTools+1
	tools[[numTools]]=args[i]
}

ylimMax=0.8

if( !file.exists(pdfToPlot)){
		auPRCvalues=list()
		for(i in 1:numTools){
			list1=read.table(tools[[i]])
			auPRCvalues[[i]]=c(list1[seq(1,10),1])
		}
		#colors=c(brewer.pal(9,"Set1")[-6])
		colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]
		pdf(pdfToPlot,pointsize=3)
		boxplot(auPRCvalues[[1]],auPRCvalues[[2]],auPRCvalues[[3]],auPRCvalues[[4]],auPRCvalues[[5]],auPRCvalues[[6]],auPRCvalues[[7]],auPRCvalues[[8]],auPRCvalues[[9]],col=colors,names=c("deepSNV", "GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper", "Varscan2"), ylab=paste("Area under Precision Recall Curve with ",fdr,"% FDR",sep=""),ylim=c(0,ylimMax),cex.lab=2.5,cex.axis=2.5,las=2)
}
