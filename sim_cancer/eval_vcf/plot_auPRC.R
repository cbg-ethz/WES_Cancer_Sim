
library("RColorBrewer")

args <- commandArgs(TRUE)

for (i in 1:length(args) ){
	print(args[i])
}

pdfToPlot <- args[1]

deepSNVauPRC_5 <- args[2]
GATK_HPauPRC_5 <- args[3]
GATK_UGauPRC_5 <- args[4]
JSM2auPRC_5 <- args[5]
muTectauPRC_5 <- args[6]
SAMtauPRC_5 <- args[7]
sinvict_5 <- args[8]
SomSniperauPRC_5 <- args[9]
VarScanauPRC_5 <- args[10]

deepSNVauPRC <- args[11]
GATK_HPauPRC <- args[12]
GATK_UGauPRC <- args[13]
JSM2auPRC <- args[14]
muTectauPRC <- args[15]
SAMtauPRC <- args[16]
sinvictauPRC <- args[17]
SomSniperauPRC <- args[18]
VarScanauPRC <- args[19]

cont <- as.logical(args[20]) 	# TRUE, if the auPRC should be plotted as a function of the normal contamination
				# FALSE, if the auPRC should be plotted as a function of the coverage



ylimMax=0.8

if( !file.exists(pdfToPlot)){

	if (!cont){
		print("This is for the cov auPRC plot:")
		list1=read.table(deepSNVauPRC)
		xvals=list1[seq(1,5),1] # these are the different coverage levels, which are the same for all tools
		#colors=c(brewer.pal(9,"Set1")[-6])
		colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]

		pdf(pdfToPlot,pointsize=3)
		vals=list1[seq(1,5),2] # this is the auPRC
		plot(cbind(xvals,vals),type="l",col=colors[1],xlab="Median Coverage", ylab="Area under Precision Recall Curve",ylim=c(0,ylimMax),lwd=2,cex.lab=2.5,cex.axis=2.5,lty=2,xaxt="n")
		print("deepSNV auPRC's for 10% FDR:")
		print(cbind(xvals,vals))

		print(JSM2auPRC)
		print(SAMtauPRC)
		print(GATK_UGauPRC)
		read.table(SAMtauPRC)
		cnt=2
		for (tools in c(GATK_HPauPRC, GATK_UGauPRC, JSM2auPRC,muTectauPRC, SAMtauPRC, sinvictauPRC,SomSniperauPRC, VarScanauPRC)){ 
			print(tools)
			list_tool=read.table(tools)
			vals=list_tool[seq(1,5),2]
			print(cbind(xvals, vals))
			lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=2)
			cnt=cnt+1
		}

		cnt=1
		for (tools in c(deepSNVauPRC_5, GATK_HPauPRC_5, GATK_UGauPRC_5, JSM2auPRC_5, muTectauPRC_5, SAMtauPRC_5, sinvict_5, SomSniperauPRC_5, VarScanauPRC_5)){
			print(tools)
			print(cnt)
			list_tool=read.table(tools)
			vals=list_tool[seq(1,5),2]
			print(cbind(xvals, vals))
			lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=1)
			cnt=cnt+1
		}
		
		abline(v=xvals,lty=3)	
		#legend("bottomright", legend=tools, col=colors, lty=c(1,1))
		axis(1,at=xvals,labels=c(25,53,106,159,212),cex.lab=2.5,cex.axis=2.5) ## median coverage
		dev.off()
		print(cbind(xvals,vals))
	} else {
		list1=read.table(deepSNVauPRC)
		xvals=list1[seq(6,9),1] # these are the different contamination levels, which are the same for all tools
		#colors=c(brewer.pal(9,"Set1")[-6])
		colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]

		pdf(pdfToPlot,pointsize=3)
		vals=list1[seq(6,9),2] # this is the auPRC
		plot(cbind(xvals,vals),type="l",col=colors[1],xlab="Level of normal contamination (%)", ylab="Area under Precision Recall Curve",ylim=c(0,ylimMax),lwd=2,cex.lab=2.5,cex.axis=2.5,lty=2,xaxt="n")
		#print(JSM2auPRC)
		#print(SAMtauPRC)
		#print(GATK_UGauPRC)
		#read.table(SAMtauPRC)
		cnt=2
		for (tools in c(GATK_HPauPRC, GATK_UGauPRC, JSM2auPRC,muTectauPRC, SAMtauPRC, sinvictauPRC,SomSniperauPRC, VarScanauPRC)){
			print(tools)
			list_tool=read.table(tools)
			vals=list_tool[seq(6,9),2]
			print(cbind(xvals, vals))
			lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=2)
			cnt=cnt+1	
		}
		cnt=1
		for (tools in c(deepSNVauPRC_5, GATK_HPauPRC_5, GATK_UGauPRC_5, JSM2auPRC_5, muTectauPRC_5, SAMtauPRC_5, sinvict_5, SomSniperauPRC_5, VarScanauPRC_5)){
			print(tools)
			print(cnt)
			list_tool=read.table(tools)
			vals=list_tool[seq(6,9),2]
			print(cbind(xvals, vals))
			lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=1)
			cnt=cnt+1
		}
		abline(v=xvals,lty=3)	
		#legend("bottomright", legend=tools, col=colors, lty=c(1,1))
		axis(1,at=xvals,labels=c(10,20,40,60),cex.lab=2.5,cex.axis=2.5)
		dev.off()
		#print(cbind(xvals,vals))
	}
}

legendPlot=paste(sub(".pdf","",pdfToPlot),"_legend.pdf",sep="")
if( !file.exists(legendPlot)){
		print("This is for the legend:")
		list1=read.table(deepSNVauPRC)
		xvals=list1[seq(1,5),1] # these are the different coverage levels, which are the same for all tools
		#colors=c(brewer.pal(9,"Set1")[-6])
		#colors=colors[c(1,5,4,2,6,3,7,8)]
		colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]

		pdf(legendPlot,pointsize=3)
		vals=list1[seq(1,5),2] # this is the auPRC
		plot(cbind(xvals,vals),type="n",col=colors[1],xlab="Median Coverage", ylab="Area under Precision Recall Curve",ylim=c(0,ylimMax),lwd=2,cex.lab=2.5,cex.axis=2.5)

		#print(JSM2auPRC)
		#print(SAMtauPRC)
		#print(GATK_UGauPRC)
		#read.table(SAMtauPRC)
		cnt=2
		for (tools in c(GATK_HPauPRC, GATK_UGauPRC, JSM2auPRC, muTectauPRC, SAMtauPRC, sinvictauPRC,SomSniperauPRC, VarScanauPRC)){
			list_tool=read.table(tools)
			vals=list_tool[seq(1,5),2]
			#lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=2)
			cnt=cnt+1	
		}
		cnt=1
		for (tools in c(deepSNVauPRC_5, GATK_HPauPRC_5, GATK_UGauPRC_5, JSM2auPRC_5, muTectauPRC_5, SAMtauPRC_5, sinvict_5, SomSniperauPRC_5, VarScanauPRC_5)){
			list_tool=read.table(tools)
			vals=list_tool[seq(1,5),2]
			#lines(cbind(xvals, vals), col=colors[cnt],lwd=2,lty=1)
			cnt=cnt+1
		}
		#abline(v=xvals,lty=3)	
		legend("bottomright", legend=c("deepSNV","GATK HP","GATK UG","JointSNVMix2","MuTect","SAMtools","SiNVICT","somaticSniper","VarScan2"), col=colors, lty=1)
		legend("bottomleft",legend=c("10%", "5%"),lwd=2,title="False discovery cutoff",lty=c(2,1))
		
		dev.off()
		#print(cbind(xvals,vals))
}
