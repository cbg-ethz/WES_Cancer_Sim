#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")

args <- commandArgs(trailingOnly = F)
script.dir <- dirname(sub("--file=","",args[grep("--file",args)]))
source(paste(script.dir,"/../../paths.R",sep=""))

for (i in 1:length(args)){
print(args[i])
}


ToSubstractFromFreq=0.025
doALL=TRUE

if(doALL){
	## Figure 1A in the paper

	all_eval_dirs=eval_dir_k20_wCNVs 
	
	for (curr_eval_dir in all_eval_dirs){
# 		names=c(paste(curr_eval_dir,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN",sep=""),
# 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN",sep=""))
# 
		names=c(paste(curr_eval_dir,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""))

 
		tools=c("deepSNV", "GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper", "Varscan2")
		num = length(names)
		colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]
		x=read.table(names[1]);
		len=dim(x)[1]/5
		xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
		xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq) 	
		# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one						
		# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants 
		# whose frequency is smaller than the next highest bin 
		# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
		#myLTY=c(2,3,1)
		#myLWD=c(2,3,2)
		myLTY=c(2,1) # without 5% FDR
		myLWD=c(2,2) # without 5% FDR
		numGroundTruthSNVs=x[(len*3+1):(len*4), 3] 
		##idxToPlot=which(numGroundTruthSNVs>=6000)
		idxToPlot=which(numGroundTruthSNVs>=3000)
	
	
		SN_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5.pdf",sep="")
		if ( ! file.exists(SN_pdf) ){
			pdf(SN_pdf,pointsize=3)
			#par(xpd=T, mar=par()$mar+c(6,6,0,6))
			cnt=1

			# 10% FDR
			block=3
			x=read.table(names[1])
			vals=x[(len*(block-1)+1):(len*block), 2]
			plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[cnt],ylim=c(0,1),lwd=myLWD[cnt],cex.lab=2.5,cex.axis=2.5)
			#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])),col=colors[1])
			abline(v=xvals[idxToPlot],lty=3,col="gray") #TODO: maybe change back

			print(cnt)
			print(block)
			for (i in 2:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the
										# sensitivity which is for a certain precision cutoff
				print(tail(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot]))))
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
				#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
	                }
	
			cnt=2
			#block=5 # 1% FDR
			block=4 # 5% FDR
			print(cnt)
			print(block)
			for (i in 1:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
										# sensitivity which is for a certain precision cutoff
				print(tail(cbind(xvals[idxToPlot], vals[idxToPlot])))
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
				#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
			}
	
	# 		block=4
	# 		cnt=3
	# 		vals = x[(len*(block-1)+1):(len*block), 2] 
	# 		for (i in 1:num)
	# 		{
	# 			x=read.table(names[i]); 
	# 			vals = x[(len*(block-1)+1):(len*block), 2]
	# 			lines(cbind(xvals, vals), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
	# 		}
				
			legend("bottomright",inset=c(-1,-1), legend=tools, col=colors, lwd=3)
			
	#		legend("bottomleft",inset=c(-1,-1),legend=c("10%", "5%", "1%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY)
			legend("bottomleft",inset=c(-1,-1),legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY)
			#par(mar=c(5, 4, 4, 2) + 0.1)
			dev.off()
		}
		legend_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_legend.pdf",sep="")
		if ( ! file.exists(legend_pdf) ){
			pdf(legend_pdf,pointsize=3)
			#par(xpd=T, mar=par()$mar+c(6,6,0,6))
			cnt=1
			i=1
			block=3
			x=read.table(names[1])
			vals=x[(len*(block-1)+1):(len*block), 2]
	
			plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="n", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[cnt],ylim=c(0,1),lwd=myLWD[cnt],cex.lab=2.5,cex.axis=2.5)
	
			for (i in 2:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the
										# sensitivity which is for a certain precision cutoff
				#lines(cbind(xvals[idxToPlot], vals[idxToPlot]), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
	                }
	
			cnt=2
			block=4
			for (i in 1:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
										# sensitivity which is for a certain precision cutoff
				#lines(cbind(xvals[idxToPlot], vals[idxToPlot]), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
			}
	
	
	# 		block=4
	# 		cnt=3
	# 		vals = x[(len*(block-1)+1):(len*block), 2] 
	# 		for (i in 1:num)
	# 		{
	# 			x=read.table(names[i]); 
	# 			vals = x[(len*(block-1)+1):(len*block), 2]
	# 			#lines(cbind(xvals, vals), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
	# 		}
			
			legend("bottomright", legend=tools, col=colors, lwd=3,cex=2.5)
			
			#legend("bottomleft",legend=c("10%", "5%", "1%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY, cex=2.5)
			legend("bottomleft",legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY,cex=2.5)
			#par(mar=c(5, 4, 4, 2) + 0.1)
			dev.off()
		}
		SN_inlet_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_inlet.pdf",sep="")
		if ( ! file.exists(SN_inlet_pdf) ){
			pdf(SN_inlet_pdf,pointsize=3)
			myLWD=c(rep(6,10))
			#par(xpd=T, mar=par()$mar+c(6,6,0,6))
			cnt=1
			i=1
			plot_untdxToPlotl=10	
			numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
			#idxToPlot=which(numGroundTruthSNVs>=200 && xvals<=0.08)
			idxToPlot=which(xvals<=0.08)

			block=3
			x=read.table(names[1])
			vals=x[(len*(block-1)+1):(len*block), 2]

			#plot(cbind(c(0,xvals[c(seq(1,plot_until))]), c(0,vals[c(seq(1,plot_until))])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[cnt],ylim=c(0,0.55),lwd=myLWD[cnt],cex.lab=2.5,cex.axis=2.5,xlim=c(0.004000,0.046))
			print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
			plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[cnt],ylim=c(0,0.55),lwd=myLWD[cnt],cex.lab=2.5,cex.axis=2.5,xlim=c(0.00,0.05))
			#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])),col=colors[1])
			abline(v=xvals[idxToPlot],lty=3,lwd=myLWD[1],col="gray") 
	
			for (i in 2:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the
										# sensitivity which is for a certain precision cutoff
				#lines(cbind(xvals[c(seq(1,plot_until))], vals[c(seq(1,plot_until))]), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
				print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
				#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
	                }
	
			cnt=2
			block=4
			for (i in 1:num)
			{
				x=read.table(names[i]);
				vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
										# sensitivity which is for a certain precision cutoff
				print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=myLTY[cnt],lwd=myLWD[cnt])
				#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
			}
			dev.off()
		}
	}
	print("Done with Figure 1A")
#quit()
} 




