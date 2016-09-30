#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")

args <- commandArgs(trailingOnly = F)
script.dir <- dirname(sub("--file=","",args[grep("--file",args)]))
source(paste(script.dir,"/../../paths.R",sep=""))


ToSubstractFromFreq=0.025
ALLBLOCKS=FALSE
FIGURE1A=FALSE
COMPARETOSNK1=FALSE
COMBISOVERLAP=FALSE
ALLTUNES=FALSE
TUNING=FALSE
LOCREBINOM=FALSE
FIGURE1ACOVERAGE=TRUE

## Compare vsn_k20 and sn_k1 directly
if(COMPARETOSNK1){
	names=c(paste(eval_dir_k20,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""),
		
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_snk1,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""))

		toolsRaw=c("deepSNV","GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper","VarScan2")
		tools=c(paste(toolsRaw," MS",sep=""),paste(toolsRaw," D",sep=""))
		num = length(names)
		colors=rainbow(num)
	
		x=read.table(names[1]); 
		len=dim(x)[1]/5
		xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
		xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq) 	
		# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one						
		# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants 
		# whose frequency is smaller than the next highest bin 
		# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
		#print(xvals)
		numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
		idxToPlot=which(numGroundTruthSNVs>=200)
		block=3 # 10% FDR

		pdfToPlot=paste(eval_dir_k20,"Compare_vsn_k20_sn_k1_SN",block,".pdf",sep="")
		if ( ! file.exists(pdfToPlot) ){
			pdf(pdfToPlot)
			vals = x[(len*(block-1)+1):(len*block), 2] 	# this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the sensitivity
									# which is for a certain precision cutoff
			plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity")
			#print(names[1])
			#print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
		
			for (i in 2:num)
			{
				x=read.table(names[i]); 
				vals = x[(len*(block-1)+1):(len*block), 2]
				if(i==14 || i==6){
					print(names[i])
					print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
					#quit()	
					#break
				}
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
			}
				
			legend("bottomright", legend=tools[c(seq(1,i))], col=colors[c(seq(1,i))], lty=c(1,1))
		
			dev.off()
		}
#quit()
}


if(COMBISOVERLAP){
	## Figure 3D - combinations and overlaps - or rather 3C now
	eval_overlap_dir_specific=paste(overlap_dir,"/specific_overlaps/overlaps/eval_10000_160527/",sep="")
	eval_combinations=paste(algmt_k20_dir,"/variants_combined/eval_10000_160527/",sep="")
	eval_combinations_sinvict=paste(algmt_k20_dir,"/variants_combined_sinvict/eval_10000_160527/",sep="")

	names=c(paste(eval_overlap_dir_specific,"overlap_deepSNV_joint_varscan_muTect_sinvict_original_joint.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_overlap_dir_specific,"overlap_deepSNV_joint_varscan_muTect_samvar_1_2_original_joint.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_overlap_dir_specific,"overlap_deepSNV_joint_varscan_muTect_original_joint.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_overlap_dir_specific,"overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_gatk_SNVs_gatkHPCaller_samvar_1_2_original_joint.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_combinations_sinvict,"combinedlistbetter_9_all_9_prod.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_combinations_sinvict,"combinedlistbetter_5_deep_varsc_mutec_joint_sinv_prod.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_combinations,"combinedlistbetter_4_deepSNV_varscan2_muTect_jointSNVMix2_prod.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_combinations,"combinedlistbetter_5_deepSNV_varscan2_samvar_1_2_muTect_jointSNVMix2_prod.vcf_indel0.eval_SN_summaryFreq",sep=""))
		
	print(names)
	
	tools=c("ov deep,Joint,MuT,VarSc sinvi",
		"ov deep,Joint,MuT,SAMtools,VarSc", 
		"ov deep,Joint,MuT,VarSc",
		"ov all nine tools", 
		"comb all nine",
		"comb deep,Joint,MuT,VarSc sinvi",
		"comb deep,Joint,MuT,VarSc",
		"comb deep,Joint,MuT,SAM,VarSc")

	print(tools)
	
	stopifnot(length(names)==length(tools))
	num = length(names)
	#colors=c(rep(c(brewer.pal(8,"Accent")[6], brewer.pal(11,"PiYG")[8], brewer.pal(9,"Oranges")[4], brewer.pal(9,"GnBu")[6]),2)) 
	colors=c(rep(c(brewer.pal(8,"Accent")[6], brewer.pal(11,"PiYG")[8], brewer.pal(9,"Oranges")[4], "black"),2)) # glorious 5, glorious 4, 4+sam, all nine
	colors=colors[c(1,3,2,4,4,1,2,3)]
	stopifnot(length(colors)==length(tools))

	x=read.table(names[1]); 
	len=dim(x)[1]/5
	xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
	xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq) 	
	# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one						
	# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants 
	# whose frequency is smaller than the next highest bin 
	# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
	#print(xvals)
	numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
	idxToPlot=which(numGroundTruthSNVs>=200)
	block=3 # 10% FDR

	myLTY=c(2,2,2,2,3,3,3,3)
	myLWD=ifelse(myLTY==3,3,2)
	
	comb_overlap_pdf=paste(eval_dir_k20,"combine_overlap_all_block",block,".pdf",sep="")
	if ( ! file.exists(comb_overlap_pdf) ){
		pdf(comb_overlap_pdf,pointsize=3)
		##pdf(comb_overlap_pdf)
		vals = x[(len*(block-1)+1):(len*block), 2] 	# this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the sensitivity
									# which is for a certain precision cutoff
		plot(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),lty=myLTY[1],lwd=myLWD[1],cex.lab=2.5,cex.axis=2.5)

		abline(v=xvals[idxToPlot],lty=3,col="gray")

		print(names[1])
		print(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])))
		for (i in 2:num)
		{
			print(names[i])
			x=read.table(names[i]); 
			vals = x[(len*(block-1)+1):(len*block), 2]
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i],lty=myLTY[i],lwd=myLWD[i])
			print(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])))
		}
			
		#legend("bottomright", legend=sub("_original","",sub("pos.sorted","",sub(".vcf","",sub("bam_minvarfreq","",sub("SNVs_Raw.vcf","",sub("vcf_indel0.eval_SN_summaryFreq","",sub("txt.snp.Somatic_qual.vcf","",sub("listbetter","",basename(names))))))))), col=colors, lty=c(1,1))
		#legend("bottomright", legend=tools, col=colors, lty=myLTY, lwd=myLWD)
		dev.off()
	}

	SN_inlet_pdf=paste(eval_dir_k20,"combine_overlap_all_block",block,"_inlet.pdf",sep="")
	if ( ! file.exists(SN_inlet_pdf) ){
		pdf(SN_inlet_pdf,pointsize=3)
		#par(xpd=T, mar=par()$mar+c(6,6,0,6))
		cnt=1
		i=1
		myLWD=c(rep(6,num))
		myLWD=ifelse(myLTY==3,8,6)
		idxToPlot=which(xvals<=0.08)
		x=read.table(names[1])
		vals=x[(len*(block-1)+1):(len*block), 2]

		print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
		plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[i],ylim=c(0,0.55),lwd=myLWD[i],cex.lab=2.5,cex.axis=2.5,xlim=c(0.00,0.05))
		abline(v=xvals[idxToPlot],lty=3,lwd=myLWD[1],col="gray")

		for (i in 2:num)
		{
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2] # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the
									# sensitivity which is for a certain precision cutoff
			print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=myLTY[i],lwd=myLWD[i])
                }
		dev.off()
	}
	comb_overlap_pdf=paste(eval_dir_k20,"combine_overlap_all_block",block,"_legend.pdf",sep="")
	if ( ! file.exists(comb_overlap_pdf) ){
	#pdf(comb_overlap_pdf,pointsize=3)
	pdf(comb_overlap_pdf)
	myLTY=c(2,2,2,2,3,3,3,3)
	myLWD=ifelse(myLTY==3,3,2)
	vals = x[(len*(block-1)+1):(len*block), 2] 	# this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the sensitivity
								# which is for a certain precision cutoff
	plot(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), type="n", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),xlim=c(0,1),lty=myLTY[1],lwd=myLWD[1])
	print(names[1])
	print(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])))
	for (i in 2:num)
	{
		print(names[i])
		x=read.table(names[i]); 
		vals = x[(len*(block-1)+1):(len*block), 2]
		#lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i],lty=myLTY[i],lwd=myLWD[i])
		print(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])))
	}
		
	#legend("bottomright", legend=sub("_original","",sub("pos.sorted","",sub(".vcf","",sub("bam_minvarfreq","",sub("SNVs_Raw.vcf","",sub("vcf_indel0.eval_SN_summaryFreq","",sub("txt.snp.Somatic_qual.vcf","",sub("listbetter","",basename(names))))))))), col=colors, lty=c(1,1))
	legend("bottomright", legend=tools, col=colors, lty=myLTY, lwd=myLWD)
	dev.off()
	}
#quit()
}

if(ALLBLOCKS){
	## just one tool of choice, and then all blocks

	names=c(paste(eval_dir_k20,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq",sep=""),
		paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""))

	tools=c("deepSNV", "JointSNVMix2", "SAMtools", "GATK UG", "GATK HP", "MuTect", "SiNVICT", "SomaticSniper","Varscan2")
	myLTY=c(3,2,1,4,6,5)
	myLWD=c(3,2,2,2,2,2)
	
	
	for (tool in 1:9){
		curr_names=names[tool]
		num = length(curr_names)
		colors=rainbow(num)
	
		x=read.table(curr_names[1]);
		len=dim(x)[1]/5
		xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
		xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq) 	
		# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one						
		# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants 
		# whose frequency is smaller than the next highest bin 
		# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
		numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
		idxToPlot=which(numGroundTruthSNVs>=200)
	
		curr_pdf=paste(eval_dir_k20,"SN_allPrec_",sub(" ","",tools[tool]),".pdf",sep="")
		if (!file.exists(curr_pdf)){
			pdf(curr_pdf)
			block=1
			cnt=1
			vals = x[(len*(block-1)+1):(len*block), 2]
			plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[1], xlab="Variant Frequency", ylab="Sensitivity",lty=myLTY[cnt],lwd=myLWD[cnt],ylim=c(0,1))
			for (block in 2:5){
				cnt=cnt+1
				vals = x[(len*(block-1)+1):(len*block), 2] 
				lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[1],lty=myLTY[cnt],lwd=myLWD[cnt])
			}
			legend("bottomright", legend=c("30%","20%","10%","5%","1%"), lty=myLTY,lwd=myLWD)
			dev.off()
		}
	}
	
} 

if(FIGURE1A){
	## Figure 1A in the paper

	all_eval_dirs=eval_dir_k20

	for (curr_eval_dir in all_eval_dirs){

 		names=c(paste(curr_eval_dir,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq",sep=""),
 			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""))
# 

		tools=c("deepSNV", "GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper", "Varscan2")
		num = length(names)
		#colors=rainbow(num)
		#colors=brewer.pal(8,"Accent")
		#colors=c(brewer.pal(9,"Set1")[-6])
		#colors=colors[c(1,5,4,2,6,3,7,8)]
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
		numGroundTruthSNVs=x[(len*3+1):(len*4), 3] ## hast to be 3 if summary!!, 4 if eval_SN
		idxToPlot=which(numGroundTruthSNVs>=200)
	
	
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
			abline(v=xvals[idxToPlot],lty=3,col="gray")

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
	
			legend("bottomright",inset=c(-1,-1), legend=tools, col=colors, lwd=3)
			legend("bottomleft",inset=c(-1,-1),legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY)
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
	
			
			legend("bottomright", legend=tools, col=colors, lwd=3,cex=2.5)
			legend("bottomleft",legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY,cex=2.5)
			#par(mar=c(5, 4, 4, 2) + 0.1)
			dev.off()
		}
		SN_inlet_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_inlet.pdf",sep="")
		if ( ! file.exists(SN_inlet_pdf) ){
			pdf(SN_inlet_pdf,pointsize=3)
			#par(xpd=T, mar=par()$mar+c(6,6,0,6))
			myLWD=c(6,6,6,6,6,6,6,6,6,6)
			cnt=1
			i=1
			plot_untdxToPlotl=10	
			numGroundTruthSNVs=x[(len*3+1):(len*4), 3] ## has to be 3 for summary!
			print(idxToPlot)
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


if(ALLTUNES){
	print("plot all tunes")
	## plot all tunes 

	normalVersionEvalDir=c(paste(base_dir,"/alignments_vsn_k20/variants_default/eval_10000_160527/",sep=""))
	all_eval_dirs=eval_dir_tuned

	tools=c("SAMtools","deepSNV","JointSNVMix2","VarScan2","SiNVICT")
	
	for (curr_eval_dir in all_eval_dirs){
		for (tool in tools){
			if (tool=="SAMtools"){
				names=c(paste(curr_eval_dir,c(list.files(path=curr_eval_dir,pattern="TU.wCont20.final.RG.50perc.*\\_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$")),sep=""),paste(normalVersionEvalDir,"/TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq",sep=""))
			} else if (tool=="deepSNV"){
				names=c(paste(curr_eval_dir, c(list.files(path=curr_eval_dir,pattern="TU.wCont20.final.RG.50perc_.*deepSNV.*vcf_indel0.eval_SN_summaryFreq$")),sep=""),paste(normalVersionEvalDir,"/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq",sep=""))
			} else if (tool=="JointSNVMix2"){
				names=c(paste(curr_eval_dir,c(list.files(path=curr_eval_dir,pattern="TU.wCont20.final.RG.50perc_.*\\joint.*\\.vcf_indel0.eval_SN_summaryFreq$")),sep=""),paste(normalVersionEvalDir,"/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq",sep=""))
			} else if (tool=="VarScan2"){
				names=c(paste(curr_eval_dir,c(list.files(path=curr_eval_dir,pattern="TU.wCont20.final.RG.50perc.bam_.*\\varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$")),sep=""),paste(normalVersionEvalDir,"/TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq",sep=""))
			} else if (tool=="SiNVICT"){
				names=c(paste(curr_eval_dir,c(list.files(path=curr_eval_dir,pattern="sinvict.*\\vcf_indel0.eval_SN_summaryFreq$")),sep=""),paste(normalVersionEvalDir,"/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq",sep=""))
			}
	
			print(names)
			num_tunes=length(names)
			print(num_tunes)
			rawColors=c(brewer.pal(8,"Accent")[-c(4,7,8)],brewer.pal(8,"Purples")[7],brewer.pal(11,"RdYlBu")[8],brewer.pal(8,"YlOrRd")[7],brewer.pal(8,"BuGn")[8],brewer.pal(8,"Oranges")[6],brewer.pal(8,"Accent")[c(7,8)])
			colors=c(rawColors,rawColors)
			MyLty=c(rep(c(1,2,3,4,6),4))
			#myLWD=c(3,2,2)
	
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
			idxToPlot=which(numGroundTruthSNVs>=200)
	
			SN_pdf=paste(curr_eval_dir,"SN_tuned_",tool,".pdf",sep="")
			if ( ! file.exists(SN_pdf) ){
				pdf(SN_pdf,pointsize=3)
				block=3
				i=1
				vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
										# sensitivity which is for a certain precision cutoff
				print(vals)
				plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=2.5)
				#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])),col=colors[i])
				abline(v=xvals[idxToPlot],lty=3,col="gray")

	
				for (i in 2:num_tunes)
				{
					x=read.table(names[i]);
					vals = x[(len*(block-1)+1):(len*block), 2]
					lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=MyLty[i],lwd=2.5)
					#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
				}
	
				dev.off()
			}
	
			legend_pdf=paste(curr_eval_dir,"SN_tuned_",tool,"_legend.pdf",sep="")
			if ( ! file.exists(legend_pdf) ){
				pdf(legend_pdf,pointsize=3)
				block=3
				i=1
				vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
										# sensitivity which is for a certain precision cutoff
				plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="n", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],ylim=c(0,1),cex.lab=2.5,cex.axis=2.5)
				
				for (i in 2:num_tunes)
				{
					x=read.table(names[i]);
					vals = x[(len*(block-1)+1):(len*block), 2]
					#lines(cbind(xvals, vals), col=colors[i], lty=MyLty[i])
				}
				print("legend:")
				print(sub("_vs_NO_final.RG.50perc_somatic","",sub("jointSNVMix2_SNVs_Raw.","",sub("_indel0.eval_SN_summaryFreq","",sub("TU.wCont20.final.RG.50perc","",basename(names))))))
				legend("bottomright", legend=sub("_vs_NO_final.RG.50perc_somatic","",sub("_NO_final.RG.50perc_alternative_two.sided_","",sub("jointSNVMix2_SNVs_Raw.","",sub("_indel0.eval_SN_summaryFreq","",sub("TU.wCont20.final.RG.50perc","",basename(names)))))), col=colors, lwd=3,cex=2.5,lty=MyLty)
				dev.off()
			}
		}
	}
} # enf of if FALSE number 3


if(TUNING){	
	## Figure 3 in Paper - or rather 3B now, the tuning plot
	eval_dirs=c(paste(base_dir,"alignments_vsn_k20/variants_tuned/eval_10000_160527/",sep=""),paste(base_dir,"alignments_vsn_k20/variants_default/eval_10000_160527/",sep=""))


	names=c(paste(eval_dirs[1],c(   list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.option_B_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.option_C20_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.option_C30_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$")),sep=""),
		paste(eval_dirs[2],	list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),sep=""), 
		paste(eval_dirs[1],c(	list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc_min_base_qual_25_min_map_qual_0.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$")),sep=""),
		paste(eval_dirs[2],     list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),sep=""),
		paste(eval_dirs[1],c(	list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided_minBaseQ0.deepSNV.vcf_indel0.eval_SN_summaryFreq$")),sep=""),
		paste(eval_dirs[2],	list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq$"),sep=""),
		paste(eval_dirs[1],c(	list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.bam_strandfilter1_varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$")),sep=""),
		paste(eval_dirs[2],	list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$"),sep=""),
		paste(eval_dirs[1], c( list.files(path=eval_dirs[1],pattern="sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq$")),sep=""),
		paste(eval_dirs[2],	list.files(path=eval_dirs[2],pattern="sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq$"),sep="") )
	
	print(names)
	num_tunes=length(names)
	print(num_tunes)
	Rawcolors=c(brewer.pal(9,"Set1")[-6])
	colors=Rawcolors[c(rep(3,5),rep(2,2),rep(1,2),rep(8,3))]
	colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)][c(rep(6,5),rep(4,2),rep(1,2),rep(9,3),rep(7,2))]
	MyLty=c(2,3,4,5,1,2,1,2,1,2,3,1,2,1)
	myLWD=c(rep(2,length(MyLty)))
	myLWD=ifelse(MyLty==1, 2, 2.5)
	myLWD=ifelse(MyLty==3, 3, myLWD)
	myLWD=c(2.5,3,2.5,2.5,2,2.5,2,2.5,2,3,3,2,2.5,2)
	x=read.table(names[1])
	len=dim(x)[1]/5
	xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
	xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq)
		# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one
		# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants
		# whose frequency is smaller than the next highest bin
		# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
	numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
	idxToPlot=which(numGroundTruthSNVs>=200)
	
	tuned_SN_pdf=paste(eval_dirs[1],"SN_tuned_deepJointSAMvarsc.pdf",sep="")
	if ( ! file.exists(tuned_SN_pdf) ){
		pdf(tuned_SN_pdf,pointsize=3)
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		#print(vals)
		print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
		plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i])
		#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])),col=colors[i])
		abline(v=xvals[idxToPlot],lty=3,col="gray")
		for (i in 2:num_tunes)
		{
			print(names[i])
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			print(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])))
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
			#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
		}
		dev.off()
	}

	legend_pdf=paste(eval_dirs[1],"SN_tuned_deepJointSAMvarsc_legend.pdf",sep="")
	if ( ! file.exists(legend_pdf) ){
		pdf(legend_pdf,pointsize=3)
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		plot(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), type="n", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],ylim=c(0,1),cex.lab=2.5,cex.axis=2.5)
		for (i in 2:num_tunes)
		{
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			#lines(cbind(xvals, vals), col=colors[i], lty=MyLty[i])
		}
		legend("bottomright", legend=sub("_NO_final.RG.50perc_alternative_two.sided_","",sub("jointSNVMix2_SNVs_Raw.","",sub("_indel0.eval_SN_summaryFreq","",sub("TU.wCont20.final.RG.50perc","",basename(names))))), col=colors, lwd=myLWD,cex=2.5,lty=MyLty)
		dev.off()
	}

	## inlet for Figure 3B - showing a zoom into the low frequency area 
	
	print(names)
	num_tunes=length(names)
	print(num_tunes)
	## Color and style as above
	#Rawcolors=c(brewer.pal(9,"Set1")[-6])
	#colors=Rawcolors[c(rep(3,4),rep(2,2),rep(1,2),rep(7,3))]
	#MyLty=c(2,3,4,1,2,1,2,1,2,3,1)
	#myLWD=c(3,2,2)

	x=read.table(names[1]);
	len=dim(x)[1]/5
	xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
	xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq)
	#idxToPlot=which(numGroundTruthSNVs>=2000 && xvals<=0.08)
	idxToPlot=which(xvals<=0.08)
	#plot_until=10
	
	tuned_SN_pdf=paste(eval_dirs[1],"SN_tuned_deepJointSAMvarsc_inlet_lowFreqVariants.pdf",sep="")
	if ( ! file.exists(tuned_SN_pdf) ){
		pdf(tuned_SN_pdf,pointsize=3)
		myLWD=c(rep(6,20))
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		print(vals)
		#plot(cbind(xvals[c(seq(1,plot_until))], vals[c(seq(1,plot_until))]), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i], ylim=c(0,0.55),xlim=c(0.004000,0.046))
		plot(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i], ylim=c(0,0.55),xlim=c(0.00,0.05))
		#points(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), col=colors[i])
		abline(v=xvals[idxToPlot],lty=3,lwd=myLWD[1],col="gray")
	
		for (i in 2:num_tunes)
		{
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			#lines(cbind(xvals[c(seq(1,plot_until))], vals[c(seq(1,plot_until))]), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
			#points(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i])
		}
		dev.off()
	}
	
#quit()
print("Done with Figure 3B")
} 


## Figure 3A in Paper - locRe and binomTest
if (LOCREBINOM){
	eval_dirs=c(paste(base_dir,"alignments_vsn_k20/variants_default/eval_10000_binomTest_160527/",sep=""),paste(base_dir,"alignments_vsn_k20/variants_default/eval_10000_160527/",sep=""))

	names=c(paste(eval_dirs[1],c(   # binomTest lists:
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_cpp_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[1],pattern="TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$")),sep=""),

					# default lists:
		paste(eval_dirs[2],c(	list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$"),
					
					# local realign lists:
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign_NO_final.RG.50perc.locRealign_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="sinvict_TU.wCont20.final.RG.50perc.locRealign_vs_NO_final.RG.50perc.locRealign_somatic.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq$"),
					list.files(path=eval_dirs[2],pattern="TU.wCont20.final.RG.50perc.locRealign.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq$"))
					,sep=""))

	tools=c("deepSNV", "GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper", "VarScan2")
	print(names)
	num_tunes=length(names)
	print(num_tunes)
	#Rawcolors=c(brewer.pal(9,"Set1")[-6])
	Rawcolors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)]
	colors=c(Rawcolors,Rawcolors,Rawcolors)
	#colors=c(Rawcolors[c(1,5,4,2,6,3,7,8)],Rawcolors[c(1,5,4,2,6,3,7,8)],Rawcolors[c(1,5,4,2,6,3,7,8)])
	
	MyLty=c(rep(2,9),rep(1,9),rep(3,9))
	myLWD=c(rep(2,length(MyLty)))
	myLWD=ifelse(MyLty==1, 2, 2.5)
	myLWD=ifelse(MyLty==3, 3, myLWD)
	##myLWD=c(2.5,3,2.5,2.5,2,2.5,2,2.5,2,2.5,3,2)
	x=read.table(names[1])
	len=dim(x)[1]/5
	xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
	xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq)
		# right now, the xvals go from 0.05-1.05, but the bins mean all variants whose frequency is smaller than the specified one
		# for plotting, we subtract 0.025 to have the range 0.025, 0.05, ..., 1.025, and the sensitivity in each bin is for all variants
		# whose frequency is smaller than the next highest bin
		# the special case 1.025 is set to 1.0, because for this bin, we know it contains all variants whose frequency is exactly 1.0
	numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
	idxToPlot=which(numGroundTruthSNVs>=200)

	tuned_SN_pdf=paste(eval_dirs[2],"SN_locRe_binom_allTools.pdf",sep="")
	if ( ! file.exists(tuned_SN_pdf) ){
		pdf(tuned_SN_pdf,pointsize=3)
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		print(vals)
		#plot(cbind(xvals, vals), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i])
		plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",ylim=c(0,1),lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i])
		abline(v=xvals[idxToPlot],lty=3,col="gray")

		for (i in 2:num_tunes)
		{
			print(names[i])
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			#lines(cbind(xvals, vals), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
		}
		dev.off()
	}
	
	legend_pdf=paste(eval_dirs[2],"SN_locRe_binom_allTools_legend.pdf",sep="")
	if ( ! file.exists(legend_pdf) ){
		pdf(legend_pdf,pointsize=3)
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		#plot(cbind(xvals, vals), type="n", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],ylim=c(0,1),cex.lab=2.5,cex.axis=2.5)
		plot(cbind(c(0,xvals[idxToPlot]),c(0,vals[idxToPlot])), type="n", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],ylim=c(0,1),cex.lab=2.5,cex.axis=2.5)
		for (i in 2:num_tunes)
		{
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			#lines(cbind(xvals, vals), col=colors[i], lty=MyLty[i])
		}
				
		#legend("bottomright", legend=tools, col=colors, lwd=3,cex=2.5)
	
		legend("bottomright", legend=c("binomial test","standard","local realignment"),lwd=myLWD[c(1,10,19)],title="Pipeline modifications",lty=MyLty[c(1,10,19)],cex=2.5)
	
		#legend("bottomright", legend=sub("SOMATIC_paired","",sub("jointSNVMix2_SNVs_Raw.","",sub("_indel0.eval_SN","",sub("TU.wCont20.final.RG.50perc","",basename(names))))), col=colors, lwd=myLWD,cex=2.5,lty=MyLty)
		dev.off()
	}

	## inlet for Figure 3A - showing a zoom into the low frequency area 
	
	print(names)
	num_tunes=length(names)
	print(num_tunes)
	## Color and style as above
	#Rawcolors=c(brewer.pal(9,"Set1")[-6])
	#colors=Rawcolors[c(rep(3,4),rep(2,2),rep(1,2),rep(7,3))]
	#MyLty=c(2,3,4,1,2,1,2,1,2,3,1)
	#myLWD=c(3,2,2)

	x=read.table(names[1]);
	len=dim(x)[1]/5
	xvals= x[(len*3+1):(len*4), 1] # this is the frequency, which is in the first column of the SN files
	xvals=ifelse(xvals > 1, 1.0, xvals-ToSubstractFromFreq)
	#idxToPlot=which(numGroundTruthSNVs>=2000 && xvals<=0.08)
	idxToPlot=which(xvals<=0.08)
	#plot_until=10
	
	tuned_SN_pdf=paste(eval_dirs[2],"SN_locRe_binom_allTools_inlet_lowFreqVariants.pdf",sep="")
	if ( ! file.exists(tuned_SN_pdf) ){
		pdf(tuned_SN_pdf,pointsize=3)
		myLWD=c(rep(6,num_tunes))
		block=3
		i=1
		vals = x[(len*(block-1)+1):(len*block), 2]      # this is the sensitivity, which is in the second column of the SN files; Here we take the block "block" of the 
								# sensitivity which is for a certain precision cutoff
		print(vals)
		#plot(cbind(xvals[c(seq(1,plot_until))], vals[c(seq(1,plot_until))]), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i], ylim=c(0,0.55),xlim=c(0.004000,0.046))
		plot(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), type="l", col=colors[i], xlab="Variant Frequency", ylab="Sensitivity",lty=MyLty[i],cex.lab=2.5,cex.axis=2.5,lwd=myLWD[i], ylim=c(0,0.55),xlim=c(0.00,0.05))
		abline(v=xvals[idxToPlot],lty=3,lwd=myLWD[1],col="gray")

		for (i in 2:num_tunes)
		{
			x=read.table(names[i]);
			vals = x[(len*(block-1)+1):(len*block), 2]
			#lines(cbind(xvals[c(seq(1,plot_until))], vals[c(seq(1,plot_until))]), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
			lines(cbind(c(0,xvals[idxToPlot]), c(0,vals[idxToPlot])), col=colors[i], lty=MyLty[i],lwd=myLWD[i])
		}
		dev.off()
	}
	
#quit()
}



if(FIGURE1ACOVERAGE){
	## Like Figure 1A in the paper but taking into account only ground truth SNVs where coverage in normal and tumor is sufficient
	all_eval_dirs=eval_dir_k20
	
	for (curr_eval_dir in all_eval_dirs){
		names=c(paste(curr_eval_dir,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""),
			paste(curr_eval_dir,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_summaryFreq_coveredRegions",sep=""))

		tools=c("deepSNV", "GATK HP", "GATK UG", "JointSNVMix2", "MuTect", "SAMtools", "SiNVICT", "SomaticSniper", "Varscan2")
		num = length(names)
		#colors=rainbow(num)
		#colors=brewer.pal(8,"Accent")
		colors=c(brewer.pal(9,"Set1")[-6])
		colors=colors[c(1,5,4,2,6,3,7,8)]
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
		idxToPlot=which(numGroundTruthSNVs>=200)
	
	
		SN_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_coveredRegions.pdf",sep="")
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
			abline(v=xvals[idxToPlot],lty=3,col="gray")

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
	
			legend("bottomright",inset=c(-1,-1), legend=tools, col=colors, lwd=3)
			
			legend("bottomleft",inset=c(-1,-1),legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY)
			#par(mar=c(5, 4, 4, 2) + 0.1)
			dev.off()
		}
		legend_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_coveredRegions_legend.pdf",sep="")
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
	
	
			legend("bottomright", legend=tools, col=colors, lwd=3,cex=2.5)
			
			legend("bottomleft",legend=c("10%", "5%"),lwd=myLWD,title="False discovery cutoff",lty=myLTY,cex=2.5)
			dev.off()
		}
		SN_inlet_pdf=paste(curr_eval_dir,"SN_all_twoPrec_10_5_coveredRegions_inlet.pdf",sep="")
		if ( ! file.exists(SN_inlet_pdf) ){
			pdf(SN_inlet_pdf,pointsize=3)
			#par(xpd=T, mar=par()$mar+c(6,6,0,6))
			cnt=1
			myLWD=c(6,6,6,6,6,6,6,6,6,6,6,6,6)
			i=1
			plot_untdxToPlotl=10	
			numGroundTruthSNVs=x[(len*3+1):(len*4), 3]
			#idxToPlot=which(numGroundTruthSNVs>=2000 && xvals<=0.08)
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
	print("Done with Figure 1A coveredRegions")
#quit()
} 



