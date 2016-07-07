#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("RColorBrewer")
args <- commandArgs(trailingOnly = F)
script.dir <- dirname(sub("--file=","",args[grep("--file",args)]))
source(paste(script.dir,"/../../paths.R",sep=""))

pdfToPlot <- paste(eval_dir_k20,"/Summary_auPRC_FDR10_barplot.pdf",sep="")

names=c(paste(eval_dir_combi,"combinedlistbetter_5_deepSNV_varscan2_samvar_1_2_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_combi,"combinedlistbetter_4_deepSNV_varscan2_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_combi,"combinedlistbetter_5_deepSNV_varscan2_somaticSniper_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_binom,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_binom,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_tuned,"TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_tuned,"TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_auPRC_10",sep="")
	)
	
print(names)

if( !file.exists(pdfToPlot)){

		all_auPRC_10s=c()
		maxval=0
		for (i in 1:length(names)){
			print(names[i])
			all_auPRC_10s=c(all_auPRC_10s,as.numeric(read.table(names[i])[,1]))			
		}
		print(all_auPRC_10s)	
		maxval=max(all_auPRC_10s)

		RawColors1=c(brewer.pal(9,"Set1")[-6])
		RawColors2=c(rep(c(brewer.pal(8,"Accent")[6], brewer.pal(11,"PiYG")[8], brewer.pal(9,"Oranges")[4], brewer.pal(9,"GnBu")[6]),2))
		colors=c(RawColors2[c(1,2,3)],RawColors1[c(1,5,4,2,6,3,7,8)])
		colors=adjustcolor(colors,alpha.f =0.7)

		pdf(pdfToPlot)
		barplot(all_auPRC_10s,col=colors,main="Direct comparison combinations and individual tools", ylab="Area under Precision Recall Curve (10% FDR)",xaxt="n",ylim=c(0,0.8))
		axis(1,at=c(seq(1,length(names))),labels=sub(".vcf.gz_indel0.eval_auPRC_10","",sub(".vcf_indel0.eval_auPRC_10","",sub("TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.","",sub("TU.wCont20.final.RG.50perc.","",sub("combinedlistbetter","",names))))),cex.lab=0.5,cex.axis=0.5,las=2)
		dev.off()

		pdfToPlot=sub(".pdf","_legend.pdf",pdfToPlot)
		if( !file.exists(pdfToPlot)){
			pdf(pdfToPlot)
			plot(1,1)
			legend("topright", legend=basename(sub("_indel0.eval_auPRC_10","",sub("TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.","",sub("TU.wCont20.final.RG.50perc.","",sub("combinedlistbetter_","",names))))), fill=colors)
			dev.off()
		}
}


pdfToPlot <- paste(eval_dir_k20,"/Summary_auPRC_FDR10_barplot_comparisonDefault.pdf",sep="")

names=c(paste(eval_dir_combi,"combinedlistbetter_5_deepSNV_varscan2_samvar_1_2_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_combi,"combinedlistbetter_4_deepSNV_varscan2_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_combi,"combinedlistbetter_5_deepSNV_varscan2_somaticSniper_muTect_jointSNVMix2_prod.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_binom,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_binom,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_tuned,"TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_tuned,"TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_auPRC_10",sep=""),
	paste(eval_dir_k20,"TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_auPRC_10",sep="")
	)
	
print(names)

shadingValue=18
nonshadingValue=0.5

if( !file.exists(pdfToPlot)){

		all_auPRC_10s=c()
		maxval=0
		for (i in 1:length(names)){
			print(names[i])
			all_auPRC_10s=c(all_auPRC_10s,as.numeric(read.table(names[i])[,1]))			
		}
		print(all_auPRC_10s)	
		maxval=max(all_auPRC_10s)

		RawColors1=c(brewer.pal(9,"Set1")[-6])
		RawColors2=c(rep(c(brewer.pal(8,"Accent")[6], brewer.pal(11,"PiYG")[8], brewer.pal(9,"Oranges")[4], brewer.pal(9,"GnBu")[6]),2))
		colors=c(RawColors2[c(1,2,3)],RawColors1[c(1,1,5,4,2,6,6,3,3,7,8,8)])
		colors=adjustcolor(colors,alpha.f =0.7)

		pdf(pdfToPlot)
		barplot(all_auPRC_10s,col=colors,main="Direct comparison combinations and individual tools", ylab="Area under Precision Recall Curve (10% FDR)",xaxt="n",ylim=c(0,0.8))
		axis(1,at=c(seq(1,length(names))),labels=sub(".vcf.gz_indel0.eval_auPRC_10","",sub(".vcf_indel0.eval_auPRC_10","",sub("TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.","",sub("TU.wCont20.final.RG.50perc.","",sub("combinedlistbetter","",names))))),cex.lab=0.5,cex.axis=0.5,las=2)
		barplot(all_auPRC_10s,col="black", angle=45,density=c(rep(nonshadingValue,4),shadingValue,nonshadingValue,nonshadingValue,nonshadingValue,nonshadingValue,shadingValue,nonshadingValue,shadingValue,nonshadingValue,nonshadingValue,shadingValue), add=T)
		dev.off()

		pdfToPlot=sub(".pdf","_legend.pdf",pdfToPlot)
		if( !file.exists(pdfToPlot)){
			pdf(pdfToPlot)
			plot(1,1)
			legend("topright", legend=basename(sub("_indel0.eval_auPRC_10","",sub("TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.","",sub("TU.wCont20.final.RG.50perc.","",sub("combinedlistbetter_","",names))))), fill=colors,density=c(rep(nonshadingValue,4),shadingValue,nonshadingValue,nonshadingValue,nonshadingValue,nonshadingValue,shadingValue,nonshadingValue,shadingValue,nonshadingValue,nonshadingValue,shadingValue))
					##density=c(rep(-1,4),shadingValue,-1,-1,-1,-1,shadingValue,-1,shadingValue,-1,-1,shadingValue))
			dev.off()
		}
}
