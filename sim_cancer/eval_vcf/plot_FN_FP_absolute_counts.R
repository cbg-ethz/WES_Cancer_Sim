#!/usr/local/beerenwinkel/R-3.3.0/bin/Rscript
.libPaths()

args <- commandArgs(trailingOnly = F)
script.dir <- dirname(sub("--file=","",args[grep("--file",args)]))
source(paste(script.dir,"/../../paths.R",sep=""))

library("RColorBrewer")
library("corrplot")
library("gplots")

algmt_k20_dir=paste(base_dir,"alignments_vsn_k20/",sep="")

align=c("alignments_vsn_k20","alignments_sn_k1") 
num=13

barplot_colors=c(brewer.pal(11,"RdYlGn")[1],
	brewer.pal(11,"RdYlGn")[4],
	brewer.pal(11,"RdYlGn")[7],
	brewer.pal(11,"RdYlGn")[10],
	brewer.pal(11,"PRGn")[2],
	brewer.pal(11,"PiYG")[2],
	brewer.pal(11,"RdYlBu")[8],
	brewer.pal(9,"Set1")[6],
	brewer.pal(11,"RdYlBu")[11],
	brewer.pal(11,"RdGy")[8])


for (j in 1:length(align))
{
	eval_dir=paste(base_dir, align[j], "/variants_default/eval_10000_160527/", sep="")

	tools=c("perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf", "perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf", "perc.gatk_SNVs.raw.SOMATIC.vcf", "perc.jointSNVMix2_SNVs_Raw.vcf", "perc.muTect_SNVs_Raw.vcf", "perc.option__noE_samvar_1_2.SOMATIC.vcf.gz", "perc.somaticSniper_SNVs_Raw_qual.noComma.vcf", "perc.bam__varscan2.txt.snp.Somatic_qual.vcf")

	names=c("deepSNV","GATK HP","GATK UG","JointSNVMix2","MuTect","SAMtools", "SomaticSniper", "VarScan2")

	FP_names=c()
	FN_names=c()
	
	for (i in 1:(length(tools)))
	{
		fn_dat_fp=paste(eval_dir, "/TU.wCont20.final.RG.50", tools[i], "_indel0.eval_FP_absolute_counts", sep="") 
		fn_dat_fn=paste(eval_dir, "/TU.wCont20.final.RG.50", tools[i], "_indel0.eval_FN_absolute_counts", sep="")
		dat_fp = read.table(fn_dat_fp)
		dat_fn = read.table(fn_dat_fn)
	
		fn_dat_fp_total=paste(eval_dir, "/TU.wCont20.final.RG.50", tools[i], "_indel0.eval_FP_pie_fix", sep="") # for the total number of FP's
		fn_dat_fn_total=paste(eval_dir, "/TU.wCont20.final.RG.50", tools[i], "_indel0.eval_FN_pie_fix", sep="") # for the total number of FN's
		dat_fp_total=read.table(fn_dat_fp_total)
		dat_fn_total=read.table(fn_dat_fn_total)

		FN_idx = -1  # not the first idx, because we are interested in the high freq false negatives
		x_fn = dat_fn[, FN_idx] 
		
		#row = which(dat_fp[ ,1]==0.99)
		row = which(dat_fp[ ,1]==0.950000)
		x_fp = dat_fp[row, -1]
		colors=c(brewer.pal(12,"Paired")[1:10],brewer.pal(8,"Greys")[6]) # length 11 for low_cov, low_qual, near_indel, low_support, filtered_normal, low_qual_normal, low_cov_normal, in_normal, not_in_true, very_high_cov,  other

		labels_fp=c("low coverage", "low quality", "variable region","low quality in normal","low coverage normal", "in normal", "alignment location", "high coverage", "other")
		##  all_counts[var_near_indel] 1, all_counts[in_normal] 2, all_counts[not_in_true]3, all_counts[very_high_cov]4, all_counts[other]5,all_counts[low_qual]6,all_counts[low_qual_normal]7,all_counts[low_cov]8,all_counts[low_cov_normal]9);
		FP_idx=c(8,6,1,7,9,2,3,4,5)
		
		labels_FN=c("low frequency", "low coverage", "low quality", "variable region", "low support", "seq error in normal", "low quality in normal", "low coverage in normal", "other")
		colors_FP=colors[c(1,2,3,6,7,8,9,10,11)]
		colors_FN=colors[c(1,2,3,4,5,6,7,11)]

		FN_names = c(FN_names, paste(names[i], sum(dat_fn_total[, FN_idx]), sep=" ")) # for the total number of FN's
		FP_names = c(FP_names, paste(names[i], sum(dat_fp_total[row,-1]), sep=" ")) # for the total number of FP's


		# get the normalized number of false positives and negatives - divide the absolute number in each category by the total number of fn's/fp's
		if (sum(dat_fn_total[,FN_idx])==0)
		{
			dat_fn_norm=rep(0.0001, length(dat_fn_total[,FN_idx]))
		}
		else
		{
			dat_fn_norm = (x_fn/sum(dat_fn_total[,FN_idx]))*100
		}
		if (sum(dat_fp_total[row,-1])==0)
		{
			dat_fp_norm=rep(0.0001, length(dat_fp_total[row,-1]))
		}
		else
		{
			dat_fp_norm = (x_fp/sum(dat_fp_total[row, -1]))*100
		}


		if (i==1)
		{
			#all_dat_fp = x_fp;
			#all_dat_fn = x_fn; 
			all_dat_fp = dat_fp_norm;
			all_dat_fn = dat_fn_norm;
			total_fp=sum(dat_fp_total[row,-1])
			total_fn=sum(dat_fn_total[, FN_idx])
		}
		else
		{
			#all_dat_fp = rbind(all_dat_fp, x_fp) 
			#all_dat_fn = rbind(all_dat_fn, x_fn);
			all_dat_fp = rbind(all_dat_fp, dat_fp_norm)
			all_dat_fn = rbind(all_dat_fn, dat_fn_norm);
			total_fp=c(total_fp, sum(dat_fp_total[row,-1])) 
			total_fn=c(total_fn, sum(dat_fn_total[, FN_idx]))
		}
	}

print(all_dat_fn) # each row is a tool
	ncol=dim(all_dat_fn)[2]
print(ncol)
	dataToPlotFN=ifelse(t(matrix(unlist(all_dat_fn),ncol=ncol)) > 0, t(matrix(unlist(all_dat_fn),ncol=ncol)), 1e-10) # each column is a tool

	ncol=dim(all_dat_fp[,FP_idx])[2]# each row is a tool
	dataToPlotFP=ifelse(t(matrix(unlist(all_dat_fp[,FP_idx]),ncol=ncol))>0, t(matrix(unlist(all_dat_fp[,FP_idx]),ncol=ncol)), 1e-10) # if all entries are zero for a caller, like in the case of MuTect, it does not work with the correlation matrix # each column is a tool
	print(dataToPlotFP)
numTools=length(tools)
stopifnot(numTools==8)

	fn_pdf=paste(eval_dir, "/FN_absolute_", align[j], "_color.pdf", sep="")
	if ( ! file.exists(fn_pdf)){
		pdf(fn_pdf)
		print(fn_pdf)
		ncol=dim(all_dat_fn)[2]
		#barplot(dataToPlot, legend.text=labels_FN[FN_idx], col=colors_FN, names.arg=FN_names, horiz=FALSE, las=2, beside=TRUE,ylim=c(1,max(total_fn)+10000),log="y",space=c(0,3),yxt="n")
		barplot(dataToPlotFN, legend.text=labels_FN[FN_idx], col=colors_FN, names.arg=FN_names, horiz=FALSE, las=2, beside=TRUE,ylim=c(0,100),space=c(0,3),yxt="n")
		axis(2,at=c(0,20,40,60,80,100),labels=c(paste(c(0,20,40,60,80,100),"%",sep="")),las=2)
		abline(h=c(0,20,40,60,80,100),lty=2,col="lightgray")
		dev.off()
	}

	fn_correlation_plot=paste(eval_dir, "/FN_errorProfile_correlation_", align[j], "_color.pdf", sep="")
	if(!file.exists(fn_correlation_plot)){
		MatrixToPlot=dataToPlotFN
		print(MatrixToPlot)
		colnames(MatrixToPlot)=names
		corMatrixToPlot=cor(MatrixToPlot)
		stopifnot(numTools==dim(corMatrixToPlot)[1])
		stopifnot(numTools==dim(corMatrixToPlot)[2])
		print(corMatrixToPlot)
		if(!any(is.na(corMatrixToPlot))){
			pdf(fn_correlation_plot)
			dim_corr_matrix=dim(corMatrixToPlot)[1]
        		corrplot(corMatrixToPlot, type='lower',tl.col="black")
			dev.off()
		}
	}
		
	fn_pdf=paste(eval_dir, "/FP_absolute_", align[j], "_color.pdf", sep="")
	if( ! file.exists(fn_pdf)){
		pdf(fn_pdf)
		print(fn_pdf)
		ncol=dim(all_dat_fp[,FP_idx])[2]
		barplot(dataToPlotFP, legend.text=labels_fp, col=colors_FP, names.arg=FP_names, horiz=FALSE, las=2,beside=TRUE,ylim=c(0,100),space=c(0,3),yxt="n")
		axis(2,at=c(0,20,40,60,80,100),labels=c(paste(c(0,20,40,60,80,100),"%",sep="")),las=2)
		abline(h=c(0,20,40,60,80,100),lty=2,col="lightgray")
		dev.off()
	}

	fp_correlation_plot=paste(eval_dir, "/FP_errorProfile_correlation_", align[j], "_color.pdf", sep="")
	if(!file.exists(fp_correlation_plot)){
		MatrixToPlot=dataToPlotFP
		print(MatrixToPlot)
		colnames(MatrixToPlot)=names
		corMatrixToPlot=cor(MatrixToPlot)
		stopifnot(numTools==dim(corMatrixToPlot)[1])
		stopifnot(numTools==dim(corMatrixToPlot)[2])
		print(corMatrixToPlot)
		if(!any(is.na(corMatrixToPlot))){
			pdf(fp_correlation_plot)
			dim_corr_matrix=dim(corMatrixToPlot)[1]
			corrplot(corMatrixToPlot, type='lower',tl.col="black")
	        	dev.off()
		} else {
			cat(paste("Do not plot correlation profile ",fp_correlation_plot,"\n",sep=""))
		}
	}
}
