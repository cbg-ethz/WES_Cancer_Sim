#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

args <- commandArgs(trailingOnly = F)
script.dir <- dirname(sub("--file=","",args[grep("--file",args)]))
source(paste(script.dir,"/../../paths.R",sep=""))

library("VennDiagram")
library("RColorBrewer")

callers_selection="top"

##display.brewer.all()

overlaps_file=c(paste(overlapDir,"/SNV_lists_fdr01_1percent/overlaps/overlap_summary_",callers_selection,"5callers.txt",sep=""),paste(overlapDir,"/SNV_lists_fdr05_5percent/overlaps/overlap_summary_",callers_selection,"5callers.txt",sep=""),paste(overlapDir,"/SNV_lists_fdr10_10percent/overlaps/overlap_summary_",callers_selection,"5callers.txt",sep=""))
overlaps_file_name=c("1percent","5percent","10percent")
tools=c("deepSNV", "JointSNVMix2", "SAMtools", "GATK UG", "GATK HP" , "MuTect", "SiNVICT", "SomaticSniper", "Varscan2")

if (callers_selection == "top"){
	## c(brewer.pal(9,"Set1")[-6]) colors for deep, joint, SAMvar, gatkUG, gatkHP, muTect, somSniper, varscan2
	## now in the order of the tools that I have: ('joint' 'deep' 'gatkUG' 'somSniper' 'SiNVICT') # top 5 callers
	#cols=c(brewer.pal(9,"Set1")[-6])[c(2,1,4,7,5)]
	# colors for deep, gatkHP, gatkUG, joint, muTect, SAMvar, SiNVICT, somSniper, varscan2
	colors=c(brewer.pal(9,"Set1")[c(-2,-6)], brewer.pal(9, "Blues")[c(4,8)])[c(1,4,3,9,5,2,8,6,7)][c(4,1,3,8,7)]
	cols=colors
	tools_venn=tools[c(2,1,4,8,7)]
} 
#else if (callers_selection == "bottom") {
#	## now in the order of the tools that I have: 'somSniper' 'gatkHP' 'varscan2' 'SAMvar' 'muTect') # bottom 5 callers
#	cols=c(brewer.pal(9,"Set1")[-6])[c(7,5,8,3,6)]
#	tools_venn=tools[c(7,5,8,3,6)]	
#}

cnt=1
for (curr_overlaps_file in overlaps_file){
	overlaps <- read.table(curr_overlaps_file)
	overlaps <- overlaps[,1]
	overlaps
	length(overlaps)
	caller_names=tools_venn
print(overlaps)
	pdfToPlot=paste(dirname(curr_overlaps_file),"/VennDiagram_overlaps_summary_",overlaps_file_name[cnt],"_",callers_selection,"5callers.pdf",sep="")
	if(!exists(pdfToPlot)){
		pdf(pdfToPlot)
		draw.quintuple.venn(overlaps[1],
                    overlaps[2],
                    overlaps[3],
                    overlaps[4],
                    overlaps[5],
                    overlaps[6],
                    overlaps[7],
                    overlaps[8],
                    overlaps[9],
                    overlaps[10],
                    overlaps[11],
                    overlaps[12],
                    overlaps[13],
                    overlaps[14],
                    overlaps[15],
                    overlaps[16],
                    overlaps[17],
                    overlaps[18],
                    overlaps[19],
                    overlaps[20],
                    overlaps[21],
                    overlaps[22],
                    overlaps[23],
                    overlaps[24],
                    overlaps[25],
                    overlaps[26],
                    overlaps[27],
                    overlaps[28],
                    overlaps[29],
                    overlaps[30],
                    overlaps[31], category = caller_names,
                    lwd = rep(2, 5),
                    lty = rep("solid", 5), 
                    col = cols, fill = cols, alpha = rep(0.4, 5),
                    label.col = rep("black", 31), cex = rep(0.75, 31), fontface = rep("plain", 31),
                    fontfamily = rep("serif", 31), cat.pos = c(0, 287.5, 215, 145, 70),
                    cat.dist = rep(0.25, 5), cat.col = rep("black", 5), cat.cex = rep(1, 5),
                    cat.fontface = rep("plain", 5), cat.fontfamily = rep("serif", 5),
                    cat.just = rep(list(c(0.5, 0.5)), 5), rotation.degree = 0,
		  #print.mode="percent",
                    rotation.centre = c(0.5, 0.5), ind = TRUE)
		dev.off()
	}
	cnt=cnt+1
}

