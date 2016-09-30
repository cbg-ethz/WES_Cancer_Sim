#!/usr/local/beerenwinkel/Rbase/R-3.0.1/bin/Rscript
.libPaths(new = c("/home/arhofman/R/R_3.0.1", .libPaths()))

args <- commandArgs(TRUE)
coverageFile <- args[1]
coverage <- read.table(coverageFile);

a <- as.numeric(as.character(coverage[,1]))
a <- a[!is.na(a)]

cat(min(a), '\tminimum coverage of the gene tracks\n')
cat(max(a), '\tmaximum coverage of the gene tracks\n')
cat(mean(a), '\tmean coverage of the gene tracks\n')
cat(median(a), '\tmedian coverage of the gene tracks\n')


## find out which percentage of the gene tracks is not covered at all
zero_cov <- ( length(which(a==0))/length(a) )*100
cat(zero_cov, '\tpercentage of bases in gene tracks that are not covered at all\n')

atLeast20_cov <- ( length(which(a>=20))/length(a) )*100
cat(atLeast20_cov, '\tpercentage of bases in gene tracks that are covered with at least 20 reads\n')

for (i in 1:500){
	atLeastX_cov <- ( length(which(a>=i))/length(a) )*100
	cat(paste(atLeastX_cov,'\t',i,'\tpercentage of bases in gene tracks that are covered with at least ',i,' reads\n',sep=''))
}

