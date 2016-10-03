
# normalize the scores to be within 0-250
normalize_scores <- function(raw_scores){
	max_score=max(raw_scores)
	cat(paste("The maximum score is: ",max_score,"\n",sep=""))
	normalized_scores=(raw_scores/max_score)*250
	normalized_scores
}

args <- commandArgs(TRUE)
if(length(args)<4){
	cat(paste("Usage: combininglistsbetter_automated.R <outDir> <numTools,n={1,...,9}> <tool1> [<tool2> ... <tooln>]\n",sep=""))
	quit()
}
for(i in 1:length(args)){
	cat(paste(args[i],"\n",sep=""))
	if(i==3 || i==4){
		cat(paste(basename(args[i]),"\n",sep=""))
		print(grepl("deepSNV",basename(args[i])))
		print(grepl("varscan2",basename(args[i])))
	}
}
print(length(args))

outDir=args[1]
numTools <- as.numeric(args[2])
print(outDir)
print(numTools)
tool1=args[3]
print(tool1)
numAdditionalArgs=2
stopifnot(length(args)==numTools+numAdditionalArgs)
toolsToCombine=list()
MuTectIn=FALSE
cntTools=0
typeoftest=""
for(i in 3:length(args)){
	currTool=args[i]
	print(currTool)
	if(grepl("muTect_SNVs_Raw.vcf",basename(currTool))){
		MuTectIn=TRUE
		muTectTool=currTool
		cat("We found muTect!\n")
		typeoftest <- paste(typeoftest,"mutec_",sep="")
	} else {
		cntTools=cntTools+1
		toolsToCombine[[cntTools]]=currTool
		print(grepl("gatk_SNVs",basename(currTool)))
		print(grepl("gatkHPCaller",basename(currTool)))
		print(grepl("deepSNV",basename(currTool)))
		print(grepl("jointSNVMix2",basename(currTool)))
		print(grepl("varscan2",basename(currTool)))
		print(grepl("somaticSniper",basename(currTool)))
		print(grepl("samvar",basename(currTool)))
		if (grepl("gatk_SNVs",basename(currTool))) {
			 typeoftest <- paste(typeoftest,"UG_",sep="")
		} else if (grepl("gatkHPCaller",basename(currTool))) {
			typeoftest <- paste(typeoftest,"HP_",sep="")
		} else if (grepl("deepSNV",basename(currTool))) {
			typeoftest <- paste(typeoftest,"deep_",sep="")
		} else if (grepl("jointSNVMix2",basename(currTool))) {
			typeoftest <- paste(typeoftest,"joint_",sep="")
		} else if (grepl("varscan2",basename(currTool))) {
			typeoftest <- paste(typeoftest,"varsc_",sep="")
		} else if (grepl("somaticSniper",basename(currTool))) {
			typeoftest <- paste(typeoftest,"snipe_",sep="")
		} else if (grepl("samvar",basename(currTool))) {
			typeoftest <- paste(typeoftest,"samv_",sep="")
		} else if (grepl("sinvict", basename(currTool))) {
			typeoftest <- paste(typeoftest,"sinv_",sep="")
		} else {
			cat(basename(currTool))
			cat("\n")
			stop("Error: Cannot recognize this tool")
		}
	}
}



ActualNumTools=numTools
if(MuTectIn){
	numTools=numTools-1
}
if(ActualNumTools==9){ ## all tools - file name would become too large
	typeoftest <- "all_9_"
}

cat(paste("There were ",ActualNumTools," tools provided for combination:\n",sep=""))
for(i in 1:numTools){
	cat(toolsToCombine[[i]])
	cat("\n")
}
if(MuTectIn){
	cat(muTectTool)
	cat("\n")
}


fileToGenerate=paste(outDir,"/combinedlistbetter_",ActualNumTools,"_",typeoftest,"prod.vcf",sep="")
if(file.exists(fileToGenerate)){
	cat(paste(fileToGenerate, " already exists.\n"))
	quit()
}


filestocombine=unlist(toolsToCombine)
if(MuTectIn){
	filestocombine=c(filestocombine,muTectTool)
}

kk=numTools ## number of tools without MuTect
toolstotest=c(seq(1,numTools)) # for looping over the elements of "filestocombine"

debugginglimit<- -1 # set to negative number to do the entire files
usefulcols<-c(2,6)

#orig<-vector("list", numtools) 
orig<-vector("list", ActualNumTools)
kk=numTools
realKK=ActualNumTools


for(toolno in toolstotest){ # we take the results from the different tools
 # score<-read.table(paste(diroflists,commonstart,filestocombine[toolno],sep=""),nrows=debugginglimit)
  score <- read.table(filestocombine[toolno],nrows=debugginglimit)
  print(dim(score))
  scoresimp<-score[,usefulcols] # select the useful columns
  scoresimp[,1]<-paste(score[,1],score[,2],score[,4],score[,5],sep="/") # define a unique id
  scoresimp[,2]<-rank(-score[,usefulcols[2]]) # replace the scores by a rank
  names(scoresimp)<-c("id","rank") # rename the columns
  orig[[toolno]]<-scoresimp # store the result
  print(paste("tool",toolno,"done"))
}

#if(kk==2){
if(MuTectIn){
	#toolno<-1 # do the muTect case
	toolno<-ActualNumTools # muTect will be the last tool in the vector
	usefulcols<-c(2,7) # column seven has the accept/reject value
	score<-read.table(filestocombine[toolno],nrows=debugginglimit)
	scoresimp<-score[,usefulcols] # select the useful columns
	scoresimp[,1]<-paste(score[,1],score[,2],score[,4],score[,5],sep="/") # define a unique id
	print(levels(score[,usefulcols[2]]))
	scoresimp[,2]<-rank(score[,usefulcols[2]]=="REJECT") # replace the scores by a rank
	names(scoresimp)<-c("id","rank") # rename the columns
	orig[[toolno]]<-scoresimp # store the result
	print(paste("tool",toolno,"done"))
}


#if(kk==2){
if(MuTectIn){
	comby<-orig[ActualNumTools]
	for(toolno in (ActualNumTools-1):1){
		comby<-merge(comby,orig[[toolstotest[toolno]]],by="id",all=TRUE,suffixes=c("",toolno)) # join data frames
		print(paste("tool",toolstotest[toolno],"merged"))
	}
} else {
	comby<-orig[toolstotest[1]]
	for(toolno in 2:kk){
		comby<-merge(comby,orig[[toolstotest[toolno]]],by="id",all=TRUE,suffixes=c("",toolno)) # join data frames
		print(paste("tool",toolstotest[toolno],"merged"))
	}
}

NN<-nrow(comby)+1
comby[is.na(comby)] <- NN # set NAs to end of the list

print("data combined")

# transform to a cube

rescalecube<-comby
for(i in 1:(realKK)){ # transform according to equation in notes
  rescalecube[,i+1]<-comby[,i+1]/(2*(NN))
}

print("transformed to cube")

# transform to a normal

rescaley<-rescalecube
for(i in 1:(realKK)){ # treat the p-type values as the CDF of a normal 
  rescaley[,i+1]<-qnorm(rescalecube[,i+1])
}

print("transformed to normal")

# correlations and covariances
covcororig<-cov.wt(rescaley[,1:(realKK)+1],cor=TRUE,center=FALSE) # the center=FALSE option creates the matrix of second moments
covmat<-covcororig$cov
cormat<-covcororig$cor

print("The original matrix of second moments")
print(covmat)

decomp<-svd(covmat) # provides a
Ddiag <- decomp$d # diagonal vector
O <- decomp$u # and an ortogonal matrix

newmat<-cormat
newmat[,]<-(sum(cormat)/(realKK)-1)/(realKK-1)
diag(newmat)<-1
desiredcovmat<-t(t(newmat*sqrt(diag(covmat)))*sqrt(diag(covmat)))

print("A matrix of second moments with equal correlations")
print(desiredcovmat)

decompd<-svd(desiredcovmat) # provides a
Ddiagd <- decompd$d # diagonal vector
Od <- decompd$u # and an ortogonal matrix
sqrtcovmat<-O%*%diag(sqrt(Ddiag))%*%t(O)
sqrtdesiredcovmat<-Od%*%diag(sqrt(Ddiagd))%*%t(Od)
transmat<-solve(sqrtcovmat)%*%sqrtdesiredcovmat

print("A possible transformation matrix")
print(transmat)# a transformation matrix with balanced correlations

# the problem is that the direction of each column of O is undefined

centreofmass<-rep(0,realKK) # find the center of mass
for(i in 1:(realKK)){
  centreofmass[i]<-mean(rescaley[,i+1])
}

transmatnew<-transmat
signtest<-t(as.matrix(centreofmass))%*%transmat #check whether the new center of mass is in the lower quadrant
for(i in 1:(realKK)){
  if(signtest[i]>0){
    transmatnew[,i]<--transmat[,i] # if not swap it
  }
}
print("Check that the transformed mean is in the lower quadrant")
print(t(as.matrix(centreofmass))%*%transmatnew)

decorry<-as.matrix(rescaley[,1:(realKK)+1])%*%transmatnew # this removes the correlation 

print("data decorrelated")

print("new correlations")
newcovs<-cov.wt(decorry,center=FALSE,cor=TRUE)
print(newcovs)

save(covcororig,newcovs,file=paste(outDir,"/covariances",sub("_$","",typeoftest),".rData",sep=""))

# going back to the cube space

decorcube<-comby
for(i in 1:(realKK)){ # now we go back to the CDF of a normal
  decorcube[,i+1]<-pnorm(decorry[,i])
}

print("transformed to cube")

sumrank<-rank(-rowSums(decorcube[,-1])) # this is the sum of the transformed scores then ranked
prodrank<-rank(-rowSums(log(decorcube[,-1]))) # this is effectively the product of the transformed scores then ranked

decorcube$prodrank<-prodrank

print(dim(decorcube))

if(ActualNumTools==3){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,2:5]) # just a matrix with 9 columns, because we want the VCF format

} else if(ActualNumTools==2){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,2:4],decorcube[,3:4]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==4){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,3:5]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==5){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,4:5]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==6){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,3]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==7){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==8){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-decorcube[,-2] # just a matrix with 9 columns, because we want the VCF format
} else if(ActualNumTools==9){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-decorcube[,c(-2,-3)] # just a matrix with 9 columns, because we want the VCF format
} else {
	stop("Error: This number of tools is not implemented!")
}

names(toprintasVCFfile)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") # rename the columns
## the intermediate vcf file is just for checking and not for evaluating!
write.table(toprintasVCFfile,file=paste(outDir,"/combinedlistbetter_",ActualNumTools,"_",typeoftest,"intermediate.vcf",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

print("sum data saved")

toprintasVCFfile[,c(3,7:9)]<-"." ## now we make the vcf format by setting the 3rd, 7th,8th, and 9th column to "."
toprintasVCFfile[,6]<-normalize_scores(sumrank)

toprintasVCFfile[,c(1:2,4:5)]<-matrix(unlist(strsplit(comby[,1], "/")),ncol=4,byrow=TRUE) ## the id is of the format "chr10/100003638/T/C", and here we split it up and assign the different fields to column 1, 2, 4 and 5
names(toprintasVCFfile)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") # rename the columns

write.table(toprintasVCFfile,file=paste(outDir,"/combinedlistbetter_",ActualNumTools,"_",typeoftest,"sum.vcf",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

print("sum data saved")

toprintasVCFfile[,6]<-normalize_scores(prodrank)

write.table(toprintasVCFfile,file=paste(outDir,"/combinedlistbetter_",ActualNumTools,"_",typeoftest,"prod.vcf",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

print("prod data saved")


sampley<-sample.int((NN-1),1000)
rescalecubetest<-rescalecube[sampley,]
save(rescalecubetest,file=paste(outDir,"/rescalecubebettertest_",sub("_$","",typeoftest),".rData",sep=""))
rescaleytest<-rescaley[sampley,]
save(rescaleytest,file=paste(outDir,"/rescaleybettertest_",sub("_$","",typeoftest),".rData",sep=""))
decorrytest<-decorry[sampley,]
save(decorrytest,file=paste(outDir,"/decorrybettertest_",sub("_$","",typeoftest),".rData",sep=""))
decorcubetest<-decorcube[sampley,]
save(decorcubetest,file=paste(outDir,"/decorcubebettertest_",sub("_$","",typeoftest),".rData",sep=""))

print("test data saved and finished")

