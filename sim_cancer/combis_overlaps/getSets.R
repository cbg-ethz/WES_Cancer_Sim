#!/usr/local/beerenwinkel/Rbase/R-3.1.2/bin/Rscript
.libPaths(c("/usr/local/beerenwinkel/R/x86_64-redhat-linux-gnu-library/3.1.2", .libPaths()))
.libPaths()

library("VennDiagram")
library("RColorBrewer")

args <- commandArgs(TRUE)
allPos=list()
namesLists=c()
numLists=5
outDir=dirname(args[1])

for (i in 1:numLists){
	allPos[[i]] <- read.table(args[i])
	names(allPos)[i] <- sub("_pos.txt","",basename(args[i]))
	namesLists=c(namesLists,sub("_pos.txt","",basename(args[i])))
}


# intersectPairs=list()
# cnt=1
# for (i in 1:numLists){
# 	if(i==numLists){
# 		break
# 	}
# 	for (j in (i+1):numLists){
# 		intersectPairs[[cnt]] <- intersect(unlist(allPos[[i]]),unlist(allPos[[j]]))
# 		names(intersectPairs)[cnt] <- paste(namesLists[i],namesLists[j],sep="_")
# 		cnt=cnt+1
# 	}
# }
# print(names(intersectPairs))


# intersectAll:
intersectAll=intersect(unlist(allPos[[1]]), unlist(allPos[[2]]))
for (i in 3:numLists){
	intersectAll=intersect(intersectAll, unlist(allPos[[i]]))
}
print(length(intersectAll))


## get specific overlaps/intersections

## one tool private
for (tool in namesLists){
	indexTool=which(namesLists==tool)
	print(namesLists[indexTool])
	# union of all others
	allOtherButTool=c()
	for (i in 1:numLists){
		if (i!=indexTool){
			allOtherButTool=union(allOtherButTool, unlist(allPos[[i]]))
		}
	}
	print(head(allOtherButTool))
	private_tool <- setdiff(unlist(allPos[[indexTool]]),allOtherButTool)
	print(length(private_tool))
	write.table(private_tool, file=paste(outDir,"/",namesLists[indexTool],"_private_pos.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
}


## two tools private
for (i in 1:numLists){
	if(i==numLists){
		break
	}
	for (j  in (i+1):numLists){
		# intersection of the two
		intersectTwo <- intersect(unlist(allPos[[i]]), unlist(allPos[[j]]))

		# union of all others
		allOtherButTool=c()
		for (k in 1:numLists){
			if (k!=i && k!=j){
				allOtherButTool=union(allOtherButTool, unlist(allPos[[k]]))
			}
		}
		private_tool <- setdiff(intersectTwo,allOtherButTool)
		print(length(private_tool))
		write.table(private_tool, file=paste(outDir,"/",namesLists[i],"_overlap_",namesLists[j],"_private_pos.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
	}
}

## three tools private
for (i in 1:numLists){
	if(i==numLists){
		break
	}
	for (j  in (i+1):numLists){
		if(j==numLists){
			break
		}
		for (k in (j+1):numLists){
			# intersection of the three
			intersectThree <- intersect(unlist(allPos[[i]]), intersect(unlist(allPos[[j]]),unlist(allPos[[k]])) )


			# union of all others
			allOtherButTool=c()
			for (m in 1:numLists){
				if (m!=i && m!=j && m!=k){
					allOtherButTool=union(allOtherButTool, unlist(allPos[[m]]))
				}
			}
			private_tool <- setdiff(intersectThree,allOtherButTool)
			print(length(private_tool))
			write.table(private_tool, file=paste(outDir,"/",namesLists[i],"_overlap_",namesLists[j],"_overlap_",namesLists[k],"_private_pos.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
		}
	}
}



## four tools private
for (i in 1:numLists){
        if(i==numLists){
                break
        }
        for (j  in (i+1):numLists){
                if(j==numLists){
                        break
                }
                for (k in (j+1):numLists){
			if(k==numLists){
				break
			}
			for(l in (k+1):numLists){
	                        # intersection of the four
	                        intersectFour <- intersect(unlist(allPos[[i]]), intersect(unlist(allPos[[j]]),intersect(unlist(allPos[[k]]),unlist(allPos[[l]])) ))

	                        # union of all others
	                        allOtherButTool=c()
	                        for (m in 1:numLists){
	                                if (m!=i && m!=j && m!=k && m!=l){
	                                        allOtherButTool=union(allOtherButTool, unlist(allPos[[m]]))
	                                }
	                        }
	                        private_tool <- setdiff(intersectFour,allOtherButTool)
	                        print(length(private_tool))
				write.table(private_tool, file=paste(outDir,"/",namesLists[i],"_overlap_",namesLists[j],"_overlap_",namesLists[k],"_overlap_",namesLists[l],"_private_pos.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
			}
                }
        }
}

## intersect all five:
intersectFive <- intersect(unlist(allPos[[1]]), intersect(unlist(allPos[[2]]),intersect(unlist(allPos[[3]]),intersect(unlist(allPos[[4]]), unlist(allPos[[5]]))) ))
print(length(intersectFive))
write.table(intersectFive, file=paste(outDir,"/intersectAll.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)




