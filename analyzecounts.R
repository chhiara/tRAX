

library("DESeq2")



colgetlogname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colgetavgname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colrename =  function(currtable, newname){

newtable = currtable[,c(5),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}


args = commandArgs(trailingOnly = TRUE)

experimentname = args[1]
inputtable = args[2]
samplefile = args[3]


readcounts = read.table(inputtable,check.names=FALSE)
sampledata = read.table(samplefile,check.names=FALSE)

#sampleinfo = as.character(sampledata[colnames(readcounts) == sampledata[,1],2])
sampleinfo = as.character(sampledata[colnames(readcounts) == gsub("-", ".", sampledata[,1]) ,2])


#gsub("-", ".", sampledata[,1]) 
#sampledata[,2]
#colnames(readcounts)
#sampledata[,1]

samplenames = unique(sampleinfo)
#combn(unique(sampleinfo),2,simplify = FALSE)
#read.table(args[4], stringsAsFactors = FALSE)
#length(args)
#args[4]

if (length(args) > 3){

pairtable = read.table(args[4], stringsAsFactors = FALSE)
pairreduce = pairtable[pairtable[,1] %in% samplenames & pairtable[,2] %in% samplenames,]
comparisons <- apply(pairreduce, 1, list)
comparisons <- lapply(comparisons,unlist)
print(comparisons)


}else{
comparisons = combn(unique(sampleinfo),2,simplify = FALSE)
}


coldata = data.frame(condition=factor(sampleinfo))

cds = DESeqDataSetFromMatrix(countData = readcounts,coldata  ,design = ~ condition)
cds = DESeq(cds,betaPrior=TRUE)





names = lapply(comparisons, function(currcompare){ })

compareresults = lapply(comparisons, function(currcompare){ list(paste(currcompare[[1]],currcompare[[2]] ,sep= ":"),results( cds, contrast=c("condition", currcompare[[1]] ,currcompare[[2]]),cooksCutoff  =TRUE))})


reslist = lapply(compareresults, function(currresult){colrename(currresult[[2]],currresult[[1]])})

resloglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})

resavglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})                                


dds = cds



#print adjusted p-values
allprobs = Reduce(function(x,y) cbind(x,y), reslist)
write.table(allprobs,paste(experimentname,"/",experimentname,"-padjs.txt", sep = ""),sep="	")
                                                                   
#Print log values
alllogvals = Reduce(function(x,y) cbind(x,y), resloglist)
write.table(alllogvals,paste(experimentname,"/",experimentname,"-logvals.txt", sep = ""),sep="	")

#Print out the size factors
write.table(rbind(colnames(readcounts),dds$sizeFactor),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)

#get deseq normalized  raw counts
normalizedrnas = sweep(readcounts,2,dds$sizeFactor, "/" )
write.table(normalizedrnas,paste(experimentname,"/",experimentname,"-normalized.txt", sep = ""), sep = "\t")


