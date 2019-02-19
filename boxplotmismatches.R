library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)

args <- commandArgs(trailingOnly = TRUE)



#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")


spec <- matrix(c(
        'mismatch'     , 'm', 1, "character", "mismatches",
        'samples'    , 's', 1, "character", "sample file (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'directory'     , 'd', 1, "character", "output directory (required)",


        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


mismatches <- read.table(opt$mismatch, header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
sampletable <- read.table(opt$samples)
trnatable <- read.table(opt$trna)
directory <- opt$directory

outputformat <- ".pdf"

#head(mismatches)
#tRNA_name	sample	position	percentmismatch	coverage	tRNAreadstotal	actualbase	mismatchedbases	adenines	thymines	cytosines	guanines

#positions =  c("58","26","37","9","76","34","20","32","1","73","6","49","e9","35","27","5")

#posmismatches = mismatches[,c("position","tRNA_name","percentmismatch")]
totalmism = aggregate(mismatches$percentmismatch,  by=list(position = mismatches$position, tRNA_name =  mismatches$tRNA_name),FUN=max)
#head(totalmism)
posmism = aggregate(totalmism$x,  by=list(position = totalmism$position),FUN=function(mism) {sum(mism > .1)})
#head(posmism)
colnames(posmism) = c("tRNA_position","Mismatched_Transcripts")
write.table(posmism,file = "positionmismatches.txt",row.names = TRUE )
positions = posmism[posmism$Mismatched_Transcripts > 5,"tRNA_position"]

#positions = unique(totalmism[totalmism$x > 10,"position"])

aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","M")

mismatchmelt = mismatches[,c("tRNA_name","sample","percentmismatch","position")]


features = unique(mismatches$tRNA_name)
#print(features)

#relevantmismatches <- list()

relevantmismatches = data.frame("Samples", "position", "tRNA","first","second","difference")
colnames(relevantmismatches) = c("Samples", "position", "tRNA","first","second","difference")


currlistpos = 1



for (currpos in positions){
mismatchmelt = mismatches[mismatches$position == currpos,c("tRNA_name","sample","percentmismatch")]
#print(head(mismatchmelt))

#head(sampletable[match(mismatchmelt$variable,sampletable[,1]),2])


mismatchmeltagg <- aggregate(mismatchmelt$percentmismatch, by=list(Feature = mismatchmelt$tRNA_name, sample = sampletable[match(mismatchmelt$sample,sampletable[,1]),2]), FUN=mean)


#coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = sampletable[match(coveragemelt$Sample,sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)

colnames(mismatchmeltagg) <- c("Feature","sample","percentmismatch")
 
posname = paste(directory,"/",currpos,"-possamplemismatches",outputformat, sep = "")


ggplot(data = mismatchmeltagg, aes(x=sample, y=percentmismatch)) + geom_boxplot(aes(fill=sample)) + geom_jitter(aes(fill=sample), size = .1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1)
ggsave(filename=posname,width=7, height=7)

mismatchmeltaminoagg <- aggregate(mismatchmelt$percentmismatch, by=list(sample = mismatchmelt$sample, amino = trnatable[match(mismatchmelt$tRNA_name,trnatable[,1]),3]), FUN=mean)

#print(head(mismatchmeltaminoagg))
colnames(mismatchmeltaminoagg) <- c("sample","amino","percentmismatch")


posname = paste(directory,"/",currpos,"-posaminomismatches",outputformat, sep = "")



ggplot(data = mismatchmeltaminoagg, aes(x=amino, y=percentmismatch)) + geom_boxplot(aes(fill=amino)) + geom_jitter(aes(fill=amino), size = .5) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1)
ggsave(filename=posname, width=7, height=7)


mismatchbase = mismatches[mismatches$position == currpos,c("tRNA_name","sample")]
mismatchbase$attotal = (mismatches[mismatches$position == currpos,"adenines"]) / mismatches[mismatches$position == currpos,"coverage"]
mismatchbase$cgtotal = (mismatches[mismatches$position == currpos,"cytosines"] ) / mismatches[mismatches$position == currpos,"coverage"]

posname = paste(currpos,"-posbasedistribution",outputformat, sep = "")

#currplot <- qplot(data=mismatchbase,x=mismatchbase$attotal,y=mismatchbase$cgtotal,xlab = "Adenines",ylab = "Cytosines", asp=1)+ theme_bw()
#ggsave(filename=posname)

#print(head(mismatchmeltagg$Feature ))
for (currfeat in features){
#print(currfeat)
comparisons = combn(unique(mismatchmeltagg$sample),2,simplify = FALSE)

mismatchmeltfeat = mismatchmeltagg[mismatchmeltagg$Feature == currfeat,]
#print( comparisons[[1]][1])
#print( comparisons[[1]][2])
#firline = 
#secline = 

#secline = mismatchmeltfeat[mismatchmeltfeat$sample == comparisons[[2]] & mismatchmeltfeat$tRNA_name == currfeat,]
#secline = mismatchmeltfeat



#print(firline)
#print(secline)
#print("***")
for (currcompare in comparisons){
#print(head(mismatchmeltfeat))
#print(unique(mismatchmeltfeat$sample))
#print(currcompare[[1]])
#print(currcompare[[2]])


resultname =  paste(currcompare[[1]],currcompare[[2]] ,sep= ":")
firline = mismatchmeltfeat[mismatchmeltfeat$sample == currcompare[[1]],]
secline = mismatchmeltfeat[mismatchmeltfeat$sample == currcompare[[2]],]


#resultname = sapply(comparisons, function(currcompare){ paste(currcompare[[1]],currcompare[[2]] ,sep= ":")})
#firnum = sapply(comparisons, function(currcompare){ mismatchmeltfeat[mismatchmeltfeat$sample == currcompare[[1]]]})
#secnum = sapply(comparisons, function(currcompare){ mismatchmeltfeat[mismatchmeltfeat$sample == currcompare[[2]]]})
#
#cbind(resultname)


#print(resultname)
#print(currpos)
#print(firline)
#print(secline)
#
#
if(sum(abs(firline$percentmismatch - secline$percentmismatch)) > .1){

#print (c(resultname,currpos,currfeat,firline$percentmismatch,secline$percentmismatch,abs(firline$percentmismatch - secline$percentmismatch)))


#x <- c("Samples", "position", "tRNA","first","second","difference")
#colnames(relevantmismatches) <- x


#relevantmismatches[[currlistpos]] <- t(c(resultname,currpos,currfeat,firline$percentmismatch,secline$percentmismatch,abs(firline$percentmismatch - secline$percentmismatch)))
#currlistpos = currlistpos + 1
#relevantmismatches

currlist <- list(Samples = resultname,position = currpos,tRNA = currfeat,first = firline$percentmismatch,second = secline$percentmismatch,difference =abs(firline$percentmismatch - secline$percentmismatch))

#print(colnames(relevantmismatches))
relevantmismatches = rbind(relevantmismatches,currlist, stringsAsFactors=FALSE)
#print(currlist)
#print(resultname)
#print(currpos)
#print(firline)
#print(secline)
#print("******")
}



#diffrows = lapply(comparisons, function(currresult){mismatchmeltfeat[mismatchmeltfeat$sample[currresult[[2]],] > mismatchmeltfeat$sample[currresult[[2]]] | currresult[[2]]] > mismatchmeltfeat$sample[currresult[[2]]]})
#mismatchmeltfeat$sample['',]
}
}
}
#mismatchmelt = mismatches[,c("tRNA_name","sample","percentmismatch","position")]

#head(mismatchmelt$sample)
#head(mismatchmelt)
#mismatchmelt[mismatchmelt$position == "58" & mismatchmelt$sample == "M_dm_Heart_M5_minusAlkB" & mismatchmelt$tRNA_name == "tRNA-Val-TAC-1",]

mismatchresults <- do.call("rbind",relevantmismatches)
#print (head(mismatchresults))