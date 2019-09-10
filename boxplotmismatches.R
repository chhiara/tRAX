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
        'samples'    , 's', 1, "character", "Sample file (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'directory'     , 'd', 1, "character", "output directory (required)",


        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


mismatches <- read.table(opt$mismatch, header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
Sampletable <- read.table(opt$samples)
trnatable <- read.table(opt$trna)
directory <- opt$directory

outputformat <- ".pdf"

#head(mismatches)
#Feature	Sample	position	percentmismatch	coverage	tRNAreadstotal	actualbase	mismatchedbases	adenines	thymines	cytosines	guanines

#positions =  c("58","26","37","9","76","34","20","32","1","73","6","49","e9","35","27","5")

#posmismatches = mismatches[,c("position","Feature","percentmismatch")]
mismatches$percentmismatch = mismatches$mismatchedbases / (mismatches$coverage + 10)

totalmism = aggregate(mismatches$percentmismatch,  by=list(position = mismatches$position, Feature =  mismatches$Feature),FUN=max)



#head(totalmism)
posmism = aggregate(totalmism$x,  by=list(position = totalmism$position),FUN=function(mism) {sum(mism > .1)})
#head(posmism)
colnames(posmism) = c("tRNA_position","Mismatched_Transcripts")
write.table(posmism,file = "positionmismatches.txt",row.names = TRUE )
mismatchpositions = posmism[posmism$Mismatched_Transcripts > 5,"tRNA_position"]




totaldelete = aggregate(mismatches$deletions/(mismatches$deletions + mismatches$coverage + 30),  by=list(position = mismatches$position, Feature =  mismatches$Feature),FUN=max)
#head(totalmism)
posdelete= aggregate(totaldelete$x,  by=list(position = totaldelete$position),FUN=function(mism) {sum(mism > .1)})
#head(posmism)
colnames(posdelete) = c("tRNA_position","Deleted_Transcripts")
write.table(posdelete,file = "positiondeletions.txt",row.names = TRUE )
deletepositions = posdelete[posdelete$Deleted_Transcripts > 5,"tRNA_position"]

#print(mismatchpositions)
#print(deletepositions)
positions = union(mismatchpositions,deletepositions )
#print ("|*|")

#positions = unique(totalmism[totalmism$x > 10,"position"])

aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","M")

mismatchmelt = mismatches[,c("Feature","Sample","percentmismatch","position")]


features = unique(mismatches$Feature)
#print(features)

#relevantmismatches <- list()

relevantmismatches = data.frame("Samples", "position", "tRNA","first","second","difference")
colnames(relevantmismatches) = c("Samples", "position", "tRNA","first","second","difference")


currlistpos = 1


dotsize = .4

for (currpos in positions){
#print(currpos)

mismatchmelt = mismatches[mismatches$position == currpos,c("Feature","Sample","percentmismatch")]
deletionmelt = mismatches[mismatches$position == currpos,c("Feature","Sample","deletions","coverage")]
#print(head(deletionmelt))
deletionmelt$deletepercent = deletionmelt$deletions/(deletionmelt$deletions + deletionmelt$coverage + 30)
#deletionmelt$deletepercent = deletionmelt$deletions 

#print(deletionmelt$deletions / (deletionmelt$deletions+deletionmelt$coverage))
#print(deletionmelt$deletions / (deletionmelt$deletions+deletionmelt$coverage))


#print(head(deletionmelt))

#head(Sampletable[match(mismatchmelt$variable,Sampletable[,1]),2])


mismatchmeltagg <- aggregate(mismatchmelt$percentmismatch, by=list(Feature = mismatchmelt$Feature, Sample = Sampletable[match(mismatchmelt$Sample,Sampletable[,1]),2]), FUN=mean)


#coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = Sampletable[match(coveragemelt$Sample,Sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)

colnames(mismatchmeltagg) <- c("Feature","Sample","percentmismatch")
 
posname = paste(directory,"/",currpos,"-possamplemismatches",outputformat, sep = "")





ggplot(data = mismatchmeltagg, aes(x=Sample, y=percentmismatch)) + geom_boxplot(aes(fill=Sample), outlier.shape=NA) + theme_bw()+ geom_jitter(aes(fill=Sample), size = dotsize) +  ggtitle(paste("Position ",currpos," Mismatches", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Sample")+ylab("Percent Misincorporation")
ggsave(filename=posname,width=7, height=7)

mismatchmeltaminoagg <- aggregate(mismatchmelt$percentmismatch, by=list(Sample = mismatchmelt$Sample, amino = trnatable[match(mismatchmelt$Feature,trnatable[,1]),3]), FUN=mean)

#print(head(mismatchmeltaminoagg))
colnames(mismatchmeltaminoagg) <- c("Sample","amino","percentmismatch")


posname = paste(directory,"/",currpos,"-posaminomismatches",outputformat, sep = "")



ggplot(data = mismatchmeltaminoagg, aes(x=amino, y=percentmismatch)) + geom_boxplot(aes(fill=amino), outlier.shape=NA) + theme_bw()+ geom_jitter(aes(fill=amino), size = dotsize) +  ggtitle(paste("Position ",currpos," Mismatches", sep = ""))+  theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Isotype")+ylab("Percent Misincorporation")
ggsave(filename=posname, width=7, height=7)


mismatchbase = mismatches[mismatches$position == currpos,c("Feature","Sample")]
mismatchbase$attotal = (mismatches[mismatches$position == currpos,"adenines"]) / mismatches[mismatches$position == currpos,"coverage"]
mismatchbase$cgtotal = (mismatches[mismatches$position == currpos,"cytosines"] ) / mismatches[mismatches$position == currpos,"coverage"]


posname = paste(currpos,"-posbasedistribution",outputformat, sep = "")


deletemeltaminoagg <- aggregate(deletionmelt$deletepercent, by=list(Sample = deletionmelt$Sample, amino = trnatable[match(deletionmelt$Feature,trnatable[,1]),3]), FUN=mean)

colnames(deletemeltaminoagg) <- c("Feature","amino","percentdeletions")
posname = paste(directory,"/",currpos,"-posaminodeletions",outputformat, sep = "")
ggplot(data = deletemeltaminoagg, aes(x=amino, y=percentdeletions)) + geom_boxplot(aes(fill=amino), outlier.shape=NA) + theme_bw()+ geom_jitter(aes(fill=amino), size = dotsize) +  ggtitle(paste("Position ",currpos," Deletions", sep = ""))+  theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1)   + xlab("Isotype")+ylab("Percent Skipped")
ggsave(filename=posname,width=7, height=7)



deletemeltSampleagg <- aggregate(deletionmelt$deletepercent, by=list(Feature = deletionmelt$Feature, Sample = Sampletable[match(deletionmelt$Sample,Sampletable[,1]),2]), FUN=mean)
colnames(deletemeltSampleagg) <- c("Feature","Sample","percentdeletions")
posname = paste(directory,"/",currpos,"-possampledeletions",outputformat, sep = "")
ggplot(data = deletemeltSampleagg, aes(x=Sample, y=percentdeletions)) + geom_boxplot(aes(fill=Sample), outlier.shape=NA) + theme_bw()+ geom_jitter(aes(fill=Sample), size = dotsize) +  ggtitle(paste("Position ",currpos," Deletions", sep = ""))+  theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Sample")+ylab("Percent Skipped")
ggsave(filename=posname,width=7, height=7)

}

#print("***||")
deletionmelt = mismatches[,c("Feature","position","Sample","coverage","deletions")]
maxskips <- aggregate(deletionmelt$deletions/(deletionmelt$deletions + deletionmelt$coverage + 30), by=list(position = deletionmelt$position, Feature = deletionmelt$Feature,Sample = deletionmelt$Sample), FUN=max)
#maxskips <- aggregate(deletionmelt$deletions/deletionmelt$coverage, by=list(position = deletionmelt$position, Feature = deletionmelt$Feature), FUN=max)

print (maxskips[order(-maxskips$x),])
#head(maxskips)
#head()
