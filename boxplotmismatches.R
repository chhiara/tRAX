library(ggplot2)
library(RColorBrewer)

library(reshape2)
library(ggsci)
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
        'comparisons'     , 'c', 1, "character", "comparisons file",
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


positionorder = c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76')

aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","M")


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



mismatchmelt = mismatches[,c("Feature","Sample","percentmismatch","position")]


features = unique(mismatches$Feature)
#print(features)

#relevantmismatches <- list()

relevantmismatches = data.frame("Samples", "position", "tRNA","first","second","difference")
colnames(relevantmismatches) = c("Samples", "position", "tRNA","first","second","difference")


currlistpos = 1


dotsize = .4
aminodotsize = .8

mismatchmeltfilter = mismatches[	mismatches$adenines+mismatches$thymines+mismatches$cytosines+mismatches$guanines > 50,c("Feature","Sample","percentmismatch","position")]
#head(mismatchmelt$adenines+mismatchmelt$thymines+mismatchmelt$cytosines+mismatchmelt$guanines)

mismatchmeltposagg <- aggregate(mismatchmeltfilter$percentmismatch, by=list(position = mismatchmeltfilter$position, Feature = mismatchmeltfilter$Feature), FUN=max)


mismatchmeltposagg <- mismatchmeltposagg[mismatchmeltposagg$position %in% positionorder,]
mismatchmeltposagg$position = factor(mismatchmeltposagg$position, levels = positionorder)

colnames(mismatchmeltposagg) <- c("Position","Feature","percentmismatch")

mismatchmeltposagg$amino = trnatable[match(mismatchmeltposagg$Feature,trnatable[,1]),3]
mismatchmeltposagg$anticodon = trnatable[match(mismatchmeltposagg$Feature,trnatable[,1]),4]



posname = paste(directory,"/trnapositionmismatches",outputformat, sep = "")
# geom_boxplot(aes(fill=position), outlier.shape=NA) 
ggplot(data = mismatchmeltposagg, aes(x=Position, y=percentmismatch)) + theme_bw()+ geom_jitter(aes(color=amino), size = 1,width = 0.25) +  ggtitle(paste("Position Mismatches", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Position")+ylab("Maximum percent Misincorporation")
ggsave(filename=posname,width=20, height=7)

fiveprimemeltfilter = mismatches[	mismatches$tRNAreadstotal > 30,c("Feature","Sample","position","readstarts","tRNAreadstotal","coverage")]
fiveprimemeltfilter$percentstart = fiveprimemeltfilter$readstarts/fiveprimemeltfilter$tRNAreadstotal

threeprimemeltfilter = mismatches[	mismatches$tRNAreadstotal > 30,c("Feature","Sample","position","readends","tRNAreadstotal","coverage")]
threeprimemeltfilter$percentstart = threeprimemeltfilter$readends/threeprimemeltfilter$tRNAreadstotal

mismatchfilter = mismatches[mismatches$tRNAreadstotal > 30,c("Feature","Sample","position","percentmismatch","coverage")]

#print(length(fiveprimemeltfilter$percentstart))
#print(length(fiveprimemeltfilter$position))
#print(length(fiveprimemeltfilter$Feature))
fiveprimeposagg <- aggregate(fiveprimemeltfilter$percentstart, by=list(position = fiveprimemeltfilter$position, Sample = fiveprimemeltfilter$Sample), FUN=mean)
fiveprimeposagg <- fiveprimeposagg[fiveprimeposagg$position %in% positionorder,]
fiveprimeposagg$position = factor(fiveprimeposagg$position, levels = positionorder)
#print("**")
colnames(fiveprimeposagg) <- c("Position","Sample","percentstart")
print("||**")
fiveprimemeltfilter$amino = trnatable[match(fiveprimemeltfilter$Feature,trnatable[,1]),3]
threeprimemeltfilter$amino = trnatable[match(threeprimemeltfilter$Feature,trnatable[,1]),3]
mismatchfilter$amino = trnatable[match(mismatchfilter$Feature,trnatable[,1]),3]
#fiveprimeposagg$anticodon = trnatable[match(fiveprimeposagg$Feature,trnatable[,1]),4]


#print("**||")
posname = paste(directory,"/trnapositionfiveprime",outputformat, sep = "")
# geom_boxplot(aes(fill=position), outlier.shape=NA) 
ggplot(data = fiveprimeposagg, aes(x=Position, y=percentstart)) + theme_bw()+ geom_jitter(aes(color=Sample), size = dotsize,width = 0.25) +  ggtitle(paste("Position Starts", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Position")+ylab("Maximum Starts")
ggsave(filename=posname,width=20, height=7)
#print("**")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

for (curramino in unique(trnatable[,3])){

posname = paste(directory,"/",curramino,"-trnapositionmismatches",outputformat, sep = "")
# geom_boxplot(aes(fill=position), outlier.shape=NA) 
mismatchmeltposaggamino = mismatchmeltposagg[mismatchmeltposagg$amino == curramino,]
ggplot(data = mismatchmeltposaggamino, aes(x=Position, y=percentmismatch)) + theme_bw()+ geom_jitter(aes(color=anticodon), size = aminodotsize,width = 0.25) +  ggtitle(paste("Position Mismatches", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Position")+ylab("Maximum percent Misincorporation")
ggsave(filename=posname,width=20, height=7)


#positionagg <- aggregate(fiveprimemeltfilter$percentstart, by=list(Feature = fiveprimemeltfilter$Feature, Sample = Sampletable[match(fiveprimemeltfilter$Sample,Sampletable[,1]),2],Position=fiveprimemeltfilter$position ), FUN=mean)
#print("**")
#fiveprimeaminoposagg = dcast(fiveprimeposagg,Feature+Position~Sample,mean)
#fiveprimeaminoposagg = aggregate(fiveprimemeltfilter$percentstart, by=list(Feature = fiveprimemeltfilter$, Sample = Sampletable[match(fiveprimemeltfilter$Sample,Sampletable[,1]),2],Position=fiveprimemeltfilter$position ), FUN=mean)

fiveprimeamino = fiveprimemeltfilter[fiveprimemeltfilter$amino == curramino,]
threeprimeamino = threeprimemeltfilter[threeprimemeltfilter$amino == curramino,]
mismatchamino = mismatchfilter[mismatchfilter$amino == curramino,]

fiveprimeaminoposagg = aggregate(fiveprimeamino$percentstart, by=list( Sample = Sampletable[match(fiveprimeamino$Sample,Sampletable[,1]),2],Position=fiveprimeamino$position ), FUN=mean) 

fiveprimeaminoposagg <- fiveprimeaminoposagg[fiveprimeaminoposagg$Position %in% positionorder,]
fiveprimeaminoposagg$Position = factor(fiveprimeaminoposagg$Position, levels = positionorder)

colourCount = length(unique(fiveprimeaminoposagg$Position))+1

#print(head(fiveprimeaminoposagg))
currplot = ggplot(fiveprimeaminoposagg,aes(x = Sample, y = x,fill = Position, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "fill",stat="identity") +
    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Percentage of reads that start at position") + 
    labs(fill="position")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    scale_fill_manual(values = getPalette(colourCount))+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")


ggsave(paste(directory,"/",curramino ,"-fiveprimecounts",outputformat,sep= ""), currplot)
#print(head(fiveprimeamino))
#print(max(fiveprimeamino$percentstart))

fiveprimeaminoposagg = aggregate(fiveprimeamino$percentstart, by=list( Sample = Sampletable[match(fiveprimeamino$Sample,Sampletable[,1]),2],position=fiveprimeamino$position, Feature = fiveprimeamino$Feature ), FUN=mean) 

fiveprimeaminoposagg$percentstart = fiveprimeaminoposagg$x
fiveprimeaminoposagg <- fiveprimeaminoposagg[fiveprimeaminoposagg$position %in% positionorder,]
fiveprimeaminoposagg$position = factor(fiveprimeaminoposagg$position, levels = positionorder)
currplot = ggplot(fiveprimeaminoposagg,aes(x = position,y= Sample, fill = percentstart)) + geom_tile()+scale_fill_gradient(low="white", high="blue", limits = c(0,1))+ theme_bw() +facet_grid(rows = vars(Feature))+ theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
    theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=24),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=16), strip.text.y = element_text(size=24))+
    xlab("Sample") +
    ylab("Percentage of reads that start at position") + 
    labs(fill="position")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=16,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")


ggsave(paste(directory,"/",curramino ,"-fiveprimeheatmap",outputformat,sep= ""), currplot,height=.25*length(unique(fiveprimeamino$Feature))*length(unique(fiveprimeamino$Sample)),width=.25*length(unique(fiveprimeamino$position)), limitsize=FALSE)




threeprimeaminoposagg = aggregate(threeprimeamino$percentstart, by=list( Sample = Sampletable[match(threeprimeamino$Sample,Sampletable[,1]),2],position=threeprimeamino$position, Feature = threeprimeamino$Feature ), FUN=mean)
#print(head(threeprimeaminoposagg))
threeprimeaminoposagg$percentstart = threeprimeaminoposagg$x
threeprimeaminoposagg <- threeprimeaminoposagg[threeprimeaminoposagg$position %in% positionorder,]
threeprimeaminoposagg$position = factor(threeprimeaminoposagg$position, levels = positionorder)

currplot = ggplot(threeprimeaminoposagg,aes(x = position,y= Sample, fill = percentstart)) + geom_tile()+scale_fill_gradient(low="white", high="blue", limits = c(0,1))+ theme_bw() +facet_grid(rows = vars(Feature))+ theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
    theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=24),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=16), strip.text.y = element_text(size=24))+
    xlab("Sample") +
    ylab("Percentage of reads that start at position") + 
    labs(fill="position")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=16,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")


ggsave(paste(directory,"/",curramino ,"-threeprimeheatmap",outputformat,sep= ""), currplot,height=.25*length(unique(fiveprimeamino$Feature))*length(unique(fiveprimeamino$Sample)),width=.25*length(unique(fiveprimeamino$position)), limitsize=FALSE)

mismatchaminoposagg = aggregate(mismatchamino$percentmismatch, by=list( Sample = Sampletable[match(mismatchamino$Sample,Sampletable[,1]),2],position=mismatchamino$position, Feature = mismatchamino$Feature ), FUN=mean)
#print("**::")
#print(head(mismatchaminoposagg))
mismatchaminoposagg$percentmismatch = mismatchaminoposagg$x
#print("::**")
mismatchaminoposagg <- mismatchaminoposagg[mismatchaminoposagg$position %in% positionorder,]
mismatchaminoposagg$position = factor(mismatchaminoposagg$position, levels = positionorder)

currplot = ggplot(mismatchaminoposagg,aes(x = position,y= Sample, fill = percentmismatch)) + geom_tile()+scale_fill_gradient(low="white", high="blue", limits = c(0,1))+ theme_bw() +facet_grid(rows = vars(Feature))+ theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
    theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=24),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=16), strip.text.y = element_text(size=24))+
    xlab("Sample") +
    ylab("Percentage of reads that start at position") + 
    labs(fill="position")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=16,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")


ggsave(paste(directory,"/",curramino ,"-mismatchheatmap",outputformat,sep= ""), currplot,height=.25*length(unique(fiveprimeamino$Feature))*length(unique(fiveprimeamino$Sample)),width=.25*length(unique(fiveprimeamino$position)), limitsize=FALSE)


}


for (currpos in positions){
#print(currpos)

mismatchmelt = mismatches[mismatches$position == currpos,c("Feature","Sample","percentmismatch")]
deletionmelt = mismatches[mismatches$position == currpos,c("Feature","Sample","deletions","coverage")]
fiveprimemelt = mismatches[mismatches$position == currpos,c("Feature","Sample","readstarts","tRNAreadstotal","coverage")]
fiveprimemelt$percentstart = fiveprimemelt$readstarts/fiveprimemelt$tRNAreadstotal

deletionmelt$deletepercent = deletionmelt$deletions/(deletionmelt$deletions + deletionmelt$coverage + 30)


fiveprimeagg <- aggregate(fiveprimemelt$percentstart, by=list(Feature = fiveprimemelt$Feature, Sample = Sampletable[match(fiveprimemelt$Sample,Sampletable[,1]),2]), FUN=mean)
#print(head(fiveprimeagg))
colnames(fiveprimeagg) <- c("Feature","Sample","percentstart")

posname = paste(directory,"/",currpos,"-possamplereadstarts",outputformat, sep = "")
ggplot(data = fiveprimeagg, aes(x=Sample, y=percentstart)) + geom_boxplot(aes(fill=Sample), outlier.shape=NA) + theme_bw()+ geom_jitter(aes(fill=Sample), size = dotsize) +  ggtitle(paste("Position ",currpos,"Read Starts", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+  ylim(0, 1) + xlab("Sample")+ylab("Percent Read Starts")
ggsave(filename=posname,width=7, height=7)


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

#print (maxskips[order(-maxskips$x),])
#head(maxskips)
#head()
#print("|||**||")
#print(head(fiveprimemeltfilter))
#print(length(Sampletable[match(fiveprimemeltfilter$Sample,Sampletable[,1]),2]))
#print(length(fiveprimemeltfilter$position))


