library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)

args <- commandArgs(trailingOnly = TRUE)



#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")


spec <- matrix(c(
        'cov'     , 'c', 1, "character", "coverage file from getcoverage.py (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'samples'    , 's', 1, "character", "sample file (required)",
        'uniquename'    , 'u', 1, "character", "file header for uniqifying files",
        'directory'    , 'f', 1, "character", "directory to place amino acid files",
        'allcov'    , 'a', 1, "character", "output coverages for all tRNAs (optional)",
        'multicov'    , 'm', 1, "character", "output coverages for all tRNAs on multiple pages(optional)",
        
        'combinecov'    , 'o', 1, "character", "output coverages for tRNAs combined",
        'modomics'    , 'd', 2, "character", "optional modomics file",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


coverages <- read.table(opt$cov, header = TRUE)
trnatable <- read.table(opt$trna)
sampletable <- read.table(opt$samples)
modomicstable <- data.frame()



outputformat <- ".pdf"





#modomicstable <- read.table("/data/scratch/importrnaseq/agingtest/sacCer3-modomics.txt", header = TRUE)

outputfile <- opt$allcov
combinedfile <- opt$combinecov
multipage <- opt$multicov
#colnames(coverages)

myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}
configurecov <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.5, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.5, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=6)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=6))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=3))

              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .1))
              
              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))
covplot<- covplot+theme(panel.margin = unit(.15, "lines"))
#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))

covplot
}



#colnames(coverages) <- c("Feature", "Sample",1:(length(colnames(coverages)) - 2))
#trnatable[coveragemelt[,1],c(3,4)]


#remove columns with too many NAs
coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8]
#coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8 | !grepl("gap",colnames(coverages), fixed = TRUE) | !grepl("intron",colnames(coverages), fixed = TRUE)]
#colnames(coverages)
#aggregate(coverages, by=list(trnatable[,3]), FUN=sum)[2]

#1:(length(colnames(coverages)) - 2)
#coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"), measure.vars = as.character(colnames(coverages)[3:length(colnames(coverages)) - 2]))
coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"))


#trnatable[coveragemelt[,1],c(3,4)]
#colnames(coveragemelt)
#aggregate(coveragemelt, by=list(trnatable[,3]), FUN=mean)[2]
#ddply(coveragemelt, ,summarise,value = sum(!is.na(value)))

#coveragemelt$value
coveragemelt[is.na(coveragemelt)] <- 0


#ddply(coveragemelt, c("Feature", "sample", "variable"),summarise,value = sum(!is.na(value)))

#aggregate(coveragemelt, data=dat, FUN = function(x) c(M=mean(x), SD=sd(x)))

#tapply(coveragemelt$value,sampletable[coveragemelt$Sample,2])
#( normalizationtable[2,coveragemelt$Sample])

#Normalization
#This normalization takes way too long
#coveragemelt$value <- as.vector(coveragemelt$value / as.vector( normalizationtable[1,coveragemelt$Sample]), mode = "numeric")
#coveragetest <- coveragemelt

#coveragemelt
#out <- as.vector(coveragetest$value / as.vector( normalizationtable[1,coveragetest$Sample]), mode = "numeric")


coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = sampletable[match(coveragemelt$Sample,sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)
colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
coveragemelt <- coveragemeltagg
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2" & coveragemelt$Sample == "dmSCDd12",]
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2" & coveragemelt$Sample == "dmSCDd1",]
#coveragemelt[coveragemelt$tRNA == "tRNA-Phe-GAA-2",]
#"41390"	"tRNA-Leu-TAA-1"	"dmMet_Amino"	"105"	546.017749321

write.table(coveragemelt,"aggtables.txt" ,sep = "\t")



#acceptorType

#coveragemeltagg <- aggregate(coveragemelt$value, by=list(Sample = coveragemelt$Sample, variable = coveragemelt$variable), FUN=mean)
#colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
#coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)

#colnames(coveragemeltagg)
#unique(coveragemelt$tRNA) 

#acceptorType
acceptorType = trnatable[match(coveragemelt$Feature, trnatable[,1]),3]
acceptorType <- factor(acceptorType, levels = sort(unique(acceptorType)))


#acceptorType <- factor(acceptorType)

#coveragemelt = coveragemelt[order(acceptorType),]

#unique(coveragemelt$variable)
#unique(acceptorType)


sortcovmelt <- coveragemelt[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature),]  
sortacceptor <- acceptorType[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature)]


covsummary <- ggplot(coveragemelt,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor))) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+theme(axis.text.y=element_text(colour="black",size=8), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=8),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8))+ ylab("Normalized Read Count") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))
ggsave(filename=combinedfile,covsummary)


smallcovsummary <- ggplot(coveragemelt,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor)),width = 2, size = 2	) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+ ylab("Normalized Read Count") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail"),expand = c(0.05, .01)) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))
																																#																																	#, expand = c(0.05, 0)
smallcovsummary<-smallcovsummary+ theme_bw()
smallcovsummary<-smallcovsummary+theme(text = element_text(size = 4))
smallcovsummary<-smallcovsummary+theme(rect = element_rect(size=.01))

smallcovsummary<-smallcovsummary+theme(legend.key.size =  unit(0.01, "in")) # Change key size in the legend
smallcovsummary<-smallcovsummary+theme(legend.key.height =  unit(0.01, "in")) # Change key size in the legend
 
smallcovsummary<-smallcovsummary+theme(legend.text = element_text(size=4)) # Change the size labels in the legend
smallcovsummary<-smallcovsummary+theme(legend.title= element_text(size=4)) # Change key title size in the legend

smallcovsummary<-smallcovsummary+theme(axis.text.y=element_text(size=3))


smallcovsummary<-smallcovsummary+theme(axis.ticks = element_line(size = .1)) 
#smallcovsummary<-smallcovsummary+theme(axis.ticks.y=element_blank())

smallcovsummary<-smallcovsummary+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=3))

smallcovsummary<-smallcovsummary+theme(strip.text.y = element_text(angle=0,size=8))
smallcovsummary<-smallcovsummary+theme(strip.text.x = element_text(angle=0,size=8))

smallcovsummary<-smallcovsummary+theme(axis.line = element_line(size = .1))
smallcovsummary<-smallcovsummary+theme(panel.grid  = element_line(size = .1))

#smallcovsummary<-smallcovsummary+theme(strip.switch.pad.grid = unit(.1, "cm"))


#smallcovsummary<-smallcovsummary+theme(strip.background  = element_rect(size=.01))
scalefactor = .5
ggsave(filename=paste(opt$unique, "-smallcombineduniquecoverages.pdf",sep= ""), smallcovsummary, width =scalefactor*(.5+.75*length(unique(coveragemelt$Sample))), height = scalefactor*1)
ggsave(filename=paste(opt$unique, "-smallcombineduniquecoverages.png",sep= ""), smallcovsummary, width =.5+.75*length(unique(coveragemelt$Sample)), height = 1)
ggsave(filename=paste(opt$unique, "-smallcombineduniquecoverages.eps",sep= ""), smallcovsummary, width =.5+.75*length(unique(coveragemelt$Sample)), height = 1)


#axis.text.y=element_text is y-axis coverage labels
#strip.text.y is trna names
#strip.text.x is experiment names

#ggsave(filename="testwrap.pdf", allcoverages,scale=2)
#dev.off()
#pdf(outputfile, height=1*length(unique(coveragemelt$Feature)),width=5*length(unique(coveragemelt$Sample)))

#reformat <- function(x,lab="\n"){ sapply(x, function(c){ paste(unlist(strsplit(as.character(c) , split="_")),collapse=lab) }) }
    
#coveragemelt$Sample <- factor(coveragemelt$Sample, labels=unique(reformat(coveragemelt$Sample))

modomicstable <- data.frame(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)

modomicstable <- read.table(text ="",col.names = c("trna", "mod", "pos"),colClasses = c("character", "character", "character")) #(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)
#colNames(modomicstable) <- c("trna", "mod","pos")
if(!is.null(opt$modomics) && file.exists(opt$modomics) ){
    modomicstable <- read.table(opt$modomics, header = TRUE)
    modomicstable$pos <- paste("X",modomicstable$pos, sep = "")
}



modomicstable <- modomicstable[as.character(modomicstable$pos) %in% unique(as.character(coveragemelt$variable)),]
modomicstable$dist <- match(modomicstable$pos,  levels(coveragemelt$variable)) #factor(modomicstable$pos, levels = levels(coveragemelt$variable))
modomicstable$Feature <- factor(modomicstable$trna,levels = levels(coveragemelt$Feature) )

stopmods = c("m1A","m2,2G","m1G","m1I","m3C")
modomicstable <- modomicstable[modomicstable$mod %in% stopmods,]

colnames(modomicstable)[colnames(modomicstable) == 'mod'] <- 'Modification'


    
    
#modomicstable$dist <- modomicstable$pos
#head(coveragemelt)
#unique(as.character(coveragemelt$variable))
#
#geom_bar(aes(fill = factor(baseMod, levels=c("m1A", "m1G", "m3C", "other bases", "not documented"))), stat="identity")
#modomicstable <- data.frame(Feature = factor(c("nmt-tRNA-Leu-TAA-1"),levels(coveragemelt$Feature)), pos =10)
allcoverages <- ggplot(coveragemelt,aes(x=variable,y=value), size = 2) + theme_bw()+ facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE)+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Normalized Read Count") + xlab("tRNA position") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail"))   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))  
#ggsave(filename=outputfile, width = 30, height = 30)
#3*length(unique(coveragemelt$Feature))
#5*length(unique(coveragemelt$Sample))
#dev.off()
scalefactor = .5

ggsave(filename=outputfile, allcoverages,height=scalefactor*1*length(unique(coveragemelt$Feature)),width=2+scalefactor*5*length(unique(coveragemelt$Sample)), limitsize=FALSE, dpi = 600)
#ggsave(filename=outputfile, allcoverages,dpi = 1000, limitsize=FALSE)
#set dpi

#unique(acceptorType)
if(!is.null(opt$directory) && FALSE){
for (curramino in unique(acceptorType)){
#

aminodata = coveragemelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]
aminocoverage <- ggplot(aminodata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE)+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail")) + guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))

aminoname = paste(opt$directory,"/",curramino,"-coverage",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=2 + scalefactor*1.5*length(unique(aminodata$Feature)),width=scalefactor*5*length(unique(aminodata$Sample)), limitsize=FALSE, dpi = 600)

}
}


if(!is.null(opt$unique)){


multamino <- read.table(paste(opt$unique, "-multaminocoverages.txt",sep= ""), header = TRUE)
multactable <- read.table(paste(opt$unique, "-multaccoverages.txt",sep= ""), header = TRUE)
multtrnas <- read.table(paste(opt$unique, "-multtrnacoverages.txt",sep= ""), header = TRUE)
uniquetable <- read.table(paste(opt$unique, "-uniquecoverages.txt",sep= ""), header = TRUE)

uniquetable$maptype = "Transcript specific"
multactable$maptype = "Isotype Specific"
multamino$maptype = "Nonspecific"
multtrnas$maptype = "Isodecoder Specific" 
#allmulttables <- rbind(uniquetable,multactable,multamino,multtrnas)

#allmulttables <- rbind(multamino,multtrnas,multactable,uniquetable)
allmulttables <- rbind(uniquetable, multtrnas,multactable,multamino)


#print(head(allmulttables))
allmulttables <- allmulttables[ , colSums(is.na(allmulttables)) < nrow(allmulttables)/8]
#allmulttables <- allmulttables[ , colSums(is.na(allmulttables)) < nrow(allmulttables)/8]
#print(colnames(allmulttables))
allmultmelt = melt(allmulttables, id.vars = c("Feature", "Sample", "maptype"))
allmultmelt[is.na(allmultmelt)] <- 0


allmultmeltagg <- aggregate(allmultmelt$value, by=list(Feature = allmultmelt$Feature, Sample = sampletable[match(allmultmelt$Sample,sampletable[,1]),2], maptype = allmultmelt$maptype,variable = allmultmelt$variable), FUN=mean)
colnames(allmultmeltagg)[colnames(allmultmeltagg) == "x"]  <- "value"
allmultmeltagg$Sample <- factor(allmultmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
allmultmelt <- allmultmeltagg

#
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Unique tRNA","Unique Decoder","Unique Acceptor","Nonunique Acceptor"))
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Nonunique Acceptor","Unique Acceptor","Unique Decoder","Unique tRNA"))
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Transcript specific","Isodecoder Specific","Isotype Specific","Nonspecific"))
allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Nonspecific","Isotype Specific","Isodecoder Specific","Transcript specific"))


allmultmelt <- allmultmelt[order(-as.numeric(allmultmelt$maptype)),]
#print(head(allmultmelt))
#allcoverages <- ggplot(coveragemelt,aes(x=variable,y=value), size = 2) +                                           theme_bw()+ facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE)+                     theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+                                            ylab("Normalized Read Count") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail"))   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))  



for (curramino in unique(acceptorType)){

aminodata = allmultmelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]
#This bit messes up the sample ordering
#aminodata$Sample <- gsub("_", " ", aminodata$Sample)  


#scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","acc-loop start","acc-loop end","D-loop start","D-loop end","tail"))
#scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail"))
aminocoverage   <- ggplot(aminodata,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","acc-loop start","acc-loop end","D-loop start","D-loop end","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))
#all formula
#aminocoverage <-ggplot(aminodata,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = 8,angle=0)) + ylab("Normalized Read Count") +   xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))

aminocoverage <- configurecov(aminocoverage)

aminoname = paste(opt$directory,"/",curramino,"-uniqcoverage",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=scalefactor*(2 + 1.5*length(unique(aminodata$Feature))),width=scalefactor*(2+5*length(unique(aminodata$Sample))), limitsize=FALSE, dpi = 600)


}



uniqcoverage <- ggplot(allmultmelt,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = 8,angle=0)) + ylab("Normalized Read Count") +   xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))


#ggsave(paste(opt$unique, "-uniquecoverages",outputformat,sep= ""), width = 15, height = 30)

print(unique(allmultmelt$Sample))
print(scalefactor*5*length(unique(allmultmelt$Sample)))

smallcovall<-uniqcoverage
configurecov <- configurecov(smallcovall)
ggsave(filename=paste(opt$unique, "-uniquecoverages.pdf",sep= ""), uniqcoverage,height=scalefactor*(2 + 1.5*length(unique(allmultmelt$Feature))),width=scalefactor*(2+5*length(unique(allmultmelt$Sample))), limitsize=FALSE, dpi = 600)



#configurecov <- configurecov(smallcovall)

#ggsave(filename=paste(opt$unique, "-smalluniquecoverages.pdf",sep= ""), smallcovall,height=2 + scalefactor*1.5*length(unique(allmultmelt$Feature)),width=2+scalefactor*5*length(unique(allmultmelt$Sample)), limitsize=FALSE, dpi = 600)

#width should be .75, same as for combined
#ggsave(filename=paste(opt$unique, "-smallalluniquecoverages.png",sep= ""), smallcovall, width =2+.75*length(unique(coveragemelt$Sample)), height = 1 + scalefactor*.75*length(unique(allmultmelt$Feature)), limitsize=FALSE)


#ggsave(filename=paste(opt$unique, "-uniquecoverages.svg",sep= ""), uniqcoverage,height=2 + scalefactor*1.5*length(unique(allmultmelt$Feature)),width=scalefactor*5*length(unique(allmultmelt$Sample)), limitsize=FALSE, dpi = 600)

#ggsave(filename=paste(opt$unique, "-uniquecoverages.png",sep= ""), uniqcoverage,height=2 + scalefactor*1.5*length(unique(allmultmelt$Feature)),width=scalefactor*5*length(unique(allmultmelt$Sample)), limitsize=FALSE, dpi = 600)


}





