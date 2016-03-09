library(ggplot2)
library(reshape2)
library(scales)
library(getopt)

#args <- commandArgs(trailingOnly = TRUE)

#args <- c("aging-coverage.txt","sacCer3-trnatable.txt", "YeastAging.txt", "aging-coverage.pdf")
#args <- c("ExosomeData-coverage.txt","/soe/holmes/pythonsource/trnatest/hgtrnadb/hg19-trnatable.txt","exosomesamples.txt","ExosomeData-SizeFactors.txt","ExosomeData-coverage.pdf")
#args

#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="Feature.basePosition", value.name="read.density")

spec <- matrix(c(
        'cov'     , 'c', 1, "character", "coverage file from getcoverage.py (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'samples'    , 's', 1, "character", "sample file (required)",
        'allcov'    , 'a', 1, "character", "output coverages for all tRNAs (optional)",
        'multicov'    , 'm', 1, "character", "output coverages for all tRNAs on multiple pages(optional)",
        
        'combinecov'    , 'o', 1, "character", "output coverages for tRNAs combined",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


coverages <- read.table(opt$cov, header = TRUE)
trnatable <- read.table(opt$trna)
sampletable <- read.table(opt$samples)

outputfile <- opt$allcov
combinedfile <- opt$combinecov
multipage <- opt$multicov
#colnames(coverages)

expand.delimited <- function(x, col1=1, col2=2, sep=",") {
  rnum <- 1
  expand_row <- function(y) {
    factr <- y[col1]
    strng <- toString(y[col2])
    expand <- strsplit(strng, sep)[[1]]
    num <- length(expand)
    factor <- rep(factr,num)
    return(as.data.frame(cbind(factor,expand),
          row.names=seq(rnum:(rnum+num)-1)))
    rnum <- (rnum+num)-1
  }
  expanded <- apply(x,1,expand_row)
  df <- do.call("rbind", expanded)
  names(df) <- c(names(x)[col1],names(x)[col2])
  return(df)
}


myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}




#colnames(coverages) <- c("Feature", "Sample",1:(length(colnames(coverages)) - 2))
#trnatable[coveragemelt[,1],c(3,4)]

#remove columns with too many NAs

#length(colSums(is.na(coverages)) < nrow(coverages)/8)
#colSums(is.na(coverages)) < nrow(coverages)/8
#length(!grepl("gap",colnames(coverages), fixed = TRUE))
#colnames(coverages)
#coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8 | !grepl("gap",colnames(coverages), fixed = TRUE) | !grepl("intron",colnames(coverages), fixed = TRUE)]
coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8]



#aggregate(coverages, by=list(trnatable[,3]), FUN=sum)[2]


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


#coveragemelt$Feature %in% strsplit(trnatable[,2], ",", fixed = TRUE)


#vector of amino acids
#acceptorType = trnatable[match(coveragemelt$Feature, trnatable[,1]),3]
locustable = expand.delimited(trnatable,3,2)
acceptorType = locustable[match(coveragemelt$Feature, locustable[,2]),1]
#locustable[1,]
#acceptorType

#coveragemeltagg <- aggregate(coveragemelt$value, by=list(Sample = coveragemelt$Sample, variable = coveragemelt$variable), FUN=mean)
#colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
#coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)

#colnames(coveragemeltagg)
#coveragemelt[1,]
#length(acceptorType)
#length(coveragemelt$Feature)
sortcovmelt <- coveragemelt[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature),]  
sortacceptor <- acceptorType[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature)]

covsummary <- ggplot(sortcovmelt,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor))) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+theme(axis.text.y=element_text(colour="black",size=8), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=8),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8))+ ylab("Read Share") +   scale_y_continuous(breaks=myBreaks, labels = c("0","1")) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))
ggsave(filename=combinedfile,covsummary)


#ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid( ~ Sample, scales="free") + geom_bar(aes(fill = factor(acceptorType)),stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=8),axis.text.x=element_text(colour="black"), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=5))+ ylab("Normalized read count") +   scale_y_continuous(breaks=myBreaks) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype") +scale_x_discrete(breaks=c("X1","X37","X38","X73","tail3"), labels=c("start", "intronstart","intronend" ,"end","tail"))#+scale_x_discrete("Position") 
#ggsave(filename=combinedfile,width=1.5*length(unique(coveragemelt$Sample)), limitsize=FALSE)
#exit()


#sessionInfo()
#head(coveragemelt[coveragemelt$value > 1, ])

locusplot <- ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=4),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0,size=6),strip.text.x = element_text(angle=0,size=8))+ ylab("read count") +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail"))  
ggsave(locusplot, filename=outputfile,height=.5*length(unique(coveragemelt$Feature)),width=2*length(unique(coveragemelt$Sample)), limitsize=FALSE)
#ggsave(filename=outputfile)



#ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.x=element_text(colour="black",size=2), strip.text.y = element_text(angle=0))+ ylab("Normalized read count")+scale_x_discrete("Position")

 #scale_x_discrete()


#theme(axis.title.x=element_blank(), axis.text.x=element_text(colour="black"), strip.text.y = element_text(angle=0))
#+ scale_x_discrete()

