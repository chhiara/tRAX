library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

#args[3]
#args[4]


log2_minor_break = function (...){
  function(x) {
    minx         = floor(min(log2(x), na.rm=T))-1;
    maxx         = ceiling(max(log2(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log2(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(2^(minor_breaks))
  }
}

outputformat <- ".pdf"
#uniqueInitials <- c("a", "A", "c", "C", "j", "J", "R")


trnasize = 6

experimentname <- args[1]
counts <- read.table(args[2], stringsAsFactors = FALSE)
trnatable <- read.table(args[3], stringsAsFactors = FALSE)

types <- read.table(args[4])

sampletable <- read.table(args[5], stringsAsFactors = FALSE)

comparisons <- read.table(args[6], stringsAsFactors = FALSE)

#psuedocounts for log-transformation
counts = counts + 1

# Rscript ~/pythonsource/trnaseq/makescatter.R compcellexo-normalized.txt SRR2084358 SRR2084360

#comparisons = strsplit(sampletable[5:length(args)], ",", fixed = TRUE)
#comparisons


#lincrnas
counts <- merge(types,counts, by.x = 1, by.y = 0)
colnames(counts)[2] <- "type"


#othertypes = !(counts[,"type"] == "snoRNA" | counts[,"type"] == "miRNA" | counts[,"type"]  == "trna_fiveprime" | counts[,"type"] == "trna_threeprime")
genetypes = c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA","snoRNA","miRNA")
othertypes = !(counts[,"type"] %in% genetypes)
#length(othertypes)
#length(counts[,1])
#levels(counts[othertypes,"type"])  <- c(levels(counts[othertypes,"type"]), "other") 




counts[,"type"] <- as.character(counts[,"type"])
counts[othertypes,"type"] <- "other"
counts[,"type"] <- as.factor(counts[,"type"])
#head(counts[grepl("tRNA", counts[,1]),])



counts[,"type"] <- factor(counts[,"type"], levels = c(genetypes,"other"))
trnatypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA")
fragtypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts")
trnagenes =  counts[,"type"] %in% trnatypes

counts <- rbind(counts[counts[,"type"] == "other",],counts[counts[,"type"] == "snoRNA",],counts[counts[,"type"] == "miRNA",], counts[trnagenes,])
#head(unlist(strsplit(as.character(counts[,1]), split='_', fixed=TRUE))[1])
#head(strsplit(counts[,1], split='_', fixed=TRUE))
#head(counts[!is.character(counts[,1]),])

#counts$amino <- trnatable[match(trnatable[1,], counts[1,]),2]
#trnatable[match(trnatable[,1], counts[,1]),3]
#na.omit(match(trnatable[,1], counts[,1])),3]
#got to find a way to do this with the fragment names
#head(trnatable)
colnames(counts)[1] <- "name"
colnames(trnatable) <- c("trnaname", "loci", "amino", "anticodon")
#as.character(counts[,"name"])
#unlist(strsplit(as.character(counts[,"name"]), split='_', fixed=TRUE))

#counts[,"type"] %in% trnatypes
counts$trnaname = ifelse( counts[,"type"] %in%  fragtypes,sub("^(.*)_[^_]*$", "\\1", as.character(counts[,"name"])),as.character(counts[,"name"]))
#head(counts[counts[,"type"] %in% trnatypes,])

trnacounts <- merge(trnatable,counts, by.x = "trnaname", by.y = "trnaname", all.y=TRUE)

#head(trnacounts[ trnacounts[,"type"] %in% trnatypes,])


#head(counts[counts$amino == "Ala",])
#comparisons
#sampletable
#comparisons[1,1]
i = 1
#trnacounts[is.na(trnacounts$amino),"amino"] = "non-tRNA"
trnacounts <- rbind(trnacounts[trnacounts[,"type"] == "other",],trnacounts[trnacounts[,"type"] == "snoRNA",],trnacounts[trnacounts[,"type"] == "miRNA",], trnacounts[trnacounts[,"type"] %in% trnatypes,])
trnacounts$dotsize = ifelse(trnacounts[,"type"] %in% trnatypes, trnasize, 1)
#tail(trnacounts)
#xaxis = sampletable[sampletable[,2] == comparisons[i,1],1][1]
#yaxis = sampletable[sampletable[,2] == comparisons[i,2],1][2]

#sampletable[,2]
#comparisons[i,2]

#xaxis
#yaxis

dashinterc = 1.5
#cor.test(log(counts[,xaxis]+1),log(counts[,yaxis]+1))

onlytrnas <- trnacounts[trnacounts[,"type"] %in% trnatypes,]

#trnacounts[,"type"]
#
trnacounts$fragtype <- factor(ifelse(as.character(trnacounts[,"type"]) %in% trnatypes, as.character(trnacounts[,"type"]), "nontRNA" ), levels = c("nontRNA",trnatypes))
#revalue(trnacounts$fragtype, c("beta"="two", "gamma"="three"))
#head(unique(trnacounts$fragtype))
#scale_shape_manual(values = c(20, 40, 41, ,21))

nlevels(trnacounts$amino)
trnacounts$amino <- as.factor(trnacounts$amino)
#levels(trnacounts$amino) <- sort(trnacounts$amino)


aminoinfo <- read.table("/projects/lowelab/users/holmes/pythonsource/trnaseq/aminotable.txt", stringsAsFactors = FALSE)
#aminoinfo <- aminoinfo[,3]
#aminoinfo[]

aminoletters <- aminoinfo[match(aminoinfo[,2], levels(trnacounts$amino)),3]
aminoletters <- aminoinfo[match(levels(trnacounts$amino),aminoinfo[,2]),3]

aminoletters[is.na(aminoletters)] <- "X"
#print(aminoletters)
#levels(trnacounts$amino)


aminoletters <-  unlist(lapply(aminoletters, utf8ToInt))
#match(aminoinfo[,2],trnacounts$amino)
#print(trnacounts$fragtype)

trnacounts$fragtype <- mapvalues(trnacounts$fragtype, from = c("trna_threeprime", "trna_fiveprime", "trna_other", "trna_wholecounts"), to = c("Three-prime fragments","Five-prime fragments","Other fragments","Whole tRNAs"))


for (i in 1:length(rownames(comparisons))){

yaxis = sampletable[sampletable[,2] == comparisons[i,1],1][1]
xaxis = sampletable[sampletable[,2] == comparisons[i,2],1][1]

yname = comparisons[i,1]
xname = comparisons[i,2]


if(comparisons[i,2] == comparisons[i,1])
{
xaxis = sampletable[sampletable[,2] == comparisons[i,2],1][2]

}


#print(xaxis)
#print(yaxis)

#maxlim = max(max(log(counts[,xaxis])), max(log(counts[,yaxis]))) #, xlim = c(0,maxlim), ylim = c(0,maxlim)
maxlim = max(c(max(counts[,xaxis]), max(counts[,yaxis]))) #, xlim = c(0,maxlim), ylim = c(0,maxlim) + 1


corr = cor.test(log(counts[,xaxis]+1),log(counts[,yaxis]+1))

sublabel = paste("Pearson Correlation: ",corr$estimate, sep = "")  



currplot <- qplot(data=counts,x=counts[,xaxis],y=counts[,yaxis],xlab = xaxis,ylab = yaxis,color=type, asp=1)+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
currplot <- arrangeGrob(currplot, sub = textGrob(sublabel, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontsize = 14)))
ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-typescatter",outputformat,sep= ""), currplot)


#+scale_shape_manual(values = c(20,15, 17,18,19)) theme(legend.position = "bottom")
#currplot <- ggplot(trnacounts, aes_string(x=xaxis, y=yaxis))+geom_point(aes(color=amino, size = dotsize, shape=fragtype))+guides(size=FALSE, ncol = 1)+scale_size_continuous(range = c(.75,2))+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#scale_shape_manual(values = aminoletters)
#scale_shape_manual(values=1:nlevels(trnacounts$amino))

currplot <- ggplot(trnacounts[trnacounts$fragtype != "nontRNA",], aes_string(x=xaxis, y=yaxis))+geom_point(data = transform(trnacounts[trnacounts$fragtype == "nontRNA",], fragtype=NULL), aes(size = dotsize))+xlab(xname)+ylab(yname)+geom_point(aes(shape=amino, color=amino, size = dotsize))+guides(size=FALSE, ncol = 1)+scale_shape_manual(values=aminoletters) +facet_wrap( ~fragtype , ncol = 2) + scale_size_continuous(range = c(.75,4))+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim),breaks = trans_breaks('log2', function(x) 2^x, n = 10)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim),breaks= trans_breaks('log2', function(x) 2^x, n = 10)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-aminoscatter",outputformat,sep= ""), currplot, height = 10, width = 12)

#currplot <- qplot(data=onlytrnas,x=onlytrnas[,xaxis],y=onlytrnas[,yaxis],xlab = xaxis,ylab = yaxis, asp=1)+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) +facet_wrap( ~amino, , ncol = 2)
#ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-tRNAscatter.pdf",sep= ""), currplot,height=4*length(unique(trnacounts$amino)),width=20, limitsize=FALSE)




}