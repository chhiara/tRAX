


args = commandArgs(trailingOnly=TRUE)
if(args[1] == "DESeq2") {library("DESeq2")}
if(args[1] == "DESeq") {library("DESeq")}
library("ggplot2")
library("reshape2")
library("scales")
library("plyr")
library("gridExtra")
library("getopt")
library("ggrepel")
#library("rbamtools") #I use this to test module checker, it's not needed for the pipeline
