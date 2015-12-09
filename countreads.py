#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools

dbname = 'sacCer3'
'''
~/pythonsource/trnaseq/countreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/KitComparison.txt --bedfile=hg19-transcripts.bed --maturetrnas=hg19-maturetRNAs.bed --trnaloci=hg19-trnas.bed

~/pythonsource/trnaseq/countreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bedfile=hg19-transcripts.bed --maturetrnas=hg19-maturetRNAs.bed --trnaloci=hg19-trnas.bed >allcounts.txt

Rscript ~/pythonsource/trnaseq/analyzecounts.R YeastAging featurecounts.txt agingshort.txt dmStat_Amino:dmAll_Amino dmStat_Amino:dmMet_Amino dmStat_Amino:dmLeu_Amino


This currently swaps when given huge bed files for features
'''

def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)
        

count = 0



parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='bed file with non-tRNA features')
parser.add_argument('--gtffile',  nargs='+', default=list(),
                   help='gtf file with non-tRNA features')
parser.add_argument('--ensemblgtf',
                   help='ensembl gtf file with tRNA features')
parser.add_argument('--gtftrnas',  nargs='+', default=list(),
                   help='gtf file with tRNA features')
parser.add_argument('--trnaloci',  nargs='+', default=list(),
                   help='bed file with tRNA features')
parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                   help='bed file with mature tRNA features')
parser.add_argument('--onlyfullpretrnas', action="store_true", default=False,
                   help='only include full pre-trnas')

parser.add_argument('--removepseudo', action="store_true", default=False,
                   help='remove psuedogenes from ensembl GTFs')
parser.add_argument('--genetypefile',
                   help='Output file with gene types')
parser.add_argument('--trnacounts',
                   help='Output file with just trna gene counts')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='disable fragment determination')
parser.add_argument('--nomultimap', action="store_true", default=False,
                   help='do not count multiply mapped reads')


args = parser.parse_args()
includebase = args.nofrag

fullpretrnasonly = args.onlyfullpretrnas
removepseudo = args.removepseudo

ensemblgtf = args.ensemblgtf

nomultimap = args.nomultimap

wholetrnas = dict()
fivefrags = dict()
threefrags = dict()
trailerfrags = dict()
otherfrags = dict()
allfrags = dict()


alltrnas = list()


samplefiles = dict()

sampledata = samplefile(args.samplefile)
samples = sampledata.getsamples()

fullpretrnathreshold = 2

#print >>sys.stderr, " ".join(bamlist)
try:
    featurelist = list()
    trnaloci = list()
    for currfile in args.bedfile:
        featurelist.extend(list(readfeatures(currfile, removepseudo = removepseudo)))
    trnalist = list()
    for currfile in args.trnaloci:
        trnaloci.extend(list(readbed(currfile)))
    for currfile in args.maturetrnas:
        trnalist.extend(list(readbed(currfile)))       
    if ensemblgtf is not None:    
        embllist = list(readgtf(ensemblgtf, filterpsuedo = removepseudo))
    else:
        embllist = list()
except IOError as e:
    print >>sys.stderr, e
    sys.exit()
#if checkduplicates(curr.name for curr in trnalist+trnaloci+featurelist):
    #sys.exit()
#featurelist = list(readbed(sys.argv[1]))
#featurelist = list(readbed(sys.argv[1]))

#./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
'''
./countcomplete.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
'''

featcount = defaultdict(int)
allfeats = featurelist+trnaloci+trnalist

if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
    #print >>sys.stderr, list(curr.name for curr in featurelist )
    #print >>sys.stderr, len(set(curr.name for curr in featurelist))
    print >>sys.stderr, "Duplicate names in feature list:"
    #print >>sys.stderr, ",".join(getdupes(curr.name for curr in allfeats))
    #currname
    #sys.exit(1)


#featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
alltrnas = list(curr.name for curr in featurelist)
#print >>sys.stderr, "***"

#lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
counts = defaultdict(lambda: defaultdict(int))
trnacounts = defaultdict(lambda: defaultdict(int))
trnawholecounts = defaultdict(lambda: defaultdict(int))
trnafivecounts = defaultdict(lambda: defaultdict(int))
trnathreecounts = defaultdict(lambda: defaultdict(int))

trnalocuscounts = defaultdict(lambda: defaultdict(int))
trnalocustrailercounts = defaultdict(lambda: defaultdict(int))
partialtrnalocuscounts = defaultdict(lambda: defaultdict(int))
fulltrnalocuscounts  = defaultdict(lambda: defaultdict(int))

genetypes = dict()


maxoffset = 10
minmapq = 0
if nomultimap:
    minmapq = 2
minreads = 20
#include antisense tRNAs
#also fix thing where tRNAs don't quite overlap enough
#print >>sys.stderr, "multi"
#print >>sys.stderr, nomultimap
for currsample in samples:
    
    currbam = sampledata.getbam(currsample)
    
    try:
        #print >>sys.stderr, currbam
        if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
    
    
    for currfeat in featurelist:
        #try catch is to account for weird chromosomes and the like
        #means that if I can't find a feature, I record no counts for it rather than bailing
        try:
            for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap):
                if currfeat.coverage(currread) > 10:
                    counts[currsample][currfeat.name] += 1
        except ValueError:
            pass
    for genename, featset in itertools.groupby(embllist,lambda x: x.data["genename"]):
        try:
            allreads =set()
            for currfeat in list(featset):
                for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap):
                    allreads.add(currread)
            for currread in allreads:
                if currfeat.coverage(currread) > 10:
                    counts[currsample][genename] += 1 
                    genetypes[genename] = currfeat.data["source"]
                    #print >>sys.stderr, currfeat.bedstring()

        except ValueError:
            pass
    for currfeat in trnaloci:
        for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap):
            if currfeat.coverage(currread) > 10:
                trnalocuscounts[currsample][currfeat.name] += 1
                if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                    fulltrnalocuscounts[currsample][currfeat.name] += 1
                else:# currread.start + fullpretrnathreshold <  currfeat.start or currread.end - fullpretrnathreshold +3 >  currfeat.end:
                    partialtrnalocuscounts[currsample][currfeat.name] += 1
            if currfeat.getdownstream(30).coverage(currread) > 10:
                trnalocustrailercounts[currsample][currfeat.name] += 1
    
    for currfeat in trnalist:
        for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap):
            
            if not currfeat.strand == currread.strand:
                continue
            if not currfeat.coverage(currread) > 10:
                continue
                #print >>sys.stderr, currsample
                #print >>sys.stderr, currread.bedstring()
                #print >>sys.stderr, currfeat.bedstring()
                #print >>sys.stderr, "********"
                
            trnacounts[currsample][currfeat.name] += 1
                
            fragtype = getfragtype(currfeat, currread)
            if fragtype == "Whole":
                trnawholecounts[currsample][currfeat.name] += 1
            elif fragtype == "Fiveprime":
                trnafivecounts[currsample][currfeat.name] += 1
            elif fragtype == "Threeprime":
                trnathreecounts[currsample][currfeat.name] += 1
              
                    #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                    #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                    #print >>sys.stderr, "****"
                        


print "\t".join(samples)

trnanames = set()
for currfeat in trnalist:
    if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
        continue
    if includebase:
        print currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
    else:
        print currfeat.name+"_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currfeat.name]) for currsample in samples)
        print currfeat.name+"_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currfeat.name]) for currsample in samples)
        print currfeat.name+"_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currfeat.name]) for currsample in samples)
        print currfeat.name+"_other\t"+"\t".join(str(trnacounts[currsample][currfeat.name] - (trnathreecounts[currsample][currfeat.name] + trnafivecounts[currsample][currfeat.name] + trnawholecounts[currsample][currfeat.name])) for currsample in samples)
    
    

for currfeat in trnaloci:
    if max(trnalocuscounts[currsample][currfeat.name] for currsample in samples) < minreads:
        continue
    if includebase:
        print currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
    else:
        print currfeat.name+"_wholeprecounts\t"+"\t".join(str(fulltrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        print currfeat.name+"_partialprecounts\t"+"\t".join(str(partialtrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        print currfeat.name+"_trailercounts\t"+"\t".join(str(trnalocustrailercounts[currsample][currfeat.name]) for currsample in samples)

    
    
    
    
for currfeat in featurelist :
    if currfeat.name in trnanames:
        continue
    trnanames.add(currfeat.name)
    if max(counts[currsample][currfeat.name] for currsample in samples) > minreads:
        print currfeat.name+"\t"+"\t".join(str(counts[currsample][currfeat.name]) for currsample in samples)

if args.genetypefile is not None:
    genetypefile = open(args.genetypefile, "w")
    for currfeat in trnaloci:
        print >>genetypefile, currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts"
        print >>genetypefile, currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"
        print >>genetypefile, currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"
        print >>genetypefile, currfeat.name+""+"\t"+"tRNA_locus"
    for currfeat in trnalist:
        print >>genetypefile, currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"
        print >>genetypefile, currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"
        print >>genetypefile, currfeat.name+"_threeprime"+"\t"+"trna_threeprime"
        print >>genetypefile, currfeat.name+"_other"+"\t"+"trna_other"
        print >>genetypefile, currfeat.name+""+"\t"+"tRNA"

    
        
for currfeat in embllist:
    genename = currfeat.data['genename']
    if genename in trnanames:
        continue
    trnanames.add(genename)
    if max(counts[currsample][genename] for currsample in samples) > minreads:
        print genename+"\t"+"\t".join(str(counts[currsample][genename]) for currsample in samples)
        if args.genetypefile is not None:
            print >>genetypefile, genename+"\t"+genetypes[genename]
                
if args.genetypefile is not None:
    genetypefile.close()    
        
        
        
if args.trnacounts is not None:
    trnacountfile = open(args.trnacounts, "w")
    
    for currfeat in trnaloci:
        if max(trnalocuscounts[currsample][currfeat.name] for currsample in samples) < minreads:
            continue
        print >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
    for currfeat in trnalist:
        if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
            continue
        print  >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
    trnacountfile.close()            