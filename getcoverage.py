#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *


gapchars = set("-._~")
class readcoverage:
    def __init__(self, region):
        self.region = region
        self.samplereads = 0
        self.coverage = list()
        self.length = region.length()
        self.totalreads = 0
        for i in range(0,region.length()):
            self.coverage.append(0)
    def addread(self, read):
        self.totalreads += 1
        if self.region.strand == "+":
            start = max([0, read.start - self.region.start])
            end = min([self.length, read.end - self.region.start])
            
        else:
            start = max([0, self.region.end - read.end])
            end = min([self.length, self.region.end - read.start])
        
        for currpos in range(self.length):
            if start <= currpos <= end - 1:
                self.coverage[currpos] += 1
                
    def coveragelist(self):
        return self.coverage
    def coveragealign(self, alignment, gapoutput = "NA",sizefactor = 1):
        if len(self.coverage) != len(string.translate(alignment, None, str(gapchars))):
            print >>sys.stderr, "Alignment length does not match bed length"            
        i = 0
        for curr in alignment:
            #print >>sys.stderr, curr
            if curr in gapchars:
                yield gapoutput
            else:
                yield self.coverage[i]/sizefactor
                i += 1
def nasum(operands, naval = "NA"):
    if sum("NA"== curr for curr in operands) == len(operands):
        return "NA"
    elif not any("NA"== curr for curr in operands):
        return sum(operands)
    else:
        print >>sys.stderr, "Trying to add incompatible alignments"
        sys.exit(1)
    #return ",".join(str(curr) for curr in operands)
def sumsamples(coverage,sampledata, repname, currfeat, sizefactors = defaultdict(lambda: 1)):
    return (nasum(curr) for curr in itertools.izip(*(allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]) for currsample in sampledata.getrepsamples(repname))))
    
count = 0

positions = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'-',18,19,20,'-','-',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
#print >>sys.stderr, len(positions)
#this gets the tRNA numbers by the sprinzel numbering system
def gettnanums(trnaalign, margin = 0):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    for i, struct in enumerate(trnaalign.consensus):
        if currcount >= len(positions):
            trnanum.append('gap'+str(gapnum))
            gapnum += 1
            currcount += 1
        elif struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=':
                currcount = 1
            if positions[currcount] == 'e':
                trnanum.append('e'+str(enum))
                enum += 1
                currcount += 1
            elif positions[currcount] == '-':
                trnanum.append('gap'+str(gapnum))
                gapnum += 1
                currcount += 1
            else:
                trnanum.append(str(positions[currcount]))
                currcount += 1
        else:
            #if intron

            if positions[currcount] == 38:
                trnanum.append('intron'+str(gapnum))
                intronnum += 1
            else:
                
                trnanum.append('gap'+str(gapnum))
                gapnum += 1
    for i in range(margin):
        trnanum.append('tail'+str(i+1))
    return trnanum
def main(**argdict):
    #print >>sys.stderr, argdict
    argdict = defaultdict(lambda: None, argdict)
    if "edgemargin" not in  argdict:
        edgemargin = 0
    else:
        edgemargin = int(argdict["edgemargin"])
    #currently crashes if set to zero
    if "mincoverage" not in  argdict:
        mincoverage = 10
    else:
        mincoverage = int(argdict["mincoverage"])  
    
        
    sampledata = samplefile(argdict["samplefile"])

    maxmismatches = argdict["maxmismatches"]
    uniquename = argdict["uniquename"]
    uniquegenome = argdict["uniquegenome"]
    trnastk = list(readrnastk(open(argdict["stkfile"], "r")))[0]
    bedfile = argdict["bedfile"]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print >>sys.stderr, "Size factor file "+argdict["sizefactors"]+" missing "+currsample
                sys.exit(1)
    combinereps = argdict["combinereps"]
    allcoveragefile = None
    if "allcoverage" not in argdict or argdict["allcoverage"] == "stdout":
        allcoveragefile = sys.stdout
    else:
        allcoveragefile = open(argdict["allcoverage"],"w")
    samples = sampledata.getsamples()
    minextend = None
    if argdict["minextend"]:
        minextend = int(argdict["minextend"])

    alltrnas = list()
    #gettnanums
    
    positionnums = gettnanums(trnastk, margin = edgemargin)
    trnastk = trnastk.addmargin(edgemargin)
    
    try:
        basetrnas = list()
        for currfile in bedfile:
            basetrnas.extend(list(currbed for currbed in readbed(currfile)))       
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    
    trnalist = list(curr.addmargin(edgemargin) for curr in basetrnas)

    featcount = defaultdict(int)

    maxoffset = 10

    allcoverages = dict()
    multaminocoverages = dict()
    multaccoverages = dict()
    multtrnacoverages = dict()
    uniquecoverages = dict()    
    uniquegenomecoverages = dict()  
    multigenomecoverages = dict()
    for currsample in samples:
        currbam = sampledata.getbam(currsample)
        allcoverages[currsample] = dict()
        multaminocoverages[currsample] = dict()
        multaccoverages[currsample] = dict()
        multtrnacoverages[currsample] = dict()
        uniquecoverages[currsample] = dict()
        uniquegenomecoverages[currsample] = dict()
        multigenomecoverages[currsample] = dict()
        try:
            #print >>sys.stderr, currbam
            if not os.path.isfile(currbam+".bai"):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit()
            
        for i, currfeat in enumerate(basetrnas):
            allcoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaminocoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaccoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multtrnacoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquegenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multigenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            if trnalist[i].name == "tRNA-Leu-AAG-1-2":
                pass
                #print >>sys.stderr, "**|"
            for currread in getbamrange(bamfile, trnalist[i],maxmismatches = maxmismatches,allowindels = False):

                if  trnalist[i].coverage(currread) > 10:
                    if minextend is not None and not (currread.start + minextend <= trnalist[i].start or currread.end - minextend >= trnalist[i].end):
                        continue
                    if trnalist[i].name == "tRNA-Ile-AAT-5-5":
                        pass
                        #print >>sys.stderr, "**"
                    if currfeat.coverage(currread) > 10:
                        #continue
                        pass

                    
                    allcoverages[currsample][trnalist[i].name].addread(currread)
                    if not isuniqueaminomapping(currread):
                        multaminocoverages[currsample][trnalist[i].name].addread(currread)
                    elif not isuniqueacmapping(currread):
                        multaccoverages[currsample][trnalist[i].name].addread(currread)
                    elif not isuniquetrnamapping(currread):
                        multtrnacoverages[currsample][trnalist[i].name].addread(currread)
                    else:
                        uniquecoverages[currsample][trnalist[i].name].addread(currread)
                    if issinglemapped(currread):
                        uniquegenomecoverages[currsample][trnalist[i].name].addread(currread)
                    else:
                        multigenomecoverages[currsample][trnalist[i].name].addread(currread)
                    
                else:
                    pass
    covfiles = dict()
    print >>allcoveragefile, "Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
        
    for currfeat in trnalist:
        totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)

        if totalreads < mincoverage:
            continue
        if currfeat.name == "tRNA-Leu-AAG-1-2":
            print >>sys.stderr, "total"
            print >>sys.stderr, totalreads
        if combinereps:
    
            replicates = sampledata.allreplicates()
            for currrep in replicates:
                print >>allcoveragefile, currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(allcoverages,sampledata,currrep,currfeat,sizefactor))
            
        else:
            for currsample in samples:
                print >>allcoveragefile, currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]))
    if uniquegenome:
        covfiles = {uniquegenome + '-uniquegenomecoverages.txt':uniquegenomecoverages,uniquegenome + '-multgenomecoverages.txt':multigenomecoverages}

        for filename, currcoverage in covfiles.iteritems():
            covfile = open(filename, "w")
            print >>covfile,"Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
                
            for currfeat in trnalist:
                totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
                if totalreads < mincoverage:
                    continue
                if combinereps:
            
                    replicates = sampledata.allreplicates()
                    for currrep in replicates:
                        print  >>covfile,currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor))
                    
                else:
                    for currsample in samples:
                        print  >>covfile,currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in currcoverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]))
            covfile.close()
    if uniquename:
        covfiles = {uniquename + '-uniquecoverages.txt':uniquecoverages,uniquename + '-multaccoverages.txt':multaccoverages,uniquename + '-multtrnacoverages.txt':multtrnacoverages,uniquename + '-multaminocoverages.txt':multaminocoverages}
        
        for filename, currcoverage in covfiles.iteritems():
            covfile = open(filename, "w")
            print >>covfile,"Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
                
            for currfeat in trnalist:

                totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)

                if totalreads < mincoverage:
                    continue
                if combinereps:
            
                    replicates = sampledata.allreplicates()
                    for currrep in replicates:
                        print  >>covfile,currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor))
                    
                else:
                    for currsample in samples:
                        print  >>covfile,currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in currcoverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]))
            covfile.close()
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--bedfile',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--stkfile',
                       help='Stockholm file')
    parser.add_argument('--sizefactors',
                       help='Optional file including size factors that will be used for normalization')
    parser.add_argument('--combinereps', action="store_true", default=False,
                       help='Sum samples that are replicates')
    parser.add_argument('--edgemargin', type=int, default=0,
                       help='margin to add to feature coordinates')
    
    parser.add_argument('--mincoverage', type=int, default=10,
                       help='Reads with less then this are filtered out (default 10)')
    parser.add_argument('--uniquename',
                       help='Name for files showing unique and non-unique tRNA reads')
    parser.add_argument('--uniquegenome',
                       help='Name for files showing unique and non-unique genome reads')
    parser.add_argument('--maxmismatches', default=None,
                       help='Set maximum number of allowable mismatches')
    '''
    parser.add_argument('--trnapositions', action="store_true", default=False,
                       help='Use tRNA positions')
    '''
    
    
    '''
    Perform check on sizefactor file to ensure it has all samples
    '''
    args = parser.parse_args()
    argdict = vars(args)
    main(**argdict)
