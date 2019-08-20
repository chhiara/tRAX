#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import os.path
import re
from trnasequtils import *
from multiprocessing import Pool
import itertools

'''

mergable eteps
variation on raw read read counts
min reads per sample
percentage of reads mapped 60%-70%
read sizes 35-75 for dm-tRNAseq 50-60% ??  15-60 for armseq ???

percentage of reads that map to tRNA  
pca comp ???
number of tRNAs with any reads in at least one sample
number of tRNA reads with reads in at least one sample
whole vs fragment types

check if replicates replicate


'''
def makehtmldoc(outputname):
    #htmlfile = open(outputname, "w")
    print >>sys.stderr, ""

def checkreadprep(prepinfo):
    pass
    

    
class mappingresults:
    def __init__(self, unmap, single, multi):
        self.unmap = unmap
        self.single = single
        self.multi = multi
        
    def totalreadscount(self, sample):
        return self.unmap[sample] + self.single[sample] + self.multi[sample]
    def getmappercent(self, sample):
        totalreads = self.totalreadscount(sample)
        return (totalreads  - self.unmap[sample]) / (1.*totalreads) 
        
minmapreads = 200000
              
minmappercent = .65
def checkreadsmapping(samplename, sampleinfo):
    mapresults = getreadmapping(samplename, sampleinfo)
    samples = sampleinfo.getsamples()
    totalreads = {currsample : mapresults.totalreadscount(currsample) for currsample in samples}
    mappercent = {currsample : mapresults.getmappercent(currsample) for currsample in samples}
    
    lowcountsamples = list(currsample for currsample in samples if totalreads[currsample] < minmapreads)
    if len(lowcountsamples) > 0:
        print >>sys.stderr, str(len(lowcountsamples)) +" contain fewer mappable reads than recommended minimum ("+str(minmapreads)+") [" +",".join(currsample+":"+str(totalreads[currsample]) for currsample in lowcountsamples)+"]"

    lowmapsamples = list(currsample for currsample in samples if mappercent[currsample] < minmappercent)
    if len(lowmapsamples) > 0:
        print >>sys.stderr, str(len(lowmapsamples)) +" have lower mappable read percentage than recommended minimum ("+str(100*minmappercent)+"%) [" +",".join(currsample+":"+str(mappercent[currsample]) for currsample in lowmapsamples)+"]"
def getmapfile(samplename):
    return samplename+"/"+samplename+"-mapinfo.txt"
def getreadmapping(samplename, sampleinfo):
    mappingcounts = dict()
    mapresults = open(getmapfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    multimaps = dict()
    singlemaps = dict()
    unmaps = dict()
    for i, currline in enumerate(mapresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print >>sys.stderr, runsamples
                print >>sys.stderr, allsamples
                print >>sys.stderr, "QAError"
            continue
            
        if len(fields) != len(allsamples) + 1:
            continue
        for j in range(0, len(runsamples)):
            if fields[0] == "unmap":
                unmaps[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "single":
                singlemaps[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "multi":
                multimaps[runsamples[j]] = int(fields[j + 1])
                
    return mappingresults(unmaps, singlemaps, multimaps)

def gettypefile(samplename):
    return samplename+"/"+samplename+"-typerealcounts.txt"
def getreadlengthfile(samplename):
    return samplename+"/"+samplename+"-readlengths.txt"
def getfragtypefile(samplename):
    return samplename+"/"+samplename+"-fragtypes.txt"
def gettrnacountfile(samplename):
    return samplename+"/"+samplename+"-trnacounts.txt"
def getsizefactorfile(samplename):
    return samplename+"/"+samplename+"-SizeFactors.txt"
    
class typecount:
    def __init__(self, typecounts):
        self.typecounts = typecounts
    def gettotal(self, sample):
        return sum(self.typecounts[sample].values()) 
    def gettrnapercent(self, sample):
        return (self.typecounts[sample]["tRNA"] + self.typecounts[sample]["pretRNA"] )/ (1.*self.gettotal(sample))
    def getrrnapercent(self, sample):
        return self.typecounts[sample]["rRNA"] / (1.*self.gettotal(sample))
    def getcountablepercent(self, sample):
        return (self.gettotal(sample) - self.typecounts[sample]["other"]) / (1.*self.gettotal(sample))
def gettypecounts(samplename, sampleinfo):
    typeresults = open(gettypefile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    typecounts = defaultdict(dict)
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print >>sys.stderr, runsamples
                print >>sys.stderr, allsamples
                print >>sys.stderr, "QAError"
            continue
           
        if len(fields) != len(allsamples) + 1:
            print >>sys.stderr, runsamples
            print >>sys.stderr, fields
            print >>sys.stderr, "QAError"
            continue
        for j in range(0, len(runsamples)):
            typecounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return typecount(typecounts)
    
'''
Length	Sample	other	trnas	pretrnas
0	nuc_rep2	0	0	0
1	nuc_rep2	0	0	0
'''
def getmeanfreq(freqtable):
    return sum(curr * freqtable[curr] for curr in freqtable.keys()) / (1.*sum(freqtable.values()))
class lengthcount:
    def __init__(self, trnalengthcounts, pretrnalengthcounts, otherlengthcounts):
        self.trnalengthcounts = trnalengthcounts
        self.pretrnalengthcounts = pretrnalengthcounts
        self.otherlengthcounts = otherlengthcounts
        self.samples = set(itertools.chain(trnalengthcounts.keys(), pretrnalengthcounts.keys(), otherlengthcounts.keys()))
        
        #print >>sys.stderr, list(itertools.chain(list(trnalengthcounts[currsample].keys() for currsample in self.samples),list(pretrnalengthcounts[currsample].keys() for currsample in self.samples),list(otherlengthcounts[currsample].keys() for currsample in self.samples)))

        
        self.maxlength = max(itertools.chain(itertools.chain.from_iterable(trnalengthcounts[currsample].keys() for currsample in self.samples),itertools.chain.from_iterable(pretrnalengthcounts[currsample].keys() for currsample in self.samples),itertools.chain.from_iterable(otherlengthcounts[currsample].keys() for currsample in self.samples)))
        
    def getalllengths(self, sample):
        return {currlength: self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength] for currlength in range(self.maxlength)}
    def getsamplemean(self, sample):
        return getmeanfreq(self.getalllengths(sample))
    def getthreshold(self, sample, minsize, maxsize):
        alllengths = self.getalllengths(sample)
        return sum(alllengths[i] for i in range(minsize, maxsize))
    def getthresholdpercent(self, sample, minsize, maxsize):
        
        return self.getthreshold(sample, minsize, maxsize) / (1.*sum(self.getalllengths(sample).values()))

def getreadlengths(samplename,sampleinfo):
    lengthresults = open(getreadlengthfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    trnalengthcounts = defaultdict(lambda: defaultdict(int))
    otherlengthcounts = defaultdict(lambda: defaultdict(int))
    pretrnalengthcounts = defaultdict(lambda: defaultdict(int))
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(lengthresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            continue
           
        if len(fields) < 4:
            continue
            pass
        trnalengthcounts[fields[1]][int(fields[0])] = int(fields[3])
        otherlengthcounts[fields[1]][int(fields[0])] = int(fields[2])
        pretrnalengthcounts[fields[1]][int(fields[0])] = int(fields[4])
    return lengthcount(trnalengthcounts, pretrnalengthcounts, otherlengthcounts)
    



def getfragtypes(samplefile, trnainfo):
    pass

trnapercentcutoff = .05
ribopercentcutoff = .35
unmappercentcutoff = .35

highmeanlength = 75
minsizethreshold = 15
maxsizethreshold = 75
percentsizethreshold = 70

def checkreadtypes(samplename, sampleinfo):
    typecounts = gettypecounts(samplename, sampleinfo)
    getfragtypes(samplename, sampleinfo)
    
    samples = sampleinfo.getsamples()
    
    trnapercent = {currsample : typecounts.gettrnapercent(currsample) for currsample in samples}
    lowtrnasamples = list(currsample for currsample in samples if trnapercent[currsample] < trnapercentcutoff)
    if len(lowtrnasamples) > 0:
        print >>sys.stderr, str(len(lowtrnasamples)) +" samples have low tRNA read percentage ( < "+str(100*trnapercentcutoff)+"%) [" +",".join(currsample+":"+str(trnapercent[currsample]) for currsample in lowtrnasamples)+"]"
        
    rrnapercent = {currsample : typecounts.getrrnapercent(currsample) for currsample in samples}
    highribosamples = list(currsample for currsample in samples if rrnapercent[currsample] > ribopercentcutoff)
    if len(highribosamples) > 0:
        print >>sys.stderr, str(len(highribosamples)) +" samples have high rRNA read percentage ( > "+str(100*ribopercentcutoff)+"%) [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in highribosamples)+"]"
        
    otherpercent = {currsample : typecounts.getcountablepercent(currsample) for currsample in samples}
    highothersamples = list(currsample for currsample in samples if otherpercent[currsample] > unmappercentcutoff)
    if len(highothersamples) > 0:
        print >>sys.stderr, str(len(highothersamples)) +" samples have many reads not mapping to annotated genes ( > "+str(100*unmappercentcutoff)+"%) [" +",".join(currsample+":"+str(otherpercent[currsample]) for currsample in highothersamples)+"]"


    allreadlength = getreadlengths(samplename, sampleinfo)

    meanreadlength = {currsample : allreadlength.getsamplemean(currsample) for currsample in samples}
    meanreadsamples = list(currsample for currsample in samples if rrnapercent[currsample] > ribopercentcutoff)
    if len(meanreadsamples) > 0:
        print >>sys.stderr, str(len(meanreadsamples)) +" samples have high read length average ( > "+str(highmeanlength)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in meanreadsamples)+"]"
        
    thresholdreadpercent = {currsample : allreadlength.getthresholdpercent(currsample, minsizethreshold, maxsizethreshold) for currsample in samples}
    badsizesamples = list(currsample for currsample in samples if rrnapercent[currsample] < percentsizethreshold)
    if len(highribosamples) > 0:
        print >>sys.stderr, str(len(badsizesamples)) +" samples have high rRNA read percentage ( > "+str(100*percentsizethreshold)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in badsizesamples)+"]"



class trnacount:
    def __init__(self, trnacounts):
        self.trnacounts = trnacounts
    def gettrnaactive(currsample, cutoff = 20):
        return sum(1 for curr in self.trnacounts[currsample].iterkeys() if curr > cutoff)
    def gettrnaactivepercent(currsample, trnainfo, cutoff = 20):
        return self.gettrnaactive(currsample, cutoff)/ (1.*trnainfo.getalltrnas())
        
        
def gettrnacounts(samplename, sampleinfo):
    typeresults = open(gettrnacountfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    trnacounts = defaultdict(dict)
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print >>sys.stderr, runsamples
                print >>sys.stderr, allsamples
                print >>sys.stderr, "QAError"
            continue
           
        if len(fields) != len(allsamples) + 1:
            print >>sys.stderr, runsamples
            print >>sys.stderr, fields
            print >>sys.stderr, "QAError"
            continue
        for j in range(0, len(runsamples)):
            trnacounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return trnacount(trnacounts)
class sizefactor:
    def __init__(self, sizefactors):
        self.sizefactors = sizefactors

def getsizefactor(samplename, sampleinfo):
    typeresults = open(getsizefactorfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    sizefactors = dict()
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print >>sys.stderr, runsamples
                print >>sys.stderr, allsamples
                print >>sys.stderr, "QAError"
            continue
           
        if len(fields) != len(allsamples) + 1:
            print >>sys.stderr, runsamples
            print >>sys.stderr, fields
            print >>sys.stderr, "QAError"
            continue
        for j in range(0, len(runsamples)):
            sizefactors[runsamples[j]] = float(fields[j + 1])
    return sizefactor(sizefactors)
    
sizefactordiff = 3
def checkgenecounts(samplename, sampleinfo, trnainfo):
    print >>sys.stderr, "**"+samplename
    readcounts = gettrnacounts(samplename, sampleinfo)
    sizefactors = getsizefactor(samplename, sampleinfo)
    
    thresholdreadpercent = {currsample : allreadlength.getthresholdpercent(currsample, minsizethreshold, maxsizethreshold) for currsample in samples}
    badsizesamples = list(currsample for currsample in samples if rrnapercent[currsample] < percentsizethreshold)
    if len(highribosamples) > 0:
        print >>sys.stderr, str(len(badsizesamples)) +" samples have low tRNA read counts ( > "+str(100*percentsizethreshold)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in badsizesamples)+"]"

    if min(sizefactors.sizefactors.values())*sizefactordiff <  min(sizefactors.sizefactors.values()):
        print >>sys.stderr, "Large DESeq2 sizefactor differences ( >"+sizefactordiff+" ) [" +",".join(currsample+":"+str(sizefactors.sizefactors[currsample]) for currsample in badsizesamples)+"]"

def checkreadprep(sampleinfo):
    pass


def checktrnamappings(sampleinfo):
    pass



def checkfragmenttypes(samplefile, trnainfo):
    pass

def main(**args):
    samplename = args["experimentname"]
    sampleinfo = samplefile(os.path.expanduser(args["samplefile"]))
    trnainfo = transcriptfile(os.path.expanduser(args["databasename"] + "-trnatable.txt"))
    readcounts = checkreadprep(sampleinfo)
    readcounts = checkreadsmapping(samplename, sampleinfo)
    readcounts = checkreadtypes(samplename, sampleinfo)
    checkgenecounts(samplename, sampleinfo, trnainfo)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')
    
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--databasename',
                       help='tRNA file in format')
    parser.add_argument('--experimentname',
                       help='Sample file in format')

    
    

    
    args = parser.parse_args()
    main(**vars(args))

