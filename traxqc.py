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
    
def percentform(innum):
    return  "{0:.2f}%".format(100 * innum)
def countform(innum):
    return "{0:.2f}".format(innum)
#print  str(len(missingtrnasamples)) +" samples have low tRNA read counts ( > "+str(100*minactivepercent)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(thresholdreadpercent[currsample]) for currsample in missingtrnasamples)+"]"
def errorline(faillist, failmessage, failcriteria, faildict, percentformat = True):
    if percentformat:
        outformat = percentform
    else:
        outformat = countform
    print "<p>"
    print  str(len(faillist)) +" samples " + failmessage +" ( "+failcriteria+") [" +",".join(currsample+":"+outformat(faildict[currsample]) for currsample in faillist)+"]"
    print "</p>"

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
    #print str(len(lowcountsamples)) +" contain fewer mappable reads than recommended minimum ("+str(minmapreads)+") [" +",".join(currsample+":"+str(totalreads[currsample]) for currsample in lowcountsamples)+"]"
    errorline(lowcountsamples, "contain fewer mappable reads than recommended minimum", " < "+str(minmapreads)+"", lowcountsamples)

    lowmapsamples = list(currsample for currsample in samples if mappercent[currsample] < minmappercent)
    errorline(lowmapsamples, "have lower mappable read percentage than recommended minimum ", " < "+str(100*minmappercent)+"%", lowmapsamples)
    #print str(len(lowmapsamples)) +" have lower mappable read percentage than recommended minimum ("+str(100*minmappercent)+"%) [" +",".join(currsample+":"+str(mappercent[currsample]) for currsample in lowmapsamples)+"]"
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
    def getotherpercent(self, sample):
        return (self.typecounts[sample]["other"]) / (1.*self.gettotal(sample))
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
percentsizethreshold = .70

def checkreadtypes(samplename, sampleinfo):
    typecounts = gettypecounts(samplename, sampleinfo)
    getfragtypes(samplename, sampleinfo)
    
    samples = sampleinfo.getsamples()
    
    trnapercent = {currsample : typecounts.gettrnapercent(currsample) for currsample in samples}
    lowtrnasamples = list(currsample for currsample in samples if trnapercent[currsample] < trnapercentcutoff)
    #print  str(len(lowtrnasamples)) +" samples have low tRNA read percentage ( < "+str(100*trnapercentcutoff)+"%) [" +",".join(currsample+":"+str(trnapercent[currsample]) for currsample in lowtrnasamples)+"]"
    errorline(lowtrnasamples, "samples have low tRNA read percentage", " < "+str(100*trnapercentcutoff)+"%", trnapercent)

    rrnapercent = {currsample : typecounts.getrrnapercent(currsample) for currsample in samples}
    highribosamples = list(currsample for currsample in samples if rrnapercent[currsample] > ribopercentcutoff)
    #print str(len(highribosamples)) +" samples have high rRNA read percentage ( > "+str(100*ribopercentcutoff)+"%) [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in highribosamples)+"]"
    errorline(highribosamples, "have high rRNA read percentage", " > "+str(100*ribopercentcutoff)+"%", rrnapercent)

    otherpercent = {currsample : typecounts.getotherpercent(currsample) for currsample in samples}
    highothersamples = list(currsample for currsample in samples if otherpercent[currsample] > unmappercentcutoff)
    #print  str(len(highothersamples)) +" samples have many reads not mapping to annotated genes ( > "+str(100*unmappercentcutoff)+"%) [" +",".join(currsample+":"+str(otherpercent[currsample]) for currsample in highothersamples)+"]"
    errorline(highothersamples, "have many reads not mapping to annotated genes", "> "+str(100*unmappercentcutoff)+"%", otherpercent)

    allreadlength = getreadlengths(samplename, sampleinfo)

    meanreadlength = {currsample : allreadlength.getsamplemean(currsample) for currsample in samples}
    meanreadsamples = list(currsample for currsample in samples if rrnapercent[currsample] > highmeanlength)
    #print  str(len(meanreadsamples)) +" samples have high read length average ( > "+str(highmeanlength)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in meanreadsamples)+"]"
    errorline(meanreadsamples, " have high read length average", " > "+str(highmeanlength), meanreadlength, percentformat = False)
    thresholdreadpercent = {currsample : allreadlength.getthresholdpercent(currsample, minsizethreshold, maxsizethreshold) for currsample in samples}
    badsizesamples = list(currsample for currsample in samples if thresholdreadpercent[currsample] < percentsizethreshold)
    #print  str(len(badsizesamples)) +" samples have high rRNA read percentage ( > "+str(100*percentsizethreshold)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in badsizesamples)+"]"
    errorline(badsizesamples, "do not fit tRNA sizes", " > "+str(100*percentsizethreshold)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold), thresholdreadpercent)




class trnacount:
    def __init__(self, trnacounts):
        self.trnacounts = trnacounts
    def gettrnaactive(self, currsample, cutoff = 20):
        return sum(1 for curr in self.trnacounts[currsample].iterkeys() if curr > cutoff)
    def gettrnaactivepercent(self, currsample, trnainfo, cutoff = 20):
        #print >>sys.stderr, self.gettrnaactive(currsample, cutoff)
        #print >>sys.stderr, (1.*len(trnainfo.gettranscripts()))
        #print >>sys.stderr, self.trnacounts[currsample].keys()
        return self.gettrnaactive(currsample, cutoff)/ (1.*len(trnainfo.gettranscripts()))
        
        
def gettrnacounts(samplename, sampleinfo, trnainfo):
    typeresults = open(gettrnacountfile(samplename))
    allsamples = sampleinfo.getsamples()
    trnatranscripts = set(trnainfo.gettranscripts())
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
        if fields[0] in trnatranscripts:    
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
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split()
        fields = list(curr.strip('"') for curr in fields)
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print >>sys.stderr, list(runsamples)
                print >>sys.stderr, list(allsamples)
                print >>sys.stderr, "QAError"
            continue
           
        if len(fields) != len(allsamples):
            print >>sys.stderr, len(runsamples)
            print >>sys.stderr, len(fields)
            print >>sys.stderr, "QAError"
            continue
        #print >>sys.stderr, "**"
        for j in range(0, len(runsamples)):
            sizefactors[runsamples[j]] = float(fields[j])
    return sizefactor(sizefactors)
    
sizefactordiff = 3

minactivepercent = .5

minreadcount = 20
def checkgenecounts(samplename, sampleinfo, trnainfo):
    #print >>sys.stderr, "**"+samplename
    readcounts = gettrnacounts(samplename, sampleinfo, trnainfo)
    sizefactors = getsizefactor(samplename, sampleinfo)
    samples = sampleinfo.getsamples()
    
    thresholdreadpercent = {currsample : readcounts.gettrnaactivepercent(currsample, trnainfo, cutoff = minreadcount) for currsample in samples}
    
    
    missingtrnasamples = list(currsample for currsample in samples if thresholdreadpercent[currsample] < minactivepercent)
    #print  str(len(missingtrnasamples)) +" samples have low tRNA read counts ( > "+str(100*minactivepercent)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(thresholdreadpercent[currsample]) for currsample in missingtrnasamples)+"]"
    errorline(missingtrnasamples, "have low tRNA read counts", " > "+str(100*minactivepercent)+"% less than "+str(minreadcount) + " reads", thresholdreadpercent)
    #if min(sizefactors.sizefactors.values())*sizefactordiff <  min(sizefactors.sizefactors.values()):
    #print  "Large DESeq2 sizefactor differences ( >"+sizefactordiff+" ) [" +",".join(currsample+":"+str(sizefactors.sizefactors[currsample]) for currsample in badsizesamples)+"]"

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
    

    print "<html>"
    print "<body>"

    checkreadprep(sampleinfo)
    checkreadsmapping(samplename, sampleinfo)
    checkreadtypes(samplename, sampleinfo)
    checkgenecounts(samplename, sampleinfo, trnainfo)
    print "<body>"
    print "</html>"



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

