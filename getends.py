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
    def coveragealign(self, alignment, gapoutput = None,sizefactor = 1):
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

def getrepreadcounts(readcounts, repname, featname, sampledata):
    return sum(readcounts[currsample][featname] for currsample in sampledata.getrepsamples(currrep))
count = 0

positions = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'-',18,19,20,'-','-',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
#99 long
#positions = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'-',18,19,20,'-','-',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,71,72,73])
#print >>sys.stderr, len(positions)
def gettnanums(trnaalign, margin = 0):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    #print >>sys.stderr, len(list(curr for curr in trnaalign.consensus if curr != "."))
    #print >>sys.stderr, "".join(list(curr for curr in trnaalign.consensus if curr != "."))
    #print >>sys.stderr, "".join(list(curr for curr in trnaalign.consensus if curr != "."))

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
    #print >>sys.stderr, mincoverage
    #print >>sys.stderr, "***"
    
    if "pairfile" not in argdict or argdict["pairfile"] is None: 
        pairfiles = list()
        #sys.exit()
    else:
        pairfiles = list()
        for currline in open(argdict["pairfile"]):
            pairfiles.append(tuple(currline.rstrip().split()))
            
    sampledata = samplefile(argdict["samplefile"])
    
    maxmismatches = argdict["maxmismatches"]
    uniquename = argdict["uniquename"]
    uniquegenome = argdict["uniquegenome"]
    mismatchfilename = argdict["mismatchfile"]
    trnastk = list(readrnastk(open(argdict["stkfile"], "r")))[0]
    bedfile = argdict["bedfile"]
    sizefactor = defaultdict(lambda:1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"])
    combinereps = argdict["combinereps"]
    allcoveragefile = None
    if "allcoverage" not in argdict or argdict["allcoverage"] == "stdout":
        allcoveragefile = open(os.devnull,"w") #sys.stdout
    else:
        allcoveragefile = open(argdict["allcoverage"],"w")
        
    if "mismatchreport" not in argdict or argdict["mismatchreport"] == "stdout":
        mismatchreport = sys.stdout
    else:
        mismatchreport = open(argdict["mismatchreport"],"w")    
    samples = sampledata.getsamples()
    minextend = None
    if argdict["minextend"]:
        minextend = int(argdict["minextend"])

    alltrnas = list()
    #gettnanums
    trnafasta = argdict["trnafasta"]
    trnaseqs = fastadict(trnafasta)
    for currname in trnaseqs.keys():
        trnaseqs[currname] = ("N"*edgemargin)+trnaseqs[currname]+("N"*edgemargin)
    #print >>sys.stderr, trnaseqs['tRNA-Pro-CGG-2']
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
    
    #./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
        
    
    featcount = defaultdict(int)
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
    #print >>sys.stderr, "***"
    #lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    
    
    maxoffset = 10
    
    '''
    
    cmbuild --enone --hand trnamature-euk.cm MaturetRNAs.stk
    
    transcript
    .(((.(.(.((..,....,..<<.<<.___..____..._>>>>,.<...<<<<.__..__.___....>>>.>>,,..............,,,.<<<<<._______>>>.>...>...))....)))))::::
    
    loci:
    ((..(..(...(.....(..(.....,....,.<<<<____.___...._>>>>.........................,.<...<.<.<<...___.__.._................................................................................................._>>>.>.....>,,<<<<<<<____.>>>>>>>,..,<<<.<<._______...>>>..>..>....))....).)...))):
    
    loop = [\.\_]+
    openstem = [\.\(\<]+?
    closestem = [\.\)\>]+?
    inter = [\.\,]
    '''
    #transcriptalign = re.compile(r"\.[\.\(]+?[]") 
    readcounts = defaultdict(lambda: defaultdict(lambda: 1))
                
                
    
    def nasum(operands, naval = "NA"):
        if sum("NA"== curr or curr is None for curr in operands) == len(operands):
            return None
        elif not any("NA"== curr for curr in operands):
            return sum(curr if curr != "NA" or curr is None else 0 for curr in operands)
        else:
            print >>sys.stderr, "Trying to add incompatible alignments"
            sys.exit(1)
        #return ",".join(str(curr) for curr in operands)
    def sumsamples(coverage,sampledata, repname, currfeat, sizefactors = defaultdict(lambda: 1)):
        return (nasum(curr) for curr in itertools.izip(*(coverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]) for currsample in sampledata.getrepsamples(repname))))
    
        
    allcoverages = dict()
    multaminocoverages = dict()
    multaccoverages = dict()
    multtrnacoverages = dict()
    uniquecoverages = dict()    
    uniquegenomecoverages = dict()  
    multigenomecoverages = dict()
    readmismatches = dict()
    
    adeninemismatches = dict()
    thyminemismatches = dict()
    cytosinemismatches = dict()
    guanosinemismatches = dict()        
    
    
    for currsample in samples:
        currbam = sampledata.getbam(currsample)
        allcoverages[currsample] = dict()
        multaminocoverages[currsample] = dict()
        multaccoverages[currsample] = dict()
        multtrnacoverages[currsample] = dict()
        uniquecoverages[currsample] = dict()
        uniquegenomecoverages[currsample] = dict()
        multigenomecoverages[currsample] = dict()
        readmismatches[currsample] = dict()
        
        adeninemismatches[currsample] = dict()
        thyminemismatches[currsample] = dict()
        cytosinemismatches[currsample] = dict()
        guanosinemismatches[currsample] = dict()        
        
        
        try:
            #print >>sys.stderr, currbam
            if not os.path.isfile(currbam+".bai"):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit()
        #print >>sys.stderr, currsample
        for i, currfeat in enumerate(basetrnas):
            #print >>sys.stderr, currfeat.name
            allcoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaminocoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaccoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multtrnacoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquegenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multigenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            readmismatches[currsample][currfeat.name] = readcoverage(trnalist[i])
            
            adeninemismatches[currsample][currfeat.name] =   readcoverage(trnalist[i])
            thyminemismatches[currsample][currfeat.name] =   readcoverage(trnalist[i])
            cytosinemismatches[currsample][currfeat.name] =  readcoverage(trnalist[i])
            guanosinemismatches[currsample][currfeat.name] = readcoverage(trnalist[i])
        
        
            for currread in getbamrange(bamfile, currfeat,maxmismatches = maxmismatches,allowindels = False):
                if  basetrnas[i].coverage(currread) > 10:

                    if minextend is not None and not (currread.start + minextend <= currfeat.start or currread.end - minextend >= currfeat.end):
                        continue
                    readcounts[currsample][trnalist[i].name] += 1
                    readcounts[sampledata.getreplicatename(currsample)][trnalist[i].name] += 1  
                    readstart = currread.getfirst(1)
                    #allcoverages[currsample][trnalist[i].name].addread(readstart)
                    if not isuniqueaminomapping(currread):
                        multaminocoverages[currsample][trnalist[i].name].addread(readstart)
                    elif not isuniqueacmapping(currread):
                        multaccoverages[currsample][trnalist[i].name].addread(readstart)
                    elif not isuniquetrnamapping(currread):
                        multtrnacoverages[currsample][trnalist[i].name].addread(readstart)
                    else:
                        uniquecoverages[currsample][trnalist[i].name].addread(readstart)
                    if issinglemapped(currread):
                        uniquegenomecoverages[currsample][trnalist[i].name].addread(readstart)
                    else:
                        multigenomecoverages[currsample][trnalist[i].name].addread(readstart)
                    for currpos in range(len(currread.data["seq"])):

                        if trnaseqs[currfeat.name][currread.start+ currpos] != currread.data["seq"][currpos]:
                            #print >>sys.stderr, readmismatches[currsample].keys()
                            #print >>sys.stderr, allcoverages[currsample].keys()
                            readmismatches[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                            if 10 < currpos < 30:
                                pass
                                #print >>sys.stderr, trnaseqs[currfeat.name][currread.start+ currpos-5:currread.start+ currpos+5]
                                #print >>sys.stderr, currread.data["seq"][currpos-5:currpos+5]
                                #print >>sys.stderr, currpos
                        allcoverages[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                        if currread.data["seq"][currpos] == "A":
                            adeninemismatches[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                        elif currread.data["seq"][currpos] == "T":
                            thyminemismatches[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                        elif currread.data["seq"][currpos] == "C":
                            cytosinemismatches[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                        elif currread.data["seq"][currpos] == "G":
                            guanosinemismatches[currsample][trnalist[i].name].addread(currread.getbase(currpos))
                else:
                    pass
                    #print >>sys.stderr, "***"
    replicates = sampledata.allreplicates()
    covfiles = dict()
    for currpair in pairfiles:
        for currfeat in trnalist:
            totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
            if totalreads < mincoverage:
                continue
            if currpair[0] in replicates:
                fircounts = list(curr/(1.*readcounts[currpair[0]][currfeat.name]) if curr is not None else 0 for curr in sumsamples(allcoverages,sampledata,currpair[0],currfeat))
            elif currpair[0] in samples:
                fircounts = list(curr/(1.*readcounts[currpair[0]][currfeat.name]) if curr is not None else 0 for curr in allcoverages[currpair[0]][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            if currpair[1] in replicates:
                seccounts = list(curr/(1.*readcounts[currpair[1]][currfeat.name]) if curr is not None else 0 for curr in sumsamples(allcoverages,sampledata,currpair[1],currfeat))
            elif currpair[1] in samples:
                seccounts = list(curr/(1.*readcounts[currpair[1]][currfeat.name]) if curr is not None else 0 for curr in allcoverages[currpair[1]][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            for i, currpos in enumerate(positionnums):
                modpos = max([i-1,0])
                if seccounts[i] > fircounts[i]*2 and seccounts[i] > .1:
                    #print currpair[0]+":"+currpair[1]
                    #print currfeat.name+":"+currpos+":"+str(modpos)+","+str(modpos)
                    pass
                elif fircounts[i] > seccounts[i]*2 and fircounts[i] > .1:
                    #print currpair[0]+":"+currpair[1]
                    #print currfeat.name+":"+currpos+":"+str(modpos)+","+str(modpos)
                    pass
    print >>mismatchreport, "tRNA_name\tsample\tposition\tpercentmismatch\tcoverage\tends\ttRNAreadstotal\tactualbase\tmismatchedbases\tadenines\tthymines\tcytosines\tguanines"                            
    for currfeat in trnalist:
        
        totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
        if totalreads < mincoverage:
            continue
        reportpositions = set()
        for currsample in samples:
            #mismatchcounts = list(curr if curr is not None else 0 for curr in sumsamples(readmismatches,sampledata,currrep,currfeat))
            #covcounts = list(curr if curr is not None else 0 for curr in sumsamples(allcoverages,sampledata,currrep,currfeat))
            
            
            covcounts = list(curr/(1.*readcounts[currsample][currfeat.name]) if curr is not None else 0 for curr in readmismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            #covcounts = list(curr/(1.*readcounts[currsample][currfeat.name]) if curr is not None else "NA" for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            mismatches =  list(curr if curr is not None else 0 for curr in readmismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            
            allstarts  = list(curr if curr is not None else 0 for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            allcovcount  = list(curr if curr is not None else 0 for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            #covcounts = list(curr/(1.*readcounts[currrep][currfeat.name]) if curr is not None else 1 for curr in sumsamples(readmismatches,sampledata,currrep,currfeat))
            #mismatchfractions = list(curr  if curr is not None else 0 for curr in sumsamples(currcoverage,sampledata,currrep,currfeat))
            #mismatchfractions = list(curr/(1.*covcounts[i]+10) if curr is not None else 0 for i, curr in enumerate(sumsamples(allcoverages,sampledata,currrep,currfeat)))
            adeninecount  = list(curr if curr is not None else 0 for curr in adeninemismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            thyminecount = list(curr if curr is not None else 0 for curr in thyminemismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            cytosinecount = list(curr if curr is not None else 0 for curr in cytosinemismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            guanosinecount  = list(curr if curr is not None else 0 for curr in guanosinemismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            for i, currcount in enumerate(covcounts):
                #modpos = max([i-1,0])
                mismatchthreshold = .1
                
                #print >>sys.stderr, mismatchcounts
                #print >>sys.stderr, covcounts
                #print >>sys.stderr, "***"
                #mismatchfraction = mismatchcounts[i]/(1.*covcounts[i]+10)
                if True: #currcount > mismatchthreshold and allcovcount[i] > 20:
                    #print >>mismatchreport, currfeat.name+"\t"+currsample+"\t"+str(positionnums[i])+"\t"+str(currcount) 
                    print >>mismatchreport, "\t".join([currfeat.name,currsample,str(positionnums[i]),str(currcount),str(allcovcount[i]),str(allstarts[i]),str(1.*readcounts[currsample][currfeat.name]),trnastk.aligns[currfeat.name][i],str(mismatches[i]),str(adeninecount[i]),str(thyminecount[i]),str(cytosinecount[i]),str(guanosinecount[i]) ])
            
    mismatchfile = open(mismatchfilename, "w")
    #print >>allcoveragefile, "Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
    print >>mismatchfile, "Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)


    for currfeat in trnalist:
        totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
        if totalreads < mincoverage:
            continue
        if combinereps:
    
            replicates = sampledata.allreplicates()
            for currrep in replicates:
                covcounts = list(curr/(1.*readcounts[currrep][currfeat.name] + 20.0) if curr is not None else 0 for curr in sumsamples(readmismatches,sampledata,currrep,currfeat))
                #print >>mismatchfile, currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in covcounts)
                #print >>allcoveragefile, currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr/(1.*covcounts[i]+10)) if curr is not None else "NA" for i, curr in enumerate(sumsamples(allcoverages,sampledata,currrep,currfeat)))

        else:
            for currsample in samples:
                allcovcount  = list(curr if curr is not None else 0 for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
                covcounts  = list(curr/(1.*allcovcount[i] + 20.0) if curr is not None else 0 for i, curr in enumerate(readmismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name])))

                #covcounts = list(curr/(1.*readcounts[currsample][currfeat.name]) if curr is not None else "NA" for curr in readmismatches[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
                print >>mismatchfile, currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in covcounts)
                #print >>allcoveragefile, currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr/(1.*covcounts[i]+10)) if curr is not None else "NA" for i, curr in enumerate(allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name])))

                
                        
    if uniquename:
        #multaminocoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
        #multaccoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
        #multtrnacoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
        #uniquecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
        covfiles = {uniquename + '-uniqueends.txt':uniquecoverages,uniquename + '-multacends.txt':multaccoverages,uniquename + '-multtrnaends.txt':multtrnacoverages,uniquename + '-multaminoends.txt':multaminocoverages}
        
        for filename, currcoverage in covfiles.iteritems():
            #print >>sys.stderr, "file"
            #print >>sys.stderr, filename
            #sys.exit()
            
            covfile = open(filename, "w")
            print >>covfile,"Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
                
            for currfeat in trnalist:
                totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
                if totalreads < mincoverage:
                    continue
                if combinereps:
            
                    replicates = sampledata.allreplicates()
                    for currrep in replicates:
                        print  >>covfile,currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr/(1.*readcounts[currrep][currfeat.name]))  if curr is not None else "NA" for curr in sumsamples(currcoverage,sampledata,currrep,currfeat))
                else:
                    for currsample in samples:
                        print  >>covfile,currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr/(1.*readcounts[currsample][currfeat.name]))  if curr is not None else "NA" for curr in currcoverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
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
