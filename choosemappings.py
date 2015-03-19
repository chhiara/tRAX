#!/usr/bin/env python

import re
import sys
import os.path
import itertools

from sequtils import *
from collections import defaultdict



readtargets = set()

hitscores = dict()
curreadname = None
readlength = None
totalmatch = re.compile(r"^(?:(?P<startclip>\d+)S)?(?P<matchlength>\d+)M(?:(?P<endclip>\d+)S)?$")
#19S17M
totalclips = 0
totalreads = 0
linkerseqs = linkerseqs | set(revcom(curr) for curr in linkerseqs)
multimaps = 0
shortened = 0
mapsremoved = 0
totalmaps = 0
#prints the sam header that samtools need to convert to bam
print getsamheader(sys.argv[1]).rstrip()
#grab and iterate through all mappings of a single read
#most mappers will ensure that read mappings are in the same order as they were in the fastq file, with all mappings of the same read together
#once the file has been sorted using "samtools sort", this doesn't work anymore.  I can't detect that here, so no error message will be output
trnareads = 0
maxreads = 0
diffreads = 0
for pairedname, allmaps in itertools.groupby(readsam(sys.argv[1], getseqs = True, getline = True),lambda x: x.name+"#"+str(x.data["firstpair"])+"/"+str(x.data["secondpair"])):
    #print >>sys.stderr, "***"
    totalreads += 1
    if totalreads > 100000:
        pass
    
    readlength = None
    hitscores = dict()
    readtargets = set()
    clipsize = 50
    mappings = 0
    currscore = None
    newset = set()
    readlength = None
    #iterate through all mappings of the current read
    for currmap in allmaps:
        totalmaps += 1
        if currmap.chrom == "*":
            continue
        readlength = len(currmap.data["seq"])
        mappings += 1
        #if this is the best mapping of this read, discard any previous mappings
        if currscore is None or currscore < currmap.data["tags"]["AS"]:
            newset = set()
            newset.add(currmap)
            currscore = currmap.data["tags"]["AS"]
        #if this mappings is the same as the previous best mapping, add it to the list
        elif currscore == currmap.data["tags"]["AS"]:
            newset.add(currmap)
        else:
            pass
    #here is where I count a bunch of thing so I can report them at the end
    if mappings > 1:
        multimaps += 1
    if len(newset) < mappings:
        #print  >>sys.stderr, pairedname
        #print >>sys.stderr, str(len(newset))+"/"+str(mappings)
        mapsremoved += mappings - len(newset)
        shortened += 1
    #print str(len(newset))+"\t"+str(readlength)
    #best mappings are printed out here
    if len(newset) >= 50:
        maxreads += 1
        #print >>sys.stderr, len(newset)
    if sum("tRNA" in curr.chrom for curr in newset) > 0:
        trnareads += 1
        diff = len(newset) - sum("tRNA" in curr.chrom for curr in newset)
        if diff > 0:
            diffreads += 1
        for curr in newset:
            if "tRNA" in curr.chrom:
                pass
                print curr.data["bamline"].rstrip()
                if curr.start < 20:
                    pass
                    #print curr.data["bamline"].rstrip()
            else:
                #print curr.data["bamline"].rstrip()
                pass
    else:
        for curr in newset:
            pass
            print curr.data["bamline"].rstrip()
        
print >>sys.stderr, str(diffreads)+"/"+str(trnareads)
print >>sys.stderr, str(trnareads)+"/"+str(totalreads)
print >>sys.stderr, str(maxreads)+"/"+str(totalreads)
print >>sys.stderr, str(multimaps)+"/"+str(totalreads)
print >>sys.stderr, str(shortened)+"/"+str(multimaps)
print >>sys.stderr, str(mapsremoved)+"/"+str(totalmaps)