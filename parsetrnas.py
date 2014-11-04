#!/usr/bin/env python

import re
import random
import math
import os
import subprocess
import tempfile
import sys
import getopt
import itertools
import codecs
import collections
import time
import string
import gzip
import os.path
from collections import defaultdict
from sequtils import *
from rnautils import *



class tRNAlocus:
    def __init__(self, loc, seq, score, amino, anticodon, intronseq):
        self.name = loc.name
        self.loc = loc
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.intronseq = intronseq

class tRNAtranscript:
    def __init__(self, seq, score, amino, anticodon, loci, intronseq, name = None):
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.loci = tuple(loci)
        self.intronseqs = intronseq
        self.name = name

def getuniquetRNAs(trnalist):
    sequencedict = defaultdict(list)
    scoredict = defaultdict(set)
    for curr in trnalist:
        sequencedict[curr.seq].append(curr)
    for currtrans in sequencedict.iterkeys():
        scores = set(curr.score for curr in sequencedict[currtrans])
        anticodon = set(curr.anticodon for curr in sequencedict[currtrans])
        amino = set(curr.amino for curr in sequencedict[currtrans])
        introns = set(curr.intronseq for curr in sequencedict[currtrans])
        #remove psuedo if there's a real one somewheres
        if len(anticodon) > 1:
            anticodon.discard('Xxx')
            amino.discard('X')
        if introns == set([""]):
            introns = set()
        loci = list(curr for curr in sequencedict[currtrans])
        if len(scores) > 1:
            #print >>sys.stderr, "Multiple scores"
            pass
        if len(anticodon) > 1:
            print >>sys.stderr, "Multiple ACs"
            break
        yield tRNAtranscript(currtrans, scores,list(amino)[0],list(anticodon)[0],loci, introns)
        
        

def readrnacentral(scanfile, mode = 'locus'):
    orgname = "genome"
   
    #print reftochrom
    trnafile = open(scanfile)
    transcriptinfo = defaultdict(list)
    for line in trnafile:
        fields = line.split(",")
        #print len(fields)
        #print "**"
        if fields[0] == 'Entry number':
            continue
        #print fields[1]
        if True:
            name = fields[2]
            amino = fields[11]
            anticodon = fields[12]
            sequence = fields[21]
            score = float(fields[18])
            #print >>sys.stderr, score
            #print name+":"+amino+":"+anticodon
            #print sequence
            chrom = re.sub(r'\.\d+$', '', fields[7])
            start = int(fields[8].split('-')[0]) - 1
            end = fields[8].split('-')[1]
            transcriptname = re.sub(r'-\d+$', '', fields[2])
            if fields[9] == "yes":
                strand = "-"
                
            elif fields[9] == "no":
                strand = "+"
            currtRNA = GenomeRange(orgname, chrom,start,end, name = name,strand = strand,orderstrand = True)
            intronnums = set()
            intronseqs = ""
            for currintron in fields[14:16]:
                intronmatch  = re.search(r'(\d+)\.\.(\d+)',currintron)
                if intronmatch:
                    intronstart = int(intronmatch.group(1))-1
                    intronend = int(intronmatch.group(2))
                    intronnums |= set(range(intronstart, intronend))
                    intronseqs = sequence[intronstart:intronend]
            
            newseq = ''
            for i in range(len(sequence)):
                if i in set(intronnums):
                    newseq += '-'
                else:
                    newseq += sequence[i]
                    
            finalseq = newseq.replace('-','') + "CCA"
            if amino == "His":
                finalseq = "G" +finalseq 
                #print >>sys.stderr, "****"
            currlocus = tRNAlocus(currtRNA,finalseq, score,amino,anticodon,intronseqs)
            transcriptinfo[transcriptname].append(currlocus)
            if mode == 'locus':
                yield currlocus
        allseqs = dict()
    print >>sys.stderr, len(transcriptinfo.keys())
    for currtrans in transcriptinfo.iterkeys():
        if len(set(curr.seq for curr in transcriptinfo[currtrans])) > 1:
            print >>sys.stderr, "multiple"
        #print transcriptinfo[currtrans][0].seq
        if transcriptinfo[currtrans][0].seq in allseqs:
            print  >>sys.stderr, "duplicate:" + currtrans + ":"+allseqs[transcriptinfo[currtrans][0].seq]
        allseqs[transcriptinfo[currtrans][0].seq] = currtrans
        if mode == 'transcript':
            #print >>sys.stderr,currtrans 
            yield tRNAtranscript( transcriptinfo[currtrans][0].seq, set(curr.score for curr in transcriptinfo[currtrans]), transcriptinfo[currtrans][0].amino, transcriptinfo[currtrans][0].anticodon, set(transcriptinfo[currtrans]), transcriptinfo[currtrans][0].intronseq, name = currtrans)
        

def readtRNAscan(scanfile, genomefile, mode = None):
    #mode = 'gtRNAdb'
    trnalist = list()
    orgname = "genome"
    if  hasattr(scanfile ,'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile)
    trnascore = dict()
    trnaanticodon = dict()
    trnaamino = dict()
    tRNAintron = dict()
    trnas = dict()
    for currline in trnascan:
        if not currline.startswith("chr"):
            continue
            pass
        fields = currline.split()
        
        if mode == "gtRNAdb":
            #print >>sys.stderr, fields[6:8]
            del fields[6:8]
        if  (fields[4] == "Pseudo"):
            #fields[4] = "Pseudo"
            if fields[5] == "???":
                fields[5] = 'Xxx'
            else:
                continue
        if fields[0] == 'chrM' and not includemito:
            continue
        if fields[2] > fields[3]:
            fields[3] = int(fields[3]) - 1
            
        else:
            fields[3] = int(fields[3])
            fields[2] = int(fields[2]) - 1
        curramino = fields[4]
        currac = fields[5]
        
        currtRNA = GenomeRange(orgname, fields[0],fields[2],fields[3], name = fields[0]+"."+"tRNA"+fields[1]+"-"+fields[4]+fields[5],strand = "+",orderstrand = True)
        currtrans = currtRNA

        trnaamino[currtrans.name] = fields[4]
        trnaanticodon[currtrans.name] = fields[5]
        trnascore[currtrans.name] =  float(fields[8])
        trnas[currtrans.name] =  currtrans
    
    
        
        currtRNA.fastafile = genomefile
        trnalist.append(currtRNA)
        if int(fields[6]) != 0:
            if currtRNA.strand ==  "-":
                intronstart = int(fields[2]) - int(fields[6]) 
                intronend = int(fields[2]) - int(fields[7]) +1
            else:
                intronstart = int(fields[6]) - int(fields[2]) - 1
                intronend = int(fields[7]) - int(fields[2])
            tRNAintron[currtRNA.name] = tuple([intronstart, intronend])
    trnaseqs = getseqdict(trnalist, faifiles = {orgname:genomefile+".fai"})
    intronseq = defaultdict(str)
    for curr in trnaseqs.iterkeys():
        if curr in tRNAintron:
            start = tRNAintron[curr][0]
            end = tRNAintron[curr][1]
            intronseq[curr] = trnaseqs[curr][start:end]
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
        trnaseqs[curr] += "CCA"
        #print >>sys.stderr,trnaamino[currtrans.name]
        if trnaamino[curr] == "His":
            trnaseqs[curr] = "G"+trnaseqs[curr]
            #print >>sys.stderr, "***"
        
        yield tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr],trnaamino[curr],trnaanticodon[curr],intronseq[curr])

