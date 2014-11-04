#!/usr/bin/env python

import re
import random
import math
import os
import subprocess
import tempfile
import sys
import itertools
import configfile
from rnautils import *
from orthologtools import *
from sequtils import *
from collections import defaultdict
import argparse
from parsetrnas import *


#This program gets the mature tRNA sequences

parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--trnascan',  nargs='+', default=list(),
                   help='tRNAscan-SE file')
parser.add_argument('--rnacentral', nargs='+', default=list(),
                   help='RNAcentral tRNA file')
parser.add_argument('--genome', 
                   help='fasta sequence of genome')
parser.add_argument('--chromtranslate', 
                   help='translation table of chromosome names')

parser.add_argument('--tag', nargs='1',
                   help='tag to be added to tRNA name')

args = parser.parse_args()


#/projects/lowelab/users/holmes/pythonsource/trnastuff/ucscToINSDC.txt

'''
./getmaturetrnas.py --rnacentral hg19-trnacentral.txt  --genome /scratch/encodeshortrna/hg19.fa --chromtranslate /projects/lowelab/users/holmes/pythonsource/trnastuff/ucscToINSDC.txt
'''
#These two things are artifacts of the 
#centralspecies = 'Hsapi19/Hsapi19-summary'
#chromnames = {'Hsapi19/Hsapi19-summary' : '/projects/lowelab/users/holmes/pythonsource/trnastuff/ucscToINSDC.txt', 'Scere3/Scere3-summary' :'/projects/lowelab/users/holmes/pythonsource/trnastuff/sacCer3-gentochrom.txt'}




alltrnas = list()
trnascantrnas = list()
for currfile in args.trnascan:
    trnascantrnas.extend(readtRNAscan(currfile, args.genome, mode = 'gtRNAdb'))
trnacentraltrnas = list()
for currfile in args.rnacentral:
    trnacentraltrnas.extend(readrnacentral(currfile,mode = 'transcript'))
    
alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas
                


anticodoncount = defaultdict(int)
trnanames = dict()
trnalist = list()
for currtrans in alltrnas:
    if currtrans.name is None:
        name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
        anticodoncount[currtrans.anticodon] += 1
    else:
        #print >>sys.stderr, "***"
        name = currtrans.name
    #trnanames[name] = currtrans
    trnalist.append(name)
    #print >>sys.stderr, name
    print ">"+name
    print currtrans.seq
#sys.exit()
trnamods = dict()
allmods = set()
nomatches = 0
trnamismatches = dict()

