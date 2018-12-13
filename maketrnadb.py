#!/usr/bin/env python

import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import getmaturetrnas
import aligntrnalocus
from distutils.spawn import find_executable
from distutils.version import LooseVersion, StrictVersion
import time

def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
        


parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--databasename',required=True,
                   help='database name to be used')
parser.add_argument('--genomefile',required=True,
                   help='Fasta file containing genome sequence')
parser.add_argument('--trnascanfile',required=True,
                   help='output from tRNAscan-SE run')
parser.add_argument('--gtrnafafile',
                   help='Fasta file of tRNA sequences from gtRNAdb')
parser.add_argument('--namemapfile',
                   help='Name mapping from gtRNAdb')


def shellcall(shellcommand,failquit = False):
    retcode = subprocess.call(shellcommand, shell=True)
    if retcode > 0 and failquit:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
        print >>sys.stderr, "quiting program..."
        sys.exit(1)
    elif retcode > 0:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
    return retcode
        
        
args = parser.parse_args()
dbname = args.databasename
scanfile = args.trnascanfile
genomefile = args.genomefile
gtrnafafile = args.gtrnafafile
namemapfile = args.namemapfile



runtime = time.time()
loctime = time.localtime(runtime)


#$1 is database name
#trnascanfile is trnascan file
#genomefile is fasta file of genome

#test command line programs


    
get_location("samtools")
get_location("bowtie2-build")

if not os.path.isfile(genomefile+".fai"):
    shellcall("samtools faidx "+genomefile)
    
    


getmaturetrnas.main(trnascan=[scanfile], genome=genomefile,gtrnafa=gtrnafafile,namemap=namemapfile, bedfile=dbname+"-maturetRNAs.bed",maturetrnatable=dbname+"-trnatable.txt",trnaalignment=dbname+"-trnaalign.stk",locibed=dbname+"-trnaloci.bed",maturetrnafa=dbname+"-maturetRNAs.fa")
aligntrnalocus.main(genomefile=genomefile,stkfile=dbname+"-trnaloci.stk",trnaloci=dbname+"-trnaloci.bed")


shellcall("cat "+dbname+"-maturetRNAs.fa "+genomefile+" >"+dbname+"-tRNAgenome.fa", failquit = True)
    
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
gitversion, githash = getgithash(scriptdir)
    
    
dbinfo = open(dbname+ "-dbinfo.txt","w")
print >>dbinfo, "time\t"+str(runtime)+"("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "genomefile\t"+str(genomefile)
print >>dbinfo, "trnascanfile\t"+str(scanfile)
print >>dbinfo, "git version\t"+str(gitversion)
dbinfo.close()
    
shellcall("bowtie2-build "+dbname+"-tRNAgenome.fa "+dbname+"-tRNAgenome", failquit = True)




