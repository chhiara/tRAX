#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *
import subprocess
import tempfile

from distutils.spawn import find_executable



def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
        
        
def convertbam(dbname,inputbam, outputbam, scriptdir, force = False, logfile = sys.stderr):
    if not os.path.isfile(outputbam) or force:
   
        #print >>sys.stderr, tempfile.gettempprefix()
        tempprefix = inputbam.split(".")[0]
        bamconvertcommand =scriptdir+'/convertbam.py '+inputbam+' '+dbname+' | samtools sort -T '+tempfile.gettempdir()+'/convert'+tempprefix+' - -o '+outputbam
        #print >>sys.stderr, bamconvertcommand
        #sys.exit(1)
        if logfile:
            print >>logfile,  bamconvertcommand
            pass
        bamconvertrun = None
        logfile.flush()
        bamconvertrun = subprocess.Popen(bamconvertcommand, shell = True, stderr = subprocess.PIPE)
        
        output = bamconvertrun.communicate()
        errinfo = output[1]
        if logfile is not None:
            pass
            print >>logfile, errinfo 
        logfile.flush()
        if bamconvertrun.returncode:
            print >>sys.stderr, "Failure to convert bam to genome space"
            print >>sys.stderr, "check logfile"
            logfile.close()
            sys.exit(1)

def samtoolsmerge(bamfiles, outbam, force = False):
    samtoolsloc = get_location("samtools")
    if not os.path.isfile(outbam) or force:
        samcommand = [samtoolsloc,"merge","-f", outbam]
        samcommand.extend(bamfiles)
        
        samtoolsjob = subprocess.Popen(samcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
        print >>sys.stderr, " ".join(samcommand)
        samtoolsresults = samtoolsjob.communicate()[0]
        print >>sys.stderr, samtoolsresults
    
    
def createtrackdb(allreps, expname):

    trackdb = open (expname+"/trackhub/trackdb.txt", "w")
    currpriority = 2.3
    
    print  >>trackdb, "track "+expname+"tracks            "
    print  >>trackdb, "compositeTrack on                     "
    print  >>trackdb, "shortLabel Data from "+expname+"   "
    print  >>trackdb, "longLabel Data from "+expname+"    "
    print  >>trackdb, "type bigWig                           "
    print  >>trackdb, "dragAndDrop on                        "
    print  >>trackdb, "autoScale on                          "
    print  >>trackdb, "alwaysZero on                         "
    print  >>trackdb, "maxHeightPixels 100:32:8              "
    print  >>trackdb, "\n"
    
    for currrep in allreps:
        

        print >>trackdb, "track "+currrep+"plus                                             "
        print  >>trackdb, "parent "+expname+"tracks                                   "
        print  >>trackdb, "bigDataUrl "+currrep+".Plus.bw                                 "
        print  >>trackdb, "shortLabel Plus "+currrep+"                                    "
        print  >>trackdb, "longLabel Plus strand coverage "+currrep+" all mapped reads    "
        print  >>trackdb, "color 220,148,44                                             "
        print  >>trackdb, "type bigWig                                                  "
        print  >>trackdb, "priority "+str(currpriority)+"                                       "
        print  >>trackdb, "\n"
        
        
        print  >>trackdb, "track "+currrep+"Minus                                             "
        print  >>trackdb, "parent "+expname+"tracks                                    "
        print  >>trackdb, "bigDataUrl "+currrep+".Minus.bw                                 "
        print  >>trackdb, "shortLabel Minus "+currrep+"                                    "
        print  >>trackdb, "longLabel Minus strand coverage "+currrep+" all mapped reads    "
        print  >>trackdb, "color 112,73,18                                       "
        print  >>trackdb, "type bigWig                                                  "
        print  >>trackdb, "priority "+str(currpriority + .1)+"                                  "
        print  >>trackdb, "\n\n\n"

        currpriority += .2
def makebigwigs(bamfile, repname, faifile, directory):
    #print >>sys.stderr, 'zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 '+bamfile+' | /projects/lowelab/users/holmes/bedtools/BEDTools/bin/genomeCoverageBed -bg -ibam stdin -g '+faifile+') '+faifile+' '+directory+"/"+repname+'.Plus.bw"'
    plusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 '+bamfile+' | /projects/lowelab/users/holmes/bedtools/BEDTools/bin/genomeCoverageBed -bg -ibam stdin -g '+faifile+') '+faifile+' '+directory+"/"+repname+'.Plus.bw"', shell = True)
    minusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -f 0x10 '+bamfile+' | /projects/lowelab/users/holmes/bedtools/BEDTools/bin/genomeCoverageBed -bg -ibam stdin -g '+faifile+') '+faifile+' '+directory+"/"+repname+'.Minus.bw"', shell = True)
    plusjob.wait()
    minusjob.wait()
    pass
	
def main(**args):
    dbname = args["genomedatabase"]
    samplefilename = args["samplefile"]
    sampledata = samplefile(args["samplefile"])
    expname = args["expname"]
    trackdir = expname+"/trackhub"
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"

    if not os.path.exists(trackdir):
        os.makedirs(trackdir)
    allsamples = sampledata.getsamples()
    for currsample in allsamples:
        currbam = sampledata.getbam(currsample)
        genomebam = currsample+"-genome.bam"
        convertbam(dbname, currbam, genomebam, scriptdir, force = True)
    
    faidxjob = subprocess.Popen("samtools faidx "+dbname+"-tRNAgenome.fa",shell = True)
    faidxjob.wait()
    for currrep in sampledata.allreplicates():
        repsamples = sampledata.getrepsamples(currrep)
        samtoolsmerge(list(curr+"-genome.bam" for curr in repsamples), currrep+"-mergegenome.bam", True)
        pysam.index(currrep+"-mergegenome.bam")
        makebigwigs(currrep+"-mergegenome.bam", currrep, dbname+"-tRNAgenome.fa.fai",trackdir)
    
    createtrackdb(sampledata.allreplicates(),expname)
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert TRAX bam file into genome bam file.')
    parser.add_argument('--genomedatabase', 
                       help='fasta sequence of genome')
    parser.add_argument('--samplefile', 
                       help='sample file')
    parser.add_argument('--expname', 
                       help='experiment name')

    
    args = vars(parser.parse_args())
    main(args)    