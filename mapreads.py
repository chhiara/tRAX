#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
from trnasequtils import *

MAXMAPS = 100

'''
~/pythonsource/trnaseq/mapreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bowtiedb=hg19-trnagenome >counts.txt

~/pythonsource/trnaseq/mapreads.py --samplefile=agingshort.txt --trnafile=sacCer3-trnatable.txt --bowtiedb=sacCer3-trnagenome >aging-counts.txt

cat <(samtools view -H HBV10_1.bam) <(samtools view HBV10_1.bam | grep UNC11-SN627:301:C21RKACXX:3:1101:20489:57142) | ~/pythonsource/trnaseq/choosemappings.py hg19-trnatable.txt | samtools view -

'''

numcores = 4



def wrapbowtie2(bowtiedb, unpaired, outfile, scriptdir, trnafile, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None):
    '''
    I think the quals are irrelevant, and N should be scored as only slightly better than a mismatch
    this does both
    ignoring quals forces the mismatches to have a single score, which is necessary for precise N weighting
    
     --ignore-quals --np 5
     Very sensitive is necessary due to mismatch caused by modified base misreads
    '''

    bowtiecommand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired
    if logfile:
        print >>logfile,  bowtiecommand
    #print >>sys.stderr, bowtiecommand
    logfile.flush()
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    bowtierun = None
    if logfile is not None:
        #print >>sys.stderr, "***LOG"
        bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = logfile)
    else:                                                       
        bowtierun = subprocess.Popen(bowtiecommand, shell = True)
        
    #bowtierun = subprocess.popen(bowtiecommand, shell = True, stderr = subprocess.PIPE)
    #errorcode = bowtierun.wait()
    #bowtierun.stdout.close()
    errorcode = bowtierun.wait()
    
    
    if errorcode:
        print >>sys.stderr, "Failure to Bowtie2 map"

            
        sys.exit(1)
    
'''
./trnaseq/mapreads.py --samplefile=HumanTrnas.txt --bowtiedb=/scratch/encodeshortrna/hg19
nohup ../trnaseq/mapreads.py --samplefile=SaraShortSamples.txt --bowtiedb=/scratch/encodeshortrna/hg19



nohup ../trnaseq/mapreads.py --samplefile=SaraSamples.txt --bowtiedb=/projects/lowelab/users/holmes/pythonsource/seqqa/combinedb/hg19-pad
'''

#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)
def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    sampledata = samplefile(argdict["samplefile"])
    trnafile = argdict["trnafile"]
    logfile = argdict["logfile"]
    bowtiedb = argdict["bowtiedb"]
    forcecreate = argdict["force"]
    #sys.exit()
    
    workingdir = './'
    #samplefile = open(args.samplefile)
    
    samples = sampledata.getsamples()
    
    #scriptdir = '/projects/lowelab/users/holmes/pythonsource/trnaseq/'
    trnafile = trnafile
    
    
        
        
        
    if logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr
    for samplename in samples:
        
        #put stuff in this so that it returns the err if bowtie2 fails instead of logging it
        
        #print >>sys.stderr, sampledata.getfastq(samplename)
        
        bamfile = workingdir+samplename
        #print >>sys.stderr, bamfile
        #sys.exit()
        #print >>sys.stderr, os.path.isfile(bamfile+".bam")
        if forcecreate or not os.path.isfile(bamfile+".bam"):
            print >>logfile, "Mapping "+samplename
            print >>sys.stderr, "Mapping "+samplename
            logfile.flush()
            wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile)
        
    
        
        
        #print >>logfile, "Processing "+samplename +" mappings"
    
        
        
        #result = subprocess.call(scriptdir+'choosemappings.py '+trnafile+' <'+bamfile +' | samtools view -F 4 -b - | samtools sort - '+workingdir+samplename+'_sort', shell = True)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')
    
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--trnafile',
                       help='tRNA file in format')
    parser.add_argument('--logfile',
                       help='optional log file for error messages and mapping stats')
    parser.add_argument('--bowtiedb',
                       help='Location of Bowtie 2 database')
    parser.add_argument('--force', action="store_true", default=False,
                       help='Force remapping even if mapping results exist')
    
    args = parser.parse_args()
    main(samplefile = args.samplefile, trnafile= args.trnafile, logfile = args.logfile, bowtiedb = args.bowtiedb, force = args.force)


