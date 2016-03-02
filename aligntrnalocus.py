#!/usr/bin/env python

import pysam
import os
import sys
import argparse
import subprocess
from collections import defaultdict
from trnasequtils import *









'''
~/pythonsource/aligntrnalocus.py --genomefile= --stkfile= --trnaloci

~/pythonsource/trnaseq/aligntrnalocus.py --genomefile=hg19.fa --stkfile=hg19-trnaloci.stk --trnaloci=hg19-trnaloci.bed
'''

def main(**args):
    args = defaultdict(lambda: None, args)
    stkfile = args["stkfile"]
    genomefile = args["genomefile"]
    mitomode = args["mitomode"]
    trnalocifile = args["trnaloci"]
    
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"        
    
    if mitomode:
        trnacmfile = scriptdir+'TRNAinf.cm'
    else:
        trnacmfile = scriptdir+'TRNAinf-euk.cm'
    
    trnaloci = list(readbed(trnalocifile, orgdb = "genome", seqfile=genomefile))
    lociseqs = getseqdict(trnaloci, faifiles = {"genome":genomefile+".fai"})
    #print lociseqs
    #lociseqfile = tempmultifasta(lociseqs)
    devnull = open(os.devnull, 'w')
    seqfile = tempmultifasta(lociseqs.iteritems())
    cmcommand = ['cmalign', "-o", stkfile,"--nonbanded", "-g",trnacmfile,seqfile.name]
    #print >>sys.stderr, " ".join(cmcommand)
    cmrun = subprocess.Popen(cmcommand, stdout = devnull)
    result = cmrun.wait()
    if result:
        print >>sys.stderr, "Failure to align tRNAs"
        sys.exit(1)
    seqfile.close()
    devnull.close()
    
    #trnaalign
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--trnaloci',
                       help='bed file of tRNA loci')      
    parser.add_argument('--genomefile',
                       help='fasta file of genome')
    parser.add_argument('--stkfile',
                       help='stockholm output file')
    parser.add_argument('--mitomode', action="store_true", default=False,
                       help='Use mitochondrial models')
    args = vars(parser.parse_args())
    main(args)    

        