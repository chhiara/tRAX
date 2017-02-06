#!/usr/bin/env python

import re
import os
import sys
import itertools

from collections import defaultdict
import argparse
from parsetrnas import *
from trnasequtils import *



#This program gets the mature tRNA sequences

def main(**args):
    args = defaultdict(lambda: None, args)
    #>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-2 (tRNAscan-SE ID: chrVI.trna6) chrVI:204924-204996 (-) Ala (AGC) 73 bp Sc: 69.6
                                #nmt-tRNA-Gln-TTG-13-1
    trnanamere = re.compile(r"^\>\S+_((?:\w+\-)?tRNA\-\w+\-[\w\?]+\-\d+\-\d+)\s+\(tRNAscan\-SE\s+ID:\s+(\S+)\)")
    #trnanamere = re.compile(r"^\>\S+_(tRNA\-\w+\-\w+\-\d+\-\d+)\s+")#\((S+)\)")
    if "maturetrnafa" in args:
        maturetrnafa = open(args["maturetrnafa"], "w")
    else:
        maturetrnafa=sys.stdout
    gtrnatrans = None
    if args["gtrnafa"]:
        gtrnafa = open(args["gtrnafa"], "r")
        gtrnatrans = dict()
        for currline in gtrnafa:
            trnamatch = trnanamere.match(currline)
            #print >>sys.stderr, currline
            if trnamatch:

                print >>sys.stderr, trnamatch.group(1)+":"+trnamatch.group(2)
                gtrnatrans[trnamatch.group(2)] = trnamatch.group(1)
            elif currline.startswith(">"):
                pass
                #print >>sys.stderr, currline
        if len(gtrnatrans.keys()) == 0:
            print >>sys.stderr, "Could not extract names from gtrnadb fasta file"
            print >>sys.stderr, "must have names in format '>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-10 (tRNAscan-SE ID: chrXIII.trna9)'"
            sys.exit()
            
            
    alltrnas = list()
    trnascantrnas = list()
    trnadbtrnas = list()
    trnacentraltrnas = list()
    for currfile in args["trnascan"]:
        if gtrnatrans:
            trnadbtrnas.extend(readtRNAdb(currfile, args["genome"], gtrnatrans))
        else:
            trnascantrnas.extend(readtRNAscan(currfile, args["genome"]))
            #print >>sys.stderr, len(trnascantrnas)
            
        
    '''
    for currfile in args["trnascan"]:
        trnascantrnas.extend(readtRNAscan(currfile, args["genome"]))
    trnacentraltrnas = list()
    for currfile in args["rnacentral"]:
        trnacentraltrnas.extend(readrnacentral(currfile,args.chromtranslate,mode = 'transcript'))
    '''    
    alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas + trnadbtrnas
    mitomode = args["mitomode"]
    trnabed = None
    if args["bedfile"]:
        trnabed = open(args["bedfile"], "w")
    
    trnatable = None
    if args["maturetrnatable"]:
        trnatable = open(args["maturetrnatable"], "w")
    
    
    
    
    def readmultistk(struct):
        currrecord = ""
        structs = list()
        for line in struct.split("\n"):
            currrecord += line+"\n"
            if line == "//":
                yield currrecord
                currrecord = ""
                
    
    margin = 20
    anticodoncount = defaultdict(int)
    trnanames = dict()
    trnalist = list()
    for currtrans in alltrnas:
        if currtrans.name is None:
            name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
            currtrans.name = name
            anticodoncount[currtrans.anticodon] += 1
            
            
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    
    #cp /projects/lowelab/users/pchan/data/tRNA/stdModels/infernal-1.1/TRNAinf.cm ./
    if mitomode:
        trnacmfile = scriptdir+'TRNAMatureMitoinf.cm'
    else:
        trnacmfile = scriptdir+'trnamature-euk.cm'
    stkfile = args["trnaalignment"]
    if args["trnaalignment"]:
        devnull = open(os.devnull, 'w')
        seqfile = tempmultifasta(((currtrans.name, currtrans.getmatureseq()) for currtrans in alltrnas))
        cmcommand = ['cmalign', "-o", stkfile,"--nonbanded", "-g",trnacmfile,seqfile.name]
        #print >>sys.stderr, " ".join(cmcommand)
        cmrun = subprocess.Popen(cmcommand, stdout = devnull)
        result = cmrun.wait()
        if result:
            print >>sys.stderr, "Failure to align tRNAs"
            sys.exit(1)
        #stkout = cmrun.communicate()[0]
        #trnaalign = readrnastk(stkout.split("\n"))[0]
        seqfile.close()
        
    if args["locibed"]:
        locibed = open(args["locibed"],"w")
        
    for currtrans in alltrnas:
        name = currtrans.name
        #trnanames[name] = currtrans
        trnalist.append(name)
        #print >>sys.stderr, name
        print >>maturetrnafa, ">"+name
        print >>maturetrnafa, str("N" * margin) +currtrans.getmatureseq()+str("N" * margin)
        if trnatable is not None:
            print >>trnatable, "\t".join([name,",".join(currlocus.name for currlocus in currtrans.loci),currtrans.amino,currtrans.anticodon])
        if trnabed is not None:
            transcriptrange = GenomeRange("genome", name, margin, margin + len(currtrans.getmatureseq()), strand = "+", name = name)
            print >>trnabed, transcriptrange.bedstring()
        if args["locibed"]:
            for currlocus in currtrans.loci:
                pass
                print >>locibed, currlocus.loc.bedstring()
        
    
    
    #sys.exit()
    trnamods = dict()
    allmods = set()
    nomatches = 0
    trnamismatches = dict()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--trnascan',  nargs='+', default=list(),
                       help='tRNAscan-SE file')
    parser.add_argument('--rnacentral', nargs='+', default=list(),
                       help='RNAcentral tRNA file')
    parser.add_argument('--genome', 
                       help='fasta sequence of genome')
    parser.add_argument('--chromtranslate', 
                       help='translation table of chromosome names')
    parser.add_argument('--maturetrnafa', 
                       help='output file for mature tRNA sequences')
    parser.add_argument('--bedfile', 
                       help='Output bedfile of coordinates for mature tRNA')
    parser.add_argument('--tag', nargs='1',
                       help='tag to be added to tRNA name')
    parser.add_argument('--maturetrnatable',
                       help='Output table of mature tRNAs')
    parser.add_argument('--locibed',
                       help='Output BED format trna loci')
    parser.add_argument('--trnaalignment',
                       help='Output stockholm format alignment of mature tRNAs')
    parser.add_argument('--gtrnafa', 
                       help='fasta sequence tRNAs from tRNAscan-SE database')
    parser.add_argument('--mitomode', action="store_true", default=False,
                       help='Use mitochondrial models for alignment')
    maturetrnafa
    
    args = vars(parser.parse_args())
    main(args)    

