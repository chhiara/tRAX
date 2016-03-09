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




parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--databasename',
                   help='database name to be used')
parser.add_argument('--genomefile',
                   help='Fasta file containing genome sequence')
parser.add_argument('--trnascanfile',
                   help='output from tRNAscan-SE run')
parser.add_argument('--gtrnafafile',
                   help='Fasta file of tRNA sequences from gtRNAdb')



args = parser.parse_args()
dbname = args.databasename
scanfile = args.trnascanfile
genomefile = args.genomefile
gtrnafafile = args.gtrnafafile

#$1 is database name
#trnascanfile is trnascan file
#genomefile is fasta file of genome



#subprocess.call("samtools faidx "+genomefile, shell=True)

getmaturetrnas.main(trnascan=[scanfile], genome=genomefile,gtrnafa=gtrnafafile,bedfile=dbname+"-maturetRNAs.bed",maturetrnatable=dbname+"-trnatable.txt",trnaalignment=dbname+"-trnaalign.stk",locibed=dbname+"-trnaloci.bed",maturetrnafa=dbname+"-maturetRNAs.fa")
aligntrnalocus.main(genomefile=genomefile,stkfile=dbname+"-trnaloci.stk",trnaloci=dbname+"-trnaloci.bed")


subprocess.call("cat "+dbname+"-maturetRNAs.fa "+genomefile+" >"+dbname+"-tRNAgenome.fa", shell=True)
subprocess.call("samtools faidx "+dbname+"-tRNAgenome.fa", shell=True)
subprocess.call("bowtie2-build "+dbname+"-tRNAgenome.fa "+dbname+"-tRNAgenome", shell=True)

