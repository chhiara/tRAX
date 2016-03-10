#!/usr/bin/env python

import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import os

import mapreads
import countreads
import getcoverage
import countreadtypes



#expname is experiment name
#dbname is database name
#samplefile is sample file
#$4 is bed feature for other sRNAs




parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--experimentname',
                   help='experiment name to be used')
parser.add_argument('--databasename',
                   help='name of the tRNA database')
parser.add_argument('--samplefile',
                   help='sample file')
parser.add_argument('--ensemblgtf',
                   help='The ensembl gene list for that species')
parser.add_argument('--exppairs',
                   help='List of sample pairs to compare')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='Additional bed files for feature list')
parser.add_argument('--forceremap', action="store_true", default=False,
                   help='Force a remap of the reads, even if they are already mapped')





args = parser.parse_args()
dbname = args.databasename
expname = args.experimentname
pairfile =  args.exppairs
ensgtf = args.ensemblgtf
samplefile = args.samplefile
forceremap = args.forceremap
bedfiles= args.bedfile
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"




#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)

mapreads.main(samplefile=samplefile, trnafile=dbname+"-trnatable.txt",bowtiedb=dbname+"-tRNAgenome",logfile=expname+"/"+expname+"-mapstats.txt", force=forceremap)
countreads.main(samplefile=samplefile,ensemblgtf=ensgtf,maturetrnas=[dbname+"-maturetRNAs.bed"],trnaloci=[dbname+"-trnaloci.bed"],removepseudo=True,genetypefile=expname+"/"+expname+"-genetypes.txt" ,trnatable=dbname+"-trnatable.txt",countfile=expname+"/"+expname+"-counts.txt",bedfile=bedfiles)


if pairfile:

	subprocess.call("Rscript "+scriptdir+"/analyzecounts.R "+expname +" "+expname+"/"+expname+"-counts.txt "+samplefile+" "+ pairfile, shell=True)
	subprocess.call("Rscript "+scriptdir+"/makescatter.R "+expname +" "+expname+"/"+expname+"-normalized.txt "+ dbname+"-trnatable.txt "+expname+"/"+expname+"-genetypes.txt "+samplefile+" "+pairfile, shell=True)
else:
	
	subprocess.call("Rscript "+scriptdir+"/analyzecounts.R "+expname +" "+expname+"/"+expname+"-counts.txt "+samplefile, shell=True)
	pass

        #"$SCRIPTDIR/countreadtypes.py" --sizefactors=expname/expname-SizeFactors.txt --combinereps --samplefile=samplefile  --maturetrnas=dbname-maturetRNAs.bed --trnatable=dbname-trnatable.txt --trnaaminofile=expname/expname-aminocounts.txt --ensemblgtf $4 --trnaloci=dbname-trnaloci.bed   >expname/expname-typecounts.txt #--countfrags
#sys.exit()
countreadtypes.main(sizefactors=expname+"/"+expname+"-SizeFactors.txt",combinereps= True ,samplefile=samplefile,maturetrnas=[dbname+"-maturetRNAs.bed"],trnatable=dbname+"-trnatable.txt",trnaaminofile=expname+"/"+expname+"-aminocounts.txt",ensemblgtf=ensgtf,trnaloci=[dbname+"-trnaloci.bed"],countfile=expname+"/"+expname+"-typecounts.txt",bedfile= bedfiles)
subprocess.call("Rscript "+scriptdir+"/featuretypes.R "+expname+"/"+expname+"-typecounts.txt "+expname+"/"+expname+"-typecounts.pdf", shell=True)
subprocess.call("Rscript "+scriptdir+"/featuretypes.R "+expname+"/"+expname+"-aminocounts.txt "+expname+"/"+expname+"-aminocounts.pdf",shell=True)




#sys.exit()
#countfragsize.main(combinereps=True, samplefile=samplefile,maturetrnas=dbname+"-maturetRNAs.bed",countfrags=True,trnaloci=dbname+"-trnaloci.bed",ensemblgtf=  --trnanormfile=expname-tRNANormFactors.txt --allreadsnormfile=expname-AllreadNormFactors.txt >expname/expname-readlengths.txt
#subprocess.call("Rscript "+scriptdir+"/readlengthhistogram.R "+expname+"/"+expname+"-readlengths.txt "+expname+"/"+expname+"-readlengths.pdf", shell=True)
	
#Rscript "$SCRIPTDIR/coverageplots.R" --cov=expname/expname-coverage.txt --trna=dbname-trnatable.txt --samples=samplefile --allcov=expname/expname-coverage.pdf  --uniquename=expname/expname --modomics=dbname-modomics.txt --multicov=expname/expname-multipagecoverage.pdf --combinecov=expname/expname-combinecoverage.pdf --directory=expname


getcoverage.main(samplefile=samplefile,bedfile=[dbname+"-maturetRNAs.bed"],sizefactors=expname+"/"+expname+"-SizeFactors.txt",stkfile=dbname+"-trnaalign.stk",uniquename=expname+"/"+expname, allcoverage=expname+"/"+expname+"-coverage.txt")


subprocess.call("Rscript "+scriptdir+"/coverageplots.R --cov="+expname+"/"+expname+"-coverage.txt --trna="+dbname+"-trnatable.txt --samples="+samplefile+" --allcov="+expname+"/"+expname+"-coverage.pdf  --uniquename="+expname+"/"+expname+" --modomics="+dbname+"-modomics.txt --multicov="+expname+"/"+expname+"-multipagecoverage.pdf --combinecov="+expname+"/"+expname+"-combinecoverage.pdf --directory="+expname, shell=True)
getcoverage.main(samplefile=samplefile ,bedfile=[dbname+"-trnaloci.bed"], sizefactors=expname+"/"+expname+"-SizeFactors.txt",stkfile=dbname+"-trnaloci.stk",edgemargin=30, outfile = expname+"/"+expname+"-locicoverage.txt", allcoverage=expname+"/"+expname+"-locicoverage.txt")
subprocess.call("Rscript "+scriptdir+"/locuscoverage.R --cov="+expname+"/"+expname+"-locicoverage.txt --trna="+dbname+"-trnatable.txt --samples="+samplefile+" --allcov="+expname+"/"+expname+"-locicoverage.pdf --multicov="+expname+"/"+expname+"-locimultipagecoverage.pdf --combinecov="+expname+"/"+expname+"-locicombinecoverage.pdf", shell=True) 
