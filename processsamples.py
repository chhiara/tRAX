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
parser.add_argument('--experimentname',required=True,
                   help='experiment name to be used')
parser.add_argument('--databasename',required=True,
                   help='name of the tRNA database')
parser.add_argument('--samplefile',required=True,
                   help='sample file')
parser.add_argument('--ensemblgtf',
                   help='The ensembl gene list for that species')
parser.add_argument('--exppairs',
                   help='List of sample pairs to compare')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='Additional bed files for feature list')
parser.add_argument('--lazyremap', action="store_true", default=False,
                   help='Skip mapping reads if bam files exit')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='Omit fragment determination (Used for TGIRT mapping)')
parser.add_argument('--nosizefactors', action="store_true", default=False,
                   help='Don\'t use Deseq size factors in plotting')
parser.add_argument('--maxmismatch',
                   help='Maximum allowed mismatches')
rlogname = "Rlog.txt"
rlogfile = open(rlogname, "w")

def runrscript(*script):
    retcode = subprocess.call("Rscript "+" ".join(script), shell=True, stdout = rlogfile, stderr = subprocess.STDOUT)
    if retcode > 0:
        print >>sys.stderr, "R script "+script[0]+" failed"
        print >>sys.stderr, "Check "+rlogname+" for details"
        
        #sys.exit()
    return retcode
    

class trnadatabase:
    def __init__(self, dbname):
        self.trnatable = dbname+"-trnatable.txt"
        self.bowtiedb = dbname+"-tRNAgenome"
        self.locifile = dbname+"-trnaloci.bed"
        self.maturetrnas=dbname+"-maturetRNAs.bed"
        self.trnaalign = dbname+"-trnaalign.stk"
        self.locialign = dbname+"-trnaloci.stk"
        self.modomics = dbname+"-modomics.txt"
        
class expdatabase:
    def __init__(self, expname):
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt"
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"
        
        self.maplog = expname+"/"+expname+"-mapstats.txt"
        self.genetypes = expname+"/"+expname+"-genetypes.txt"
        self.genecounts = expname+"/"+expname+"-counts.txt"
        self.normalizedcounts = expname+"/"+expname+"-normalized.txt"
        self.sizefactors = expname+"/"+expname+"-SizeFactors.txt"

        self.genetypecounts=expname+"/"+expname+"-typecounts.txt"
        self.genetypeplot=expname+"/"+expname+"-typecounts.pdf"
        
        self.trnaaminofile=expname+"/"+expname+"-aminocounts.txt"
        self.trnaaminoplot=expname+"/"+expname+"-aminocounts.pdf"
        
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"

        self.locicoveragefile=expname+"/"+expname+"-locicoverage.txt"
        self.locicoverageplot=expname+"/"+expname+"-locicoverage.pdf"
        self.locicombinecoverageplot=expname+"/"+expname+"-locicombinecoverage.pdf"
        

def mapsamples(samplefile, trnainfo,expinfo, lazyremap):
    mapreads.main(samplefile=samplefile, trnafile=trnainfo.trnatable,bowtiedb=trnainfo.bowtiedb,logfile=expinfo.maplog,mapfile=expinfo.mapinfo, lazy=lazyremap)
def countfeatures(samplefile, trnainfo,expinfo, ensgtf, bedfiles):
    countreads.main(samplefile=samplefile,ensemblgtf=ensgtf,maturetrnas=[trnainfo.maturetrnas],trnaloci=[trnainfo.locifile],removepseudo=True,genetypefile=expinfo.genetypes,trnatable=trnainfo.trnatable,countfile=expinfo.genecounts,bedfile=bedfiles, nofrag=nofrag)
def counttypes(samplefile, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = False):
    if not ignoresizefactors:
        
        countreadtypes.main(sizefactors=expinfo.sizefactors,combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile )
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
    else:
        countreadtypes.main(combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)

def gettrnacoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],sizefactors=expinfo.sizefactors,stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allcoverage=expinfo.trnacoveragefile)
        runrscript(scriptdir+"/coverageplots.R","--cov="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allcoverage=expinfo.trnacoveragefile)
        runrscript(scriptdir+"/coverageplots.R","--cov="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--uniquename="+expname+"/"+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)

def getlocuscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],sizefactors=expinfo.sizefactors,stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile,minextend = 5)
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile,minextend = 5)
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)


args = parser.parse_args()
dbname = args.databasename
expname = args.experimentname
pairfile =  args.exppairs
ensgtf = args.ensemblgtf
samplefilename = args.samplefile
lazyremap = args.lazyremap
bedfiles= args.bedfile
nofrag= args.nofrag
nosizefactors = args.nosizefactors
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"


#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)
if not os.path.exists(expname+"/indiv"):
    os.makedirs(expname+"/indiv")

    
    
dbname = os.path.expanduser(dbname)
if ensgtf is not None:
    ensgtf = os.path.expanduser(ensgtf)
bedfiles = list(os.path.expanduser(curr) for curr in bedfiles)
trnainfo = trnadatabase(dbname)
expinfo = expdatabase(expname)
getsamples = samplefile(samplefilename)
if len(getsamples.getsamples()) == 1:
    nosizefactors = True

#Map the reads
print >>sys.stderr, "Mapping Reads"
mapsamples(samplefilename, trnainfo,expinfo, lazyremap)
#Count the reads for DEseq2 and scatter plots
print >>sys.stderr, "Counting Reads"
countfeatures(samplefilename, trnainfo,expinfo, ensgtf, bedfiles)
#Create a plot of mapped reads
print >>sys.stderr, "Analyzing counts"
runrscript(scriptdir+"/featuretypes.R",expinfo.mapinfo,expinfo.mapplot)

#Analyze counts and create scatter plots if pair file is provided
if pairfile:
	deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)
	if deseqret == 2:
	    print >>sys.stderr, "Deseq analysis failed, cannot continue"
	    sys.exit(1)
	runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
elif not nosizefactors:
	
	deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename)
	if deseqret == 2:
	    print >>sys.stderr, "Deseq analysis failed, cannot continue"
	    sys.exit(1)
#Count the reads by gene type
print >>sys.stderr, "Counting Read Types"
counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors)


#coverage plot of tRNAs
print >>sys.stderr, "Generating Read Coverage plots"      
gettrnacoverage(samplefilename, trnainfo,expinfo, ignoresizefactors = nosizefactors)

#coverage plot of pre-tRNAs
getlocuscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

