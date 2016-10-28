#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *




def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    countfrags = argdict["countfrags"]
    combinereps = argdict["combinereps"]
    ensemblgtf = argdict["ensemblgtf"]
    bamnofeature = argdict["bamnofeature"]
    trnatable = argdict["trnatable"]
    trnaaminofile = argdict["trnaaminofile"]
    sampledata = samplefile(argdict["samplefile"])
    minpretrnaextend = 5
    mitochrom = None
    if argdict["mitochrom"]:
        mitochrom = argdict["mitochrom"]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print >>sys.stderr, "Size factor file "+argdict["sizefactors"]+" missing "+currsample
                sys.exit(1)
        
    bedfiles = list()
    
    if argdict["bedfile"]  is not None:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if argdict["trnaloci"] is not None:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if argdict["maturetrnas"] is not None:
        maturetrnas = argdict["bedfile"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    
    genetypefile = argdict["genetypefile"]
    locifiles = argdict["trnaloci"]
    maturetrnafiles = argdict["maturetrnas"]    
    trnaaminofilename = argdict["trnaaminofile"]
    trnanormfile = argdict["trnanormfile"]
    allreadsnormfile = argdict["allreadsnormfile"]
    readlengthfile = argdict["readlengthfile"]
    
    if argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"],"w")
    

    
    
    wholetrnas = dict()
    fivefrags = dict()
    threefrags = dict()
    trailerfrags = dict()
    otherfrags = dict()
    allfrags = dict()
    alltrnas = list()
    
    
    ncrnaorder = defaultdict(int)
    #This decides the priority for ensembl genes, the rest I don't really care about yet
    for i, curr in enumerate(reversed(list(["snoRNA","miRNA", "rRNA","snRNA","misc_RNA","lincRNA", "protein_coding"]))):
        ncrnaorder[curr] = i + 1
    fullpretrnathreshold = 2
    
    trnainfo = transcriptfile(trnatable)
    samplefiles = dict()
    
    samples = list(sampledata.getsamples())
    
        

    #print >>sys.stderr, "**READCOUNT**"
    try:
        featurelist = dict()
        trnaloci = dict()
        trnalist = dict()
        ensembllist = dict()
        for currfile in bedfiles:
            featurelist[currfile] = RangeBin(readfeatures(currfile))
        
        for currfile in locifiles:
            trnaloci[currfile] = RangeBin(readbed(currfile))
        for currfile in maturetrnafiles:
            trnalist[currfile] = RangeBin(readbed(currfile))
        if ensemblgtf is not None:    
            embllist = RangeBin(readgtf(ensemblgtf, filtertypes = set()))
        else:
            embllist = None
    
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    
    
    #featurelist = list(readbed(sys.argv[1]))
    #featurelist = list(readbed(sys.argv[1]))
    
    #./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
    '''
    ./countcomplete.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
    '''
    
    featcount = defaultdict(int)
    
    
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
    #print >>sys.stderr, "***"
    
    #lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    counts = defaultdict(lambda: defaultdict(int))
    emblcounts = defaultdict(lambda: defaultdict(int))
    trnacounts = defaultdict(lambda: defaultdict(int))
    trnawholecounts = defaultdict(lambda: defaultdict(int))
    trnafivecounts = defaultdict(lambda: defaultdict(int))
    trnathreecounts = defaultdict(lambda: defaultdict(int))
    trnatrailercounts  = defaultdict(lambda: defaultdict(int))
    trnaantisense  = defaultdict(lambda: defaultdict(int))
    
    trnalocuscounts = defaultdict(lambda: defaultdict(int))
    trnatrailercounts = defaultdict(lambda: defaultdict(int))
    
    fulltrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    partialtrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    trnalocustrailercounts = defaultdict(lambda: defaultdict(int))
    trnaaminocounts = defaultdict(lambda: defaultdict(int))
    
    readlengths = defaultdict(lambda: defaultdict(int))
    trnareadlengths = defaultdict(lambda: defaultdict(int))
    pretrnareadlengths = defaultdict(lambda: defaultdict(int))
    
    aminos = set()
    
    othercounts = defaultdict(int)
    
    maxoffset = 10
    
    bedlist = list(featurelist.iterkeys())
    locilist = list(trnaloci.iterkeys())
    #print >>sys.stderr, bedlist
    #sys.exit()
    '''
    I use readbins rather than searching the bamfiles for regions here because I want to ensure that each read is only mapped to a single feature
    
    '''
    
    '''
    
    samtools sort SRR1508385_nofeat.bam SRR1508385_nofeatsort
    bash ~/pythonsource/tcgatest/gettranscripts.bash SRR1508385_nofeatsort.bam 30 >SRR1508385_nofeat.bed
    
    '''
    trnasamplecounts  = defaultdict(int)
    totalsamplecounts = defaultdict(int)
    emblbiotypes = set()
    
    indelaligns = defaultdict(int)
    
    for currsample in samples:
        currbam = sampledata.getbam(currsample)
        
        try:
            #print >>sys.stderr, currbam
            
            
            if not os.path.isfile(currbam+".bai") or  os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )
            if bamnofeature:
                outname = os.path.splitext(currbam)[0]+"_nofeat.bam"
                outbamnofeature =  pysam.Samfile( outname, "wb", template =  bamfile)
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit()
        
        for i, currread in enumerate(getbamrange(bamfile, primaryonly = True)):
            if i > 10000:
                pass
            isindel = False
    
            gotread = False
            totalsamplecounts[currsample] += 1
            if len(currread.data["CIGAR"]) > 1:
                
                #print >>sys.stderr, currread.bedstring()
                #print >>sys.stderr, currread.data["CIGAR"]
                #deletions
                if sum(curr[1] for curr in currread.data["CIGAR"] if curr[0] == 2):
                    indelaligns[currsample] += 1
                    pass
                #insertions
                if sum(curr[1] for curr in currread.data["CIGAR"] if curr[0] == 2):
                    #indelaligns[currsample] += 1
                    pass
                isindel = True
                continue

            else:
                pass
                #continue
            readlength = len(currread.data['seq'])
            readlengths[currsample][readlength] += 1
            for currbed in locilist:
                for currfeat in trnaloci[currbed].getbin(currread):
                    expandfeat = currfeat.addmargin(30)
                    if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                        pretrnareadlengths[currsample][readlength] += 1
                        trnalocuscounts[currsample][currbed] += 1
                        if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                            fulltrnalocuscounts[currsample][currbed] += 1
                        else:# currread.start + fullpretrnathreshold <  currfeat.start or currread.end - fullpretrnathreshold +3 >  currfeat.end:
                            partialtrnalocuscounts[currsample][currbed] += 1
                            #print >>sys.stderr, "***"
                        gotread = True
                        break
                    if currfeat.getdownstream(30).coverage(currread) > 10:
                        pretrnareadlengths[currsample][readlength] += 1
                        trnalocustrailercounts[currsample][currbed] += 1 
                        #print >>sys.stderr, currfeat.bedstring()
                        gotread = True
                        break
                    elif expandfeat.antisense().coverage(currread) > 5:
                        trnaantisense[currsample][currbed] += 1
                        gotread = True
                        break
            if gotread: 
                continue
            for currbed in trnalist:
                for currfeat in trnalist[currbed].getbin(currread):
                    if currfeat.coverage(currread) > 10:
                        trnareadlengths[currsample][readlength] += 1
                        trnasamplecounts[currsample] += 1
                        trnacounts[currsample][currbed] += 1
                        fragtype = getfragtype(currfeat, currread)
                        trnaaminocounts[currsample][trnainfo.getamino(currfeat.name)] += 1
                        aminos.add(trnainfo.getamino(currfeat.name))
                        if fragtype == "Whole":
                            
                            trnawholecounts[currsample][currbed] += 1
                        elif fragtype == "Fiveprime":
                            trnafivecounts[currsample][currbed] += 1
                        elif fragtype == "Threeprime":
                            trnathreecounts[currsample][currbed] += 1
                        elif fragtype == "Trailer":
                            trnatrailercounts[currsample][currbed] += 1
                        gotread = True
                        break
                            #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                            #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                            #print >>sys.stderr, "****"
                    elif currfeat.antisense().coverage(currread) > 10:
                        trnaantisense[currsample][currbed] += 1
                        gotread = True
                        break
            if gotread: 
                continue
    
            if embllist is not None:
                currtype = None
                for currfeat in embllist.getbin(currread):
                    if currfeat.coverage(currread) > 10: 
                        if currfeat.data["biotype"] == "processed_transcript":
                            #print >>sys.stderr, currfeat.bedstring()
                            
                            pass
    
                        #print >>sys.stderr, "*******"
                        if currtype is None or ncrnaorder[currfeat.data["source"]] > ncrnaorder[currtype]:
                            currtype= currfeat.data["source"]
                            if mitochrom == currread.chrom:
                                currtype = "mito"+currtype
                        
                        
                        
                if currtype is not None:
                    emblcounts[currsample][currtype] += 1
                    emblbiotypes.add(currtype)
                    gotread = True
                        #print >>sys.stderr, currbam +":"+ currbed
            if gotread: 
                continue
            for currbed in bedlist:
                
                for currfeat in featurelist[currbed].getbin(currread):
                    if currfeat.coverage(currread) > 10:
                        counts[currsample][currbed] += 1
                        gotread = True
                        break
                        #print >>sys.stderr, currbam +":"+ currbed
    
            othercounts[currsample] += 1
            if not gotread and embllist is not None and mitochrom == currread.chrom:
                currtype = "Mitochondrial_other"
                emblcounts[currsample][currtype] += 1
                emblbiotypes.add(currtype)
            if not gotread and bamnofeature:
                outbamnofeature.write(currread.data["bamline"])
    #currently not counting indels, but might later.
    for currsample in samples:
        pass
        #print >>sys.stderr, currsample+" indels: "+ str(indelaligns[currsample])+"/"+str(totalsamplecounts[currsample])+":"+str(1.*indelaligns[currsample]/totalsamplecounts[currsample])
    def sumsamples(countdict,sampledata, repname, currfeat = None, sizefactors = defaultdict(lambda: 1)):
        if currfeat is None: #To account for the "other" counts, which don't have a feature
            return sum(countdict[currsample]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
        else:
            sum(sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
            return sum(countdict[currsample][currfeat]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
    #tRNA-Ser-AGA-1-1    
        
    if combinereps:
        replicates = list(sampledata.allreplicates())
        print >>countfile, "\t".join(replicates)
        for currbed in trnalist:
            
            if countfrags:
                #sumsamples(trnafivecounts,sampledata,currrep, currfeat)
                
                print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_other\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor) - (sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor))) for currrep in replicates)
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(sumsamples(trnaantisense,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                
                print  >>countfile, "tRNA\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)        
            
        
        for currbed in locilist:
            #print >>sys.stderr, currbed 
            if countfrags:
                print  >>countfile, "pretRNA_full\t"+"\t".join(str(sumsamples(fulltrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "pretRNA_partial\t"+"\t".join(str(sumsamples(partialtrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "pretRNA_trailer\t"+"\t".join(str(sumsamples(trnalocustrailercounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                print  >>countfile, "pretRNA\t"+"\t".join(str(sumsamples(trnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
        for currbiotype in emblbiotypes:
            print  >>countfile, currbiotype+"\t"+"\t".join(str(sumsamples(emblcounts,sampledata,currrep, currbiotype, sizefactors = sizefactor)) for currrep in replicates)
        for currbed in bedlist:  
            print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(sumsamples(counts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
        print  >>countfile, "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)
    else:
        print  >>countfile, "\t".join(samples)
        for currbed in trnalist:
            
            if countfrags:
                print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_other\t"+"\t".join(str((trnacounts[currsample][currbed] - (trnathreecounts[currsample][currbed] + trnafivecounts[currsample][currbed] + trnawholecounts[currsample][currbed]))/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(trnaantisense[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            else:
                print  >>countfile, currbed+"\t"+"\t".join(str(trnacounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            
            
        
        for currbed in locilist:
            if countfrags:
                print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(fulltrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(partialtrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(trnalocustrailercounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            else:
                print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(trnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbiotype in emblbiotypes:
            print  >>countfile, currbiotype+"\t"+"\t".join(str(emblcounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbed in bedlist:
            print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(counts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        print  >>countfile, "other"+"\t"+"\t".join(str(othercounts[currsample]/sizefactor[currsample]) for currsample in samples)
    
    #print >>sys.stderr, trnaaminocounts
    if trnaaminofilename is not None:
        trnaaminofile = open(trnaaminofilename, "w")
        if combinereps:
            print  >>trnaaminofile, "\t".join(replicates)
            for curramino in aminos:
                #print >>sys.stderr, curramino
                print >>trnaaminofile, curramino+"\t"+"\t".join(str(sumsamples(trnaaminocounts,sampledata,currrep, curramino, sizefactors = sizefactor)) for currrep in replicates)
        else:
            print  >>trnaaminofile, "\t".join(samples)
            for curramino in aminos:
                print >>trnaaminofile, curramino+"\t"+"\t".join(str(trnaaminocounts[currsample][curramino]/sizefactor[currsample]) for currsample in samples)
            
            
    if trnanormfile is not None:
        #samples trnasamplecounts.keys()
        trnanormfile = open(trnanormfile, "w")
        mean = 1.*sum(trnasamplecounts.values())/len(trnasamplecounts.values())
        print >>trnanormfile, "\t".join(samples)
        print >>trnanormfile, "\t".join(str(trnasamplecounts[currsample]/mean) for currsample in samples)
        
    if allreadsnormfile is not None:    
        allreadsnormfile = open(allreadsnormfile, "w")
        mean = 1.*sum(totalsamplecounts.values())/len(totalsamplecounts.values())
        print >>allreadsnormfile,"\t".join(samples)
        print >>allreadsnormfile,"\t".join(str(totalsamplecounts[currsample]/mean) for currsample in samples)
    if readlengthfile is not None:    
        readlengthfile = open(readlengthfile, "w")
        print >>readlengthfile, "Length\tSample\tall\ttrnas\tpretrnas"
        for currsample in readlengths.keys():
            for curr in range(0,max(readlengths[currsample].keys())+1):
                othercount = trnareadlengths[currsample][curr] + pretrnareadlengths[currsample][curr]
                print >>readlengthfile, str(curr)+"\t"+currsample+"\t"+str(readlengths[currsample][curr] - othercount)+"\t"+str(trnareadlengths[currsample][curr]) +"\t"+str(pretrnareadlengths[currsample][curr])
        
        
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--sizefactors',
                       help='Optional file including size factors that will be used for normalization')
    parser.add_argument('--bedfile',  nargs='*', default=list(),
                       help='bed file with non-tRNA features')
    parser.add_argument('--gtffile',  nargs='+', default=list(),
                       help='gtf file with non-tRNA features')
    parser.add_argument('--ensemblgtf',
                       help='ensembl gtf file with tRNA features')
    parser.add_argument('--gtftrnas',  nargs='+', default=list(),
                       help='gtf file with tRNA features')
    parser.add_argument('--trnaloci',  nargs='+', default=list(),
                       help='bed file with tRNA features')
    parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--countfrags', action="store_true", default=False,
                       help='Seperate tRNA fragment types')
    parser.add_argument('--combinereps', action="store_true", default=False,
                       help='Sum samples that are replicates')
    parser.add_argument('--trnatable',
                       help='table of tRNA features')
    parser.add_argument('--trnaaminofile',
                       help='table of tRNAs by amino acid')
    parser.add_argument('--trnanormfile',
                       help='Create normalization file to use to normalize to total tRNA reads')
    parser.add_argument('--allreadsnormfile',
                       help='Create normalization file to use to normalize to total reads')
    parser.add_argument('--readlengthfile',
                       help='optional read lengths table')
    
    
    parser.add_argument('--bamnofeature', action="store_true", default=False,
                       help='Create bam file output for reads without a feature')
    parser.add_argument('--mitochrom',
                       help='Optional name of mitochondrial chromosome in database (Used to specially label mitchondrial features)')
    
    args = parser.parse_args()
    argvars = vars(args)
    main(**argvars)
    #main(samplefile = args.samplefile, sizefactors = args.sizefactors, bedfile = args.bedfile, gtffile=args.gtffile,ensemblgtf=args.ensemblgtf,gtfrnas=args.gtfrnas,trnaloci=args.trnaloci)        
