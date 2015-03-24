
#Download and remove adapters from small RNA sequencing studies
~/.aspera/connect/bin/ascp -QT -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR029/SRR029131/SRR029131.sra ./
fastq-dump -Z SRR029131.sra | ~/pythonsource/seqqa/cutadapt-1.2.1/bin/cutadapt -m 15 --adapter='TCGTATGCCGTCTTCT' - |  gzip -c  >SRR029131.fastq.gz

~/.aspera/connect/bin/ascp -QT -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR029/SRR029124/SRR029124.sra ./
fastq-dump -Z SRR029124.sra | ~/pythonsource/seqqa/cutadapt-1.2.1/bin/cutadapt -m 15 --adapter='TCGTATGCCGTCTTCT' - |  gzip -c  >SRR029124.fastq.gz

~/.aspera/connect/bin/ascp -QT -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR207/SRR207111/SRR207111.sra ./
fastq-dump -Z SRR207111.sra | ~/pythonsource/seqqa/cutadapt-1.2.1/bin/cutadapt -m 15 --adapter='CGTATGCCGTCT' - |  gzip -c  >SRR207111.fastq.gz 

~/.aspera/connect/bin/ascp -QT -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR207/SRR207116/SRR207116.sra ./
fastq-dump -Z SRR207116.sra | ~/pythonsource/seqqa/cutadapt-1.2.1/bin/cutadapt -m 15 --adapter='CGTATGCCGTCT' - |  gzip -c  >SRR207116.fastq.gz

#Download and combine hg19 chromosomes 
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvf chromFa.tar.gz -O > hg19.fa

#Download ensembl GTF, change chromosome names

wget -q -O - ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz | gzip -cd | grep -v '^#' | awk '{print "chr" $0;}' | grep -e Mt_rRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA  >hg19-genes.gtf


#scan for nuclear tRNAs
tRNAscan-SE hg19.fa | grep -v chrM >hg19-tRNAs.txt

#Create the tRNA database
maketrnadb.bash hg19 hg19-tRNAs.txt hg19.fa


#Map the tRNAreads
mapreads.bash TestTrnas hg19 TrnaSamples.txt hg19-genes.gtf