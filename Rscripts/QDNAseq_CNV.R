BiocManager::install("QDNAseq")
#https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf

BiocManager::install("BSgenome")
BiocManager::install("GenomicFeatures")




# Generate bin annotations
library(QDNAseq)
library(Biobase)
library(BSgenome)
library(GenomicFeatures)
library (AnnotationDbi)
# Convert the gff3 file to usable format

C.harengus <- GenomicFeatures::makeTxDbFromGFF("Data/annotations/Clupea_harengus.Ch_v2.0.2.100.gff3",format="auto")

# Save as a database
saveDb(C.harengus,"Data/annotations/C_harengus")

#C.harengus<-laodDb("Data/annotations/C_harengus")

columns(C.harengus)
#intergenic <- gaps(C.harengus) 
#head (intergenic)

#Change the name of each chromosome files
sink("C_harengus/rename_files.sh")
for (i in 1:26){
    cat(paste0("mv -n Data/annotations/chr/",i," C_harengus/C_harengus_chr",i,".fa \n"))
}
sink(NULL)



#forging the Charengus pacakge
#myfile<-read.dcf("seed_file.txt", fields = NULL, all = FALSE, keep.white = NULL)
#write.dcf(myfile , file = "C_harengus/seed.dcf", append = FALSE, useBytes = FALSE, indent = 0.1 * getOption("width"), width = 0.9 * getOption("width"), keep.white = NULL)

unlink("C_harengusTxdb", recursive=TRUE, force=TRUE)
forgeBSgenomeDataPkg("C_harengus/seed.dcf", seqs_srcdir="C_harengus/",destdir="C_harengus/")

# Install the BSgenome.Charengus.NCBI.2.0.2_2.0.2.tar.gz pacakge

#@ Following the QDNAseq manual ##
library('BSgenome.Charengus.NCBI.2.0.2')
# Set the bin size
binsize<-15
## create bins from the reference genome
bins<-createBins(bsgenome=BSgenome.Charengus.NCBI.2.0.2, binSize=binsize)

# calculate mappabilites per bin from ENCODE mapability tracks
#1.Install genmap using conda to create wig file 
#2. Install wigToBigWig and convert to .bw file 
#3.Install bigWigAverageOverBed and create a file

#create congif file for BICseq2-norm
config<-data.frame(chr=paste0("chr",1:26))
config$fa<-paste0("chr",1:26,".fa")
config$map<-
    

bins$mappability <- calculateMappability(bins,
                        bigWigFile="Data/annotations/Charegus.bw",
                        bigWigAverageOverBed="Data/annotations/Charengus.tab")
