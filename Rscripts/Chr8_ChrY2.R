#Use maf00 files to see if we can detect more distinct differences between sexes

library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)
library(ape)
library(reshape2)
# color-blind friendly 
# Wong, B. Points of view: Color blindness. Nat Methods (2011).
bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'


# Visualize variation on chromosome 

my_vcf <- read.vcfR("Output/chr8/maf00/PWS_maf00_chr8sex.vcf.gz")

meta <- read.csv("Data/Sample_metadata_892pops.csv")
my_popmap <- meta[c("Sample","Population.Year")]

pws<-
ad <- extract.gt(my_vcf, element = "AD")
AD_frequency(ad, allele = 1L, sum_type = 0L, decreasing = 1L)

#read the variant position info for chrY of Atlantic Herring
chrY_vars<-read.csv("Data/chrY/pnas.2009925117.sd01.csv")
#overlapping positions from Atlantic Herring
poss<-getPOS(my_vcf)
sex.snps<-intersect(chrY_vars$Pos, poss)

sex.snps[1]
sex.snps[length(sex.snps)]

### group individuals based on GT (homo/hetero)
#Start with PWS
pws<-my_popmap[grep("PWS", my_popmap$Sample),]
pws1 <- genotype_plot(vcf_object  =  my_vcf,
                       popmap = pws,  
                      start = sex.snps[1],
                      end = sex.snps[length(sex.snps)],
                       cluster        = T,                           
                       snp_label_size = 50000, colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6') )
combine_genotype_plot(pws1)

#extract the genotype data
pws_gt<-pws1$genotypes[["data"]]

#Add the sample id
pws_gt$Sample<-pws$Sample

#first pull out all homozygous individuals

pwsGT<-data.frame(Sample=pws$Sample)
for (i in 1:length(pws$Sample)){
    id<-pws$Sample[i]
    df<-pws_gt[pws_gt$Sample==id,]
    #remove missing locus
    df<-df[!is.na(df$GT),]
    #count homozygous sites
    pwsGT$homo.loci[i]<-nrow(df[df$GT==0,])
    pwsGT$hetero.loci[i]<-nrow(df[df$GT==2,])
    pwsGT$homo_percetn[i]<-pwsGT$homo.loci[i]/nrow(df)*100
    pwsGT$hetero_percetn[i]<-pwsGT$hetero.loci[i]/nrow(df)*100
}


#dilute the results and don't see any pattern
