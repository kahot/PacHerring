#LD Decay Analysis
library(ggplot2)
library(tidyverse)
library(reticulate)
library(reshape2)
library(plyranges)
library(seqinr)
library(data.table)


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



## LD Decay 

#No big differences in LD decay between populations. Drops to background levels between 200bp 
#(Atlantic Herring LD drops off around 100bp).

pop_names = c("PWS91","PWS96","PWS07","PWS17")
pop_line <- c(4,3,2,1)
pop <- "PWS91"
i <- 1
ld <- read.table(paste0("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld"), header = TRUE, stringsAsFactors = FALSE)
ld$distance <- ld$BP_B-ld$BP_A

plot(ld$distance, ld$R2, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

for (i in c(2:length(pop_names))){
    ld <- read.table(paste("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    ld$distance <- ld$BP_B-ld$BP_A
    
    #ld <- ld[ld$CHR_A == "1" & ld$CHR_B == "1" ,]
    
    smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
}


#Create more zoomed in view to assess the cutoff
pdf("Output/LD/LDdecay_PWS_zoomedin.pdf", height = 4, width = 6)
plot(ld$distance, ld$R2, type="n", xlim = c(0,150), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

for (i in c(2:length(pop_names))){
    df <- read.table(paste("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    df$distance <- df$BP_B-df$BP_A
      
    smooth_ld <-smooth.spline(df$distance, df$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
}
dev.off()


pdf("Output/LD/LDdecay_PWS91_zoomedin.pdf", height = 4, width = 6)
plot(ld$distance, ld$R2,col="gray60", pch=".", xlim = c(0,200), xlab = "distance (bp)" ,ylab = "R2")
dev.off()
pdf("Output/LD/LDdecay_PWS17_zoomedin.pdf", height = 4, width = 6)
plot(df$distance, df$R2,col="gray60", pch=".", xlim = c(0,200), xlab = "distance (bp)" ,ylab = "R2")
dev.off()

### 
Plot SFS
#for (pop_name in pop_names){
sfs1 <- read.table("Output/fst_pbs/folded_PWS07_TB91.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs1 <- data.frame(t(sfs1))
sfs1<-as.vector(sfs1[,1])
barplot(sfs1[-1])

sfs2 <- read.table("Output/fst_pbs/folded_PWS07_TB91_2.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs2 <- data.frame(t(sfs2))
sfs2<-as.vector(sfs2[,1])
barplot(sfs2[-1])


### SUbset the PL files
samples<-read.table("Data/ngsLD/samplesnames.txt")
pw07<-grep("PWS07",samples[1,]) #406-543
pw17<-grep("PWS17",samples[1,]) #544-711
pw91<-grep("PWS91",samples[1,])  # 712-885
pw96<-grep("PWS96",samples[1,]) #886-1101

## how to prep the file for ngsLD ###
#1. extract the columns corresponding to the population.year
cut -d "," -f406-543 /home/ktist/ph/data/new_vcf/MD7000/beagle/MD7000_maf05_c12.BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/c12.PWS07
#2. remove the header (sample names)
sed -e '1d' < /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/c12.PWS07 > /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/c12.PWS07.noheader  
#3. zip the file 
bgzip /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/c12.PWS07.noheader  
#4. run ngsLD
./programs/ngsLD/ngsLD --geno /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/c12.PWS07.noheader.gz --max_kb_dist 1000 --n_ind 46 --n_sites 11759 --prob --pos /home/ktist/ph/data/new_vcf/MD7000/beagle/positions/pos_chr12.txt --n_threads 12 --out /home/ktist/ph/data/ngsLD/chr12.PWS07.ld 


# Read the output file from ngsLD

ldpws<-read.table("Data/ngsLD/chr12.PWS07.ld", header = F, sep="\t" )


plot(ldpws$V3, ldpws$V4, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ldpws$V3, ldpws$V4, spar = .15)
lines(smooth_ld, lwd = 2,lty = 2, col = red)

plot(ldpws$V3, ldpws$V4, col = "gray", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2", pch=".")

ggplot()+
    geom_point(data=ldpws[ldpws$V3<=1000,], aes(x=V3, y=V4),color = "gray", size=0.2)+
    
    geom_smooth(data=ldpws[ldpws$V3<=1000,], aes(x=V3, y=V4), method = "loess", se = FALSE, span = 1, color = blu)+
     theme_minimal()+
    ylim(0,1)+xlab("Distance (bp)")+ylab('LD')
ggsave("Output/LD/ngsLD_decay_curve_ch12_PWS07.pdf", width = 7, height = 5)

plot(ld1$R2, pch=".", col="gray")
ld_sm<-ld1[ld1$distance<=1000,]

plot(ld_sm$R2, pch=".", col="gray")


#heatmap
library(snpStats)
library(LDheatmap)

sample <- read.pedfile("Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.ped.gz", snps="Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.info.gz")

gdat<-sample$genotypes
dist<-read.table("Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.info")
dist<-dist[,-1]
LDheatmap(gdat, genetic.distances=dist, distances="physical",
          LDmeasure="r", title="Pairwise LD CA17", add.map=TRUE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=NULL, color=NULL, newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=NULL, text=FALSE)

sample <- read.pedfile("Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.ped", snps="Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.info")

gdat<-sample$genotypes
dist<-read.table("Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.info")
dist<-dist[,-1]

LDheatmap(gdat, genetic.distances=dist, distances="physical",
          LDmeasure="r", title="Pairwise LD PWS07 chr12", add.map=TRUE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=NULL, color=NULL, newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=NULL, text=FALSE)
