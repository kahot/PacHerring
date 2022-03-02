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
