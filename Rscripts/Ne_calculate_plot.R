
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


# Summary
#This notebook estimates long-term effective population size (Ө = 4Neµ) and plots Ne for each population. Plots short-term estimates of Ne estimated by [Vince Buffalo](https://vincebuffalo.com/). Measures linkage disequilibrium with PLINK and plots LD decay.

# Ne estimates from Watterson's theta

#Uses sfs generated in the [popgen_stats notebook](https://github.com/joemcgirr/pac_herring/blob/master/Rmarkdown/popgen_stats/popgen_stats.Rmd)

# http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/angsd-activity-sfs-fst-pbs/

pop_names = c("PWS91","PWS96","PWS07","PWS17","TB91","TB96","TB06","TB17","SS96","SS06","SS17","BC17","WA17","CA17")
pop <- "PWS07"
ne_s <- c()
for (pop in pop_names){
    x<-scan(paste("Data/sfs/downsample_41/folded/",pop,"_minQ20_minMQ30_folded.sfs", sep = ""))
    
    nSites<-sum(x)   #Number of sites where we have data
    nSeg<-sum(x[c(-1)])    #Number of segregating sites
    an <- function(n) sum(1/1:(n-1))
    thetaW <- nSeg/an(length(x[c(-1)])/2) # Wattersons Theta
    #print(pop)
    #print("effective population size")
    ne <- thetaW / 2.0e-9 / nSites / 4 # effective population size
    #print(ne)
    ne_s <- c(ne_s, ne)
}

#png("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/figs/ne.png", height = 5, width = 8, units = 'in', res = 600)

barplot(ne_s, names.arg = pop_names, ylab = expression('N '[e]), col = c(red,red,red,red,blu,blu,blu,blu,yel,yel,yel,grb,lir,org),cex.names=0.6)

#dev.off()

ne <- data.frame(pop = pop_names, Ne = ne_s)
ne <- ne %>% separate(pop, 
                      into = c("pop", "year"), 
                      sep = "(?<=[A-Za-z])(?=[0-9])")
ne$year <- sub("07", "06", ne$year)
ne$year <- sub("91", "1991", ne$year)
ne$year <- sub("96", "1996", ne$year)
ne$year <- sub("06", "2006", ne$year)
ne$year <- sub("17", "2017", ne$year)
ne <- ne[1:11,]

#png("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/figs/ne_watterson.png", height = 5, width = 8, units = 'in', res = 600)

ggplot(ne, aes(year, Ne, group = pop)) +
    geom_point(aes(color = pop)) +
    geom_line(aes(color = pop),size = 2) +
    scale_color_manual(values = c(red,yel,blu))+
    #geom_errorbar(aes(ymin = lower, ymax = upper), width=.05) +
    labs(x = "", y = "Ne\n") +
    theme_classic()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=14))
#dev.off()

```

# Ne estimates from temporal covariance

Vince measured the variance in temporal changes in allele frequency across sampling periods, and this measure is inversely proportional to effective population size.

```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

# uses vince's results from ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf

ne <- read.delim("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/vince/herring_Ne_estimates.tsv")
ne$time <- paste(ne$start, ne$end, sep = " ")
ne$time <- sub("2007", "2006", ne$time)


#png("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/figs/ne_temporal.png", height = 5, width = 8, units = 'in', res = 600)

ggplot(ne, aes(time, Ne, group = sample)) +
    geom_point(aes(color = sample)) +
    geom_line(aes(color = sample),size = 2) +
    scale_color_manual(values = c(red,yel,blu))+
    geom_errorbar(aes(ymin = lower, ymax = upper), width=.05) +
    labs(x = "", y = "Ne\n") +
    theme_classic()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=14))
#dev.off()

```

# Linkage Disequilibrium

Estimate LD with PLINK

example script

```{bash,eval=FALSE,class.source = 'fold-show'}

#!/bin/bash

#SBATCH --job-name=BC17_LD
#SBATCH --mem=8G 
#SBATCH --ntasks=4 
#SBATCH -e BC17_LD_%A_%a.err 
#SBATCH --time=1:00:00 
#SBATCH --mail-user=jamcgirr@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high 

#module load vcftools 
#vcftools --gzvcf /home/jamcgirr/ph/data/vcfs/split_pops/maf05/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --plink --out /home/jamcgirr/ph/data/plink/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05 

module load plink 
#plink --file /home/jamcgirr/ph/data/plink/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05 --indep-pairwise 100 10 0.8 --r2 --out /home/jamcgirr/ph/data/plink/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8 --threads 4 
plink --file /home/jamcgirr/ph/data/plink/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05 --r2 --out /home/jamcgirr/ph/data/plink/population_BC17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_r2 --threads 4 

#command to run: sbatch script_BC17_LD.sh

```

## Generate LD scripts

```{python,results='hide', eval = FALSE}
# LD with plink

job_name = 'LD'
vcf_name = 'ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05'

sbatch_header(job_name,'8','8','24')
script = 'script_' + job_name + '.sh'
o = io.open(script,'a+', newline='\n')


o.write('module load plink \n')
o.write('plink --file /home/jamcgirr/ph/data/vcfs/'+vcf_name+' --indep-pairwise 100 10 0.8 --r2 --out /home/jamcgirr/ph/data/plink/'+vcf_name+'_indep_pairwise_100_10_0.8 --threads 8 \n') 


#run sbatch submission 
o.write('\n\n#command to run: sbatch '+script)
o.close()

# ~2 min

# LD for each population
job_name = '_LD'
infiles = ["BC17","CA17","PWS07","PWS17","PWS91","PWS96","SS06","SS17","SS96","TB06","TB17","TB91","TB96","WA17"]
vcf_name = 'ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05'
for infile in infiles:
    script = 'script_' + infile + job_name + '.sh'
sbatch_header_loop(job_name,'8','4','1', infile)
o = io.open(script,'a+', newline='\n')


o.write('#module load vcftools \n')
o.write('#vcftools --gzvcf /home/jamcgirr/ph/data/vcfs/split_pops/maf05/population_'+infile+'_'+vcf_name+'.vcf.gz --plink --out /home/jamcgirr/ph/data/plink/population_'+infile+'_'+vcf_name+' \n\n')

o.write('module load plink \n')
o.write('#plink --file /home/jamcgirr/ph/data/plink/population_'+infile+'_'+vcf_name+' --indep-pairwise 100 10 0.8 --r2 --out /home/jamcgirr/ph/data/plink/population_'+infile+'_'+vcf_name+'_indep_pairwise_100_10_0.8 --threads 4 \n') 
o.write('plink --file /home/jamcgirr/ph/data/plink/population_'+infile+'_'+vcf_name+' --r2 --out /home/jamcgirr/ph/data/plink/population_'+infile+'_'+vcf_name+'_r2 --threads 4 \n') 

o.write('\n\n#command to run: sbatch '+script)
o.close()
# ~ 5 min
```

## LD Decay 

No big differences in LD decay between populations. Drops to background levels between 200bp (Atlantic Herring LD drops off around 100bp).

```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

pop_names = c("PWS91","PWS96","PWS07","PWS17")
pop_line <- c(4,3,2,1)
pop <- "PWS91"
i <- 1
ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
ld$distance <- ld$BP_B-ld$BP_A

plot(ld$distance, ld$R2, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

for (i in c(2:length(pop_names))){
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    ld$distance <- ld$BP_B-ld$BP_A
    
    #ld <- ld[ld$CHR_A == "1" & ld$CHR_B == "1" ,]
    
    smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
}


pop_names = c("TB91","TB96","TB06","TB17")
pop_line <- c(4,3,2,1)
pop <- "TB91"

i <- 1
ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
ld$distance <- ld$BP_B-ld$BP_A

plot(ld$distance, ld$R2, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "r2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = blu)

for (i in c(2:length(pop_names))){
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    ld$distance <- ld$BP_B-ld$BP_A
    
    #ld <- ld[ld$CHR_A == "1" & ld$CHR_B == "1" ,]
    
    smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = blu)
}


```


```{r, echo=FALSE, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7, eval = FALSE}

pop_names = c("PWS91","PWS96","PWS07","PWS17")
pop_line <- c(4,3,2,1)
pop <- "PWS91"

ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/angsd/ngsLD/chr1_thin20k_",pop,"_60k_pairs_ld.txt",sep = ''), header = FALSE, stringsAsFactors = FALSE)
colnames(ld) <- c("snpA","snpB","distance","r_squared_pearson","D_EM","D_prime_EM","r_squared_EM")
ld[mapply(is.infinite, ld)] <- NA
ld <- na.omit(ld)

#png("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/figs/population_structure/LD_PWS.png", height = 6, width = 7, units = 'in', res = 600)


plot(ld$distance, ld$r_squared_EM,xlim = c(0,1000), col = "white",ylim = c(0.1,0.6), xlab = "distance (bp)", ylab = expression(paste("r"^" 2")), main = "PWS")
legend(835, 0.6, legend=c("91", "96","06/07","17"),col="black", lty=c(4,3,2,1), cex=1,bty = "n" )

for (i in c(1:length(pop_names))){
    
    # plot LD decay with smooth spline
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/angsd/ngsLD/chr1_thin20k_",pop_names[i],"_60k_pairs_ld.txt",sep = ''), header = FALSE, stringsAsFactors = FALSE)
    colnames(ld) <- c("snpA","snpB","distance","r_squared_pearson","D_EM","D_prime_EM","r_squared_EM")
    ld[mapply(is.infinite, ld)] <- NA
    ld <- na.omit(ld)
    
    smooth_ld <-   smooth.spline(ld$distance, ld$r_squared_EM, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
    
}

#dev.off()

# TB

pop_names = c("TB91","TB96","TB06","TB17")
pop_line <- c(4,3,2,1)
pop <- "PWS91"

plot(ld$distance, ld$r_squared_EM,xlim = c(0,1000), col = "white",ylim = c(0.1,0.6), xlab = "distance (bp)", ylab = expression(paste("r"^" 2")), main = "TB")
legend(835, 0.6, legend=c("91", "96","06/07","17"),col="black", lty=c(4,3,2,1), cex=1,bty = "n" )


for (i in c(1:length(pop_names))){
    
    # plot LD decay with smooth spline
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/angsd/ngsLD/chr1_thin20k_",pop_names[i],"_60k_pairs_ld.txt",sep = ''), header = FALSE, stringsAsFactors = FALSE)
    colnames(ld) <- c("snpA","snpB","distance","r_squared_pearson","D_EM","D_prime_EM","r_squared_EM")
    ld[mapply(is.infinite, ld)] <- NA
    ld <- na.omit(ld)
    
    smooth_ld <-   smooth.spline(ld$distance, ld$r_squared_EM, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = blu)
    
}

# SS
pop_names = c("SS96","SS06","SS17")
pop_line <- c(3,2,1)
pop <- "PWS91"

#png("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/figs/population_structure/LD_SS.png", height = 6, width = 7, units = 'in', res = 600)

plot(ld$distance, ld$r_squared_EM,xlim = c(0,1000), col = "white",ylim = c(0.1,0.6), xlab = "distance (bp)", ylab = expression(paste("r"^" 2")), main = "SS")
legend(835, 0.6, legend=c("96","06/07","17"),col="black", lty=c(3,2,1), cex=1,bty = "n" )


for (i in c(1:length(pop_names))){
    
    # plot LD decay with smooth spline
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/angsd/ngsLD/chr1_thin20k_",pop_names[i],"_60k_pairs_ld.txt",sep = ''), header = FALSE, stringsAsFactors = FALSE)
    colnames(ld) <- c("snpA","snpB","distance","r_squared_pearson","D_EM","D_prime_EM","r_squared_EM")
    ld[mapply(is.infinite, ld)] <- NA
    ld <- na.omit(ld)
    
    smooth_ld <-   smooth.spline(ld$distance, ld$r_squared_EM, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = yel)
    
}

#dev.off()

# southern 2017s
pop_names = c("BC17","WA17","CA17")
pop_cols = c(grb,lir,org)
pop_line <- c(3,2,1)
pop <- "PWS91"

plot(ld$distance, ld$r_squared_EM,xlim = c(0,1000), col = "white",ylim = c(0.1,0.6), xlab = "distance (bp)", ylab = expression(paste("r"^" 2")), main = "BC WA CA")
legend(835, 0.6, legend=c("BC17","WA17","CA17"),col=c(grb,lir,org), lty=c(1,1,1), cex=1,bty = "n" )


for (i in c(1:length(pop_names))){
    
    # plot LD decay with smooth spline
    ld <- read.table(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/angsd/ngsLD/chr1_thin20k_",pop_names[i],"_60k_pairs_ld.txt",sep = ''), header = FALSE, stringsAsFactors = FALSE)
    colnames(ld) <- c("snpA","snpB","distance","r_squared_pearson","D_EM","D_prime_EM","r_squared_EM")
    ld[mapply(is.infinite, ld)] <- NA
    ld <- na.omit(ld)
    
    smooth_ld <-   smooth.spline(ld$distance, ld$r_squared_EM, spar = .15)
    lines(smooth_ld, lwd = 2,lty = 1, col = pop_cols[i])
    
}

## Stairway plot 
#https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12985
#https://github.com/popgenmethods/smcpp#quick-start-guide

```






