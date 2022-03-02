library(tidyverse)

# vcf whole genome maf 0.01
sfs <- read.table("Data/sfs/maf01/PWS07_folded.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 

barplot(sfs)

# vcf no maf, ld pruned
#I don't have one.
pop_name <- "PWS07"
sfs <- read.table(paste("Data/moments/sfs/1D/",pop_name,"_folded.sfs",sep = ""), header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 

# angsd downsample 41 whole genome
sfs <- read.table("Data/sfs/downsample_41/unfolded/SS17_minQ20_minMQ30_unfolded.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 
barplot(sfs)


# angsd downsample 41 5MB chr 1
sfs <- read.table("Data/moments/downsample/chr1_5mb/PWS17_folded.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 

sfs <- read.table("Data/moments/downsample/chr1_5mb/TB17_folded.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 

# 2d sfs


# light snp call
sfs <- read.table("Data/sfs/light_snp_call/TB17_folded.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs <- as.vector(colSums(sfs))
barplot(sfs[-c(1,length(sfs))]) #plot variable sites 





model_name <- "sym_mig_2_pop"
pops <- "PWS17_SS17"
pops <- "PWS17_TB17"
pops <- "SS17_TB17"

i <- 2

  #params <- read.delim(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/moments/moments_pipeline/models/",model_name,"/",pops,"/V4_Number_",i,".sym_mig.optimized.txt", sep = ""))
  
  params <- read.delim(paste("Data/moments/downsample/chr1_5mb/",pops,"/V5_Number_",i,".sym_mig.optimized.txt", sep = ""))
  params <- params[order(params$AIC),]
  head(params)
  params <- params %>% separate(Replicate, c("round","replicate"), remove = TRUE,sep = "_R")
  
  boxplot(params$log.likelihood~params$round)
  opt_params <- strsplit(params[1,8], ",")[[1]]
  print("optimized parameter set for empirical data:")
  print(opt_params)
  print("theta")
  print(params$theta[1])
  


# best run

#parameter ranges
params <- read.delim(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/moments/downsample/chr1_5mb/",pops,"/V5_Number_",i,".sym_mig.optimized.txt", sep = ""))
params <- params[order(params$AIC),]
params <- params %>% separate(Replicate, c("round","replicate"), remove = TRUE,sep = "_R")
#print(boxplot(params$log.likelihood~params$round))
params <- params %>% separate(optimized_params.nu1..nu2..m..T., into = c("nu1","nu2","m","t"), sep = ",")
boxplot(params$theta~params$round)
boxplot(as.numeric(params$nu1)~params$round)
boxplot(as.numeric(params$nu2)~params$round)
boxplot(as.numeric(params$m)~params$round)
boxplot(as.numeric(params$t)~params$round)

# calcualte split time and pop size  
params <- read.delim(paste("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/moments/downsample/chr1_5mb/",pops,"/V5_Number_",i,".sym_mig.optimized.txt", sep = ""))
params <- params[order(params$AIC),]
params <- params %>% separate(Replicate, c("round","replicate"), remove = TRUE,sep = "_R")
#print(boxplot(params$log.likelihood~params$round))
opt_params <- strsplit(params[1,8], ",")[[1]]

# Split into two populations, with symmetric migration.
# nu1: Size of population 1 after split.
# nu2: Size of population 2 after split.
# m: Migration rate between populations (2*Na*m)
# T: Time in the past of split (in units of 2*Na generations) 

theta <- as.numeric(params$theta[1])
nu1 <- as.numeric(opt_params[1])
nu2 <- as.numeric(opt_params[2])
m <- as.numeric(opt_params[3])
t <- as.numeric(opt_params[4])
m
t
# generation time
g <- 6
# atlantic herring mutation rate
mu <- 2.0e-9


#https://groups.google.com/g/dadi-user/c/DYrpTHCcC_I/m/nl3f2eGSAQAJ
# estimate of L needs to account for filtering
# L = (genome size) * (filtered SNPs)/(total SNPs)
#L <- 7.26e8  * (200000/8732577)

#using 5MB of chr1
L <- 5000000
#https://groups.google.com/g/dadi-user/c/5HtpVOCO2nc/m/yYbh48MLBAAJ
# effective size of current population
0.003/(4*mu)

# effective size of ancestral population
Nref <- theta / (4* mu* L)
Nref

# Time since split
t*2*Nref*4





