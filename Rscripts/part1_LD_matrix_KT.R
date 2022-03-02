# From https://github.com/EleniLPetrou/pacific_herring_RADseq
# conduct network analysis of LD

################################################################################
library(vcfR)
library(genetics)
library(ggplot2)
library(tidyr)
library(dplyr)

################################################################################
#To run this code, put all of your vcf files in a single directory

#setwd

#setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/LD_GCA900700415")
#list.files()

# set file names

input_fileName <- "Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf"
output_fileName <- "Output/LD/ph_filtered_snps_LDmat.txt"

#vcf is too big. divide them up into chromosomes
library(stringr)

sink("Data/vcfs/subset_by_chromosome.sh")
cat("#!/bin/bash\n")
#26 chromosomes
for (i in 1: 26){
    ch<-paste0("chr",i)
    phrase<-paste0('bcftools view ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz -r ',ch,' -Oz --threads 10 > ph_',ch,".vcf.gz")
    cat(phrase)
    cat("\n")
}
sink(NULL)

#Run the bash file


################################################################################
#read in vcf files in directory using vcfR, and start data processing
vfiles<-list.files("Data/vcfs", pattern="ph_chr")


for (i in 1: length(vfiles)){
    vcf_data <- read.vcfR(paste0("Data/vcfs/", vfiles[i]))
  
    #save metadata as dataframe - you will need this later for plotting
    vcf_df <- as.data.frame(vcf_data@fix)
    write.csv(vcf_df, paste0("Output/LD/gt_metadata_chr",i,".csv"))
    #use vcfR to extract the genotypes from the vcf file --> make matrix
    gt <- extract.gt(vcf_data, return.alleles = TRUE)
    #gt[1:10, 1:10] #take a peek
  
    #Prepare data for genetics package and LD calculation
    #transpose the genotype matrix and make it into a data frame.
    #some are read in as A|A rather than A/A. Change | to /
    gt<-apply(gt,1,function(x) gsub("\\|","\\/",x))
    gt_df <- data.frame(t(gt))
    
    # change "." to NAs
    gt_df[gt_df == "."] <- NA
    write.csv(gt_df, paste0("Output/LD/gt_df_chr",i,".csv"))
    
    # Now that you have a dataframe of genotypes, turn it into a genetics object
    df_genetics <- makeGenotypes(df_genetics, sep="/")
    write.csv(gt_genetics, paste0("Output/LD/gt_genetics_chr",i,".csv"))
    
    
}

out<-LD(df_g2)

for (i in 1: length(vfiles)){
    gt_genetics<-read.csv(paste0("Output/LD/gt_genetics_chr",i,".csv"), stringsAsFactors = T, row.names = 1)
    
    #How do you eliminate the loci with 1 level only? 
    #gt_genetics<-gt_genetics[, sapply(gt_genetics, nlevels) > 1]
    
    out<-LD(gt_genetics)
    
    # clear up some space in memory
  
    remove(vcf_data)
    remove(vcf_df)
    remove(gt)
    remove(gt_df)

    # Run the LD test using the R package genetic
    # The LD function computes pairwise linkage disequilibrium between genetic markers
    output<- LD(gt_genetics)
     
    #it ran over a day and did not fishin with the following message:
    #Warning message:
    #    In LD.data.frame(gt_genetics) :
    #    Non-genotype variables or genotype variables with more or less than two alleles detected. These variables will be omitted: chr1_24456, chr1_24457, chr1_149865, chr1_150325, chr1_152032, chr1_152931, 
    #    chr1_178489, chr1_178538, chr1_178541, chr1_178917, chr1_178918, chr1_234578, chr1_282544, chr1_284981, chr1_294561, chr1_297295, chr1_319818, chr1_330660, chr1_330672, chr1_330685, chr1_362770, chr1_488615, chr1_529874, chr1_529894, chr1_530313, chr1_570388, chr1_611148, chr1_611152, chr1_724476, chr1_737866, chr1_756871, chr1_758618, chr1_762501, chr1_768017, chr1_895437, chr1_981101, chr1_996707, chr1_996715, chr1_1010462, chr1_1010486, chr1_1010487, chr1_1010500, chr1_1010550, chr1_1011331, chr1_1011334, chr1_1011405, chr1_1011422, chr1_1108997, chr1_1112110, chr1_1304933, chr1_1318949, chr1_1318950, chr1_1332008, chr1_1506699, chr1_1542003, chr1_1542008, chr1_1546116, chr1_1636236, chr1_1636708, chr1_1637847, chr1_1639040, chr1_1639041, chr1_1656460, chr1_1677223, chr1_1687452, chr1_16 [... truncated]
    
     
    # write out the R2 matrix to a file for LDNA
    r2_mat <- output$`R^2`
    write.table(r2_mat, file=paste0("Output/LD/LDmatrix_chr",i,".txt"), row.names= TRUE, col.names= TRUE , sep = "\t")
 }

  