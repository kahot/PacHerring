library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)

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


chr <- "chr12"
chr <- "chr7"

my_vcf <- read.vcfR(paste("Data/vcfs/chrs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_",chr,".vcf.gz", sep = ""))

meta <- read.csv('Data/vcfs/sample_metadata.txt', comment.char = '#', sep = "\t")
my_popmap <- meta[c("Sample","Population_Year")]

# plot variation for individuals grouped by population
new_plot <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = my_popmap,                              
                          cluster        = FALSE,                           
                          snp_label_size = 5000000,                          
                          colour_scheme=c(yel,org,red))

pdf(paste0("Output/inversions/",chr,"_snps.pdf"),width=10,height=8)
combine_genotype_plot(new_plot)
dev.off()


# plot allele frequencies for populations
new_plot <- genotype_plot(vcf_object  =  my_vcf,   
                          popmap = my_popmap,                           
                          snp_label_size = 5000000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=TRUE,                    
                          invariant_filter = TRUE)                      

pdf(paste("Output/inversions/",chr,"_AF.pdf",sep = ""),width=10,height=8)

combine_genotype_plot(new_plot)

dev.off()



# Visualize allele frequency on chromosome 
chr <- "chr12"

my_vcf <- read.vcfR(paste("Data/vcfs/chrs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_",chr,".vcf.gz", sep = ""))

meta <- read.csv('Data/vcfs/sample_metadata.txt', comment.char = '#', sep = "\t")
my_popmap <- meta[c("Sample","Population_Year")]

#my_popmap <- head(my_popmap,100)
pops <- c("PWS17","CA17")

pops <- c("TB17","PWS17","SS17","BC17","WA17","CA17")

my_popmap <-  my_popmap[my_popmap$Population_Year %in% pops,]
my_popmap <-  sample_n(my_popmap[my_popmap$Population_Year %in% pops,],50)


# plot subset of chromosome
new_plot <- genotype_plot(vcf_object  =  my_vcf,   
                          #chr    = 12,                                 
                          #start  = 11700000,                           
                          #end    = 11800000,                           
                          popmap = my_popmap,                           
                          #cluster        = TRUE,                        
                          snp_label_size = 5000000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=TRUE,                    
                          invariant_filter = TRUE)                      

#dendro_with_tips <- new_plot$dendrogram +
#                    geom_text(aes(x=1:length(new_plot$dendro_labels),y=-2.5,label=new_plot$dendro_labels))
#dendro_with_tips

pdf("Output/inversions/test_plot.pdf",width=10,height=8)
pdf("Output/inversions/t2017_plot.pdf",width=10,height=8)

combine_genotype_plot(new_plot)

dev.off()

