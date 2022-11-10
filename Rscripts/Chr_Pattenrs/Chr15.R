#### Plot which individuals have freq changes in this region
library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)
bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'


vcf <- read.vcfR("Data/new_vcf/PH_MD7000_maf05_chr7.vcf.gz")
vcf <- read.vcfR("Data/new_vcf/PH_MD7000_maf05_chr15.vcf.gz")

vcf <- read.vcfR("Data/new_vcf/PWS/PWSonly_chr15/PWSonly_ch15_maf05.vcf")


pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[c("Sample","Population.Year")]
pops<-unique(pop_info$Population.Year)
pops<-pops[grep("PWS",pops)]

heights=c(1,9)
widths=c(2,3)

#for (i in 2:length(pops)){
for (i in 3:4){
    pop<-pops[i]
    
    p_map<-pop_info[pop_info$Population.Year==pop,]
    
    chr_plot <- genotype_plot(vcf_object  =  vcf,   
                              popmap = p_map,                           
                              snp_label_size = 1000000,                     
                              colour_scheme=c(yel,"red","blue"),
                              plot_allele_frequency=F,                    
                              invariant_filter = F)  
    
    #Change the na color to lighter gray
    gt_plot<-chr_plot$genotypes   
    updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                       breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
    chr_plot$genotypes<-updated    
    
    
    #pdf(paste0("Output/chr7/DP7000/", pop, "_chr7.pdf"),width=12,height=10)
     
    #plot_grid(chr_plot$positions,
    #          chr_plot$genotypes,
    #          ncol=1,nrow=2,axis="tblr",align="v",rel_heights = heights)
    
    pdf(paste0("Output/chr/DP7000/chr15/PWSonly_", pop, "_chr15.pdf"),width=12,height=10)
    combine_genotype_plot(chr_plot)
    dev.off()
    
    #gt<-chr_plot$genotypes[["data"]]
    #gt$Sample<-p_map$Sample
    #write.csv(gt, paste0("Output/chr7/DP7000/", pops[i], "_gtinfo.csv"))
 #
}

    
    
    chr7_plot <- genotype_plot(vcf_object  =  vcf,   
                           popmap = pop_info,                           
                           snp_label_size = 1000000,                     
                           colour_scheme=c(yel,"red","blue"),
                           plot_allele_frequency=F,                    
                           invariant_filter = F)                      

pdf(paste0("Output/AF/Chr7/chr7.pdf"),width=10,height=8)
combine_genotype_plot(chr7_plot)
dev.off()

pw07_map<-pws_popmap[pws_popmap$Population.Year=="PWS07",]

pw07_plot <- genotype_plot(vcf_object  =  pws_vcf,   
                           popmap = pw07_map,    
                           start=17000000,
                           end=29000000,
                           snp_label_size = 1000000,                     
                           colour_scheme=c(yel,org,red),
                           plot_allele_frequency=F,                    
                           invariant_filter = TRUE)                      

pdf(paste0("Output/AF/Chr7/chr7_SNPs_PWS07_focused.pdf"),width=10,height=8)
combine_genotype_plot(pw07_plot)
dev.off()

pw07_plot2 <- genotype_plot(vcf_object  =  pws_vcf,   
                            popmap = pw07_map,    
                            start=17000000,
                            end=29000000,
                            snp_label_size = 1000000,                     
                            colour_scheme=c(yel,org,red),
                            plot_allele_frequency=F,                    
                            invariant_filter = TRUE,
                            cluster=T
)      

pdf(paste0("Output/AF/Chr7/chr7_SNPs_PWS07_focused2.pdf"),width=10,height=8)
combine_genotype_plot(pw07_plot2)
dev.off()

#extract the genotype data
pws_gt<-pw07_plot$genotypes[["data"]]

#Add the sample id
pws_gt$Sample<-pw07_map$Sample

#Find which sample has inversion between 18-28M

pws_alt<-pws_gt[pws_gt$GT==2,]
pws_alt<-pws_alt[!is.na(pws_alt$snp),]

table(pws_alt$Sample)

#narrow the loci
pws_alt<-pws_alt[pws_alt$snp>=18000000&pws_alt$snp<=28000000,]
table(pws_alt$Sample)
#0435_PWS07 909

inv<-data.frame(table(pws_alt$Sample))
inv<-inv[order(inv$Freq, decreasing = T),]
#        Var1 Freq
#25 0435_PWS07  909
#5  0357_PWS07  258
#19 0412_PWS07  166

## Look at all populations 
vcf <- read.vcfR("Output/AF/Chr7/chr7_bigchang_loci2.vcf")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
popmap<-pop_info[c("Sample","Population.Year")]
pops<-unique(popmap$Population.Year)

for (i in 3:length(pops)){
    for (i in 8:8){
        print(pops[i])
        p_map<-popmap[popmap$Population.Year==pops[i],]
        chr7_plot <- genotype_plot(vcf_object  =  vcf,   
                                   popmap = p_map,                           
                                   snp_label_size = 1000000,                     
                                   colour_scheme=c(yel,org,"red"),
                                   plot_allele_frequency=F,                    
                                   invariant_filter = TRUE)  
        
        pdf(paste0("Output/AF/Chr7/", pops[i], "_chr7_18-28M.pdf"),width=10,height=8)
        combine_genotype_plot(chr7_plot)
        dev.off()
        
        gt<-chr7_plot$genotypes[["data"]]
        gt$Sample<-p_map$Sample
        write.csv(gt, paste0("Output/AF/Chr7/", pops[i], "_chr7_18-28M_gtinfo.csv"))
    }
    
    
    pw07_plot2 <- genotype_plot(vcf_object  =  pws_vcf,   
                                popmap = pw07_map,    
                                start=17000000,
                                end=29000000,
                                snp_label_size = 1000000,                     
                                colour_scheme=c(yel,org,red),
                                plot_allele_frequency=F,                    
                                invariant_filter = TRUE,
                                cluster=T
    )      
    
    
    ad <- extract.gt(pws_vcf, element = "AD")
    
    AD_frequency(ad, allele = 1L, sum_type = 0L, decreasing = 1L)
    
    
    
    