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

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[c("Sample","Population.Year")]
pops<-unique(pop_info$Population.Year)

Plots<-list()
for (i in 2:length(pops)){
    pop<-pops[i]
    
    p_map<-pop_info[pop_info$Population.Year==pop,]
    
    chr_plot <- genotype_plot(vcf_object  =  vcf,   
                              popmap = p_map,                           
                              #start=14000000,
                              #end=30000000,
                              snp_label_size = 1000000,                     
                              colour_scheme=c(yel,"red","blue"),
                              plot_allele_frequency=F,                    
                              invariant_filter = F)  
    
    #Change the na color to lighter gray
    gt_plot<-chr_plot$genotypes   
    updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                       breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
    chr_plot$genotypes<-updated    
    
    
    pdf(paste0("Output/chr7/DP7000/", pop, "_chr7.pdf"),width=12,height=10)
    print(combine_genotype_plot(chr_plot))
    dev.off()
    
    gt<-chr_plot$genotypes[["data"]]
    gt$Sample<-p_map$Sample
    write.csv(gt, paste0("Output/chr7/DP7000/", pops[i], "_gtinfo.csv"))
 
}


#Plot CA Allele frequency
ca<-pop_info[pop_info$Population.Year=="CA17",]
    
chr7_plot <- genotype_plot(vcf_object  =  vcf,   
                           popmap = ca,                           
                           snp_label_size = 1000000,                     
                           colour_scheme=c(yel,"red","blue"),
                           plot_allele_frequency=TRUE,                    
                           invariant_filter = F)                      
pdf("Output/chr7/chr7_CA17_AFplot.pdf", width = 10, height = 5)
combine_genotype_plot(chr7_plot)
dev.off()

af<-chr7_plot$genotypes[["data"]]
positions<-data.frame(getFIX(vcf))

af$pos<-as.integer(positions$POS)
ggplot(af, aes(x=pos, y=AF))+
    geom_point(size=0.2, color=red)

#breakpoint 1
af1<-af[af$snp_pos>18000000 & af$snp_pos<28000000,]
ggplot(af1, aes(x=snp_pos, y=AF))+
    geom_point(size=0.2, color=red)

bk<-af1[af1$AF<0.6,]
ggplot(bk, aes(x=pos, y=AF))+
    geom_point(size=0.2, color=red)
bk1<-bk[bk$pos<20000000&bk$pos>18500000,]
ggplot(bk1, aes(x=pos, y=AF))+
    geom_point(size=0.2, color=red)
bk2<-bk[bk$pos>25000000&bk$pos<27000000,]
ggplot(bk2, aes(x=pos, y=AF))+
    geom_point(size=0.2, color=red)+
    scale_x_continuous(breaks=seq(25000000,27000000, 100000))+
    theme(axis.text.x = element_text(angle=90))


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


############################
# Find individuals with inversions    

ca17<-read.csv("Output/AF/Chr7/CA17_chr7_18-28M_gtinfo.csv", row.names = 1)

#select snps in inversion regions 
#ca2<-ca17[ca17$snp>17500000&ca17$snp<25800000,]
#number of snps 
length(unique(ca17$snp))
#3063

ca_gt<-data.frame(table(ca17$GT, ca17$Sample))
ind<-unique(ca_gt$Var2)
ca_alt<-data.frame(Sample=ind)
for (i in 1: length(ind)){
    df<-ca_gt[ca_gt$Var2==ind[i],]
    ca_alt$Alt_homo[i]<-df$Freq[3]/sum(df$Freq)
}

#Plot frequency of alternative homozygous 
plot(ca_alt$Alt_homo, pch=16)
#extract the individual ids that have higher alt_homo alleles.
ca_alho<-ca_alt[ca_alt$Alt_homo>0.2,]
#38 individuals
write.csv(ca_alho,"Output/chr7/CA_samples_with_inversions_chr7.csv")


#per loci alt freq across individuals in CA pop
ca2$snp<-as.integer(ca2$snp)
snp<-data.frame(table(ca2$GT, ca2$snp))

for (i in 1: nrow(snp)){
    if (snp$Var1[i]==0) snp$percent[i]<-snp$Freq[i]/sum(snp$Freq[i:(i+2)])
    if (snp$Var1[i]==1) snp$percent[i]<-snp$Freq[i]/sum(snp$Freq[(i-1):(i+1)])
    if (snp$Var1[i]==2) snp$percent[i]<-snp$Freq[i]/sum(snp$Freq[(i-2):(i)])
    
}
snp$pos<-as.integer(as.character(snp$Var2))
ca2$snp>=17800000&ca2$snp<=17830000
ggplot(snp[snp$pos>=17600000&snp$pos<=17900000,],aes(x=Var2, y=Freq, fill=Var1))+
    geom_bar(position="fill",stat="identity")+
    theme(axis.text.x = element_text(angle=90, size=5))+
    scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                      breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))+xlab('')
ggsave("Output/AF/chr12/breakpoint1_CA.pdf", width = 8, height = 6)


ggplot(snp[snp$pos>=25700000&snp$pos<=25900000,],aes(x=Var2, y=Freq, fill=Var1))+
    geom_bar(position="fill",stat="identity")+
    theme(axis.text.x = element_text(angle=90, size=5))+
    scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                      breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
ggsave("Output/AF/chr12/breakpoint2_CA.pdf", width = 8, height = 6)





