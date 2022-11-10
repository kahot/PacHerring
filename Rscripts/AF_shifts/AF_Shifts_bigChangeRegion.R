library(ggplot2)
library(tidyverse)
library(reticulate)
library(reshape2)
library(plyranges)
library(seqinr)
library(windowscanr)
library(ggpubr)
library(venn)
library(gridExtra)

bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'

pops<-c("PW","SS","TB")

win <- "50kb"
step <- "5kb"

all_pops<-read.csv("Output/freqs/freq_shifts_all_populations_100kb_10kb.csv", row.names = 1)

#################################
# allele frequency shifts

pop<-"PW"
freqs_win <- read.csv(paste0("Output/freqs/PW_shifts_50kb_win_5kb_step.csv"), row.names = 1)
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

freq<-melt(freqs_win[,c("chr_num","win_mid","zt01_mean" ,"zt12_mean","zt12_mean")], id.vars=c("chr_num", "win_mid"))

#ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
#    theme_minimal()+
#    theme(axis.text.x=element_blank())+
#    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
#    geom_line(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
#    geom_line(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
#    geom_line(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
#    facet_wrap(~chr_num, ncol = 9)+
#    scale_color_manual(values=c(lir,grb,org),name = "contrast",
#                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
#    ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
#    ggtitle(pop)
#

ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/5)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/5)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/5)+
    facet_wrap(~chr_num, ncol = 9)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
    ggtitle(pop)

# Chr 7, 12,16,21, 9 shows fluctuation





## Within population chromosomes with putative inversions (7,12,15,16,20)

# chr 7
pop <- "PW"
chrs <- c(4,7,12,16,18,21)
win= "50kb"
step="5kb"

df <- freqs_win[freqs_win$chr_num %in% chrs,]
ggplot(df, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/3)+
    facet_wrap(~chr_num, ncol = 3)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+
    ggtitle(pop)
ggsave("Output/AF/Chr7/PWS_fluctuatingChr4.7.12.16.18.21.smooth0.33.pdf", width = 8, height = 5)


##########################
## Identify the regions with big/opposite allele freq changes over time
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

#start with chr 7
df <- freqs_win[freqs_win$chr_num ==7,]
require(scales)
ggplot(df, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/3)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+
    ggtitle(pop)+
    scale_x_continuous(breaks = seq(0, max(df$win_mid), by = 1000000),labels = comma)+
    theme(axis.text.x = element_text(angle=45))
ggsave("Output/AF/chr7/PWS_chr7_smooth0.33.pdf", width=6, height=4)
#between 18 to 28

af7 <- freqs_win[freqs_win$chr_num ==7 & freqs_win$start>=16000000& freqs_win$start<30000000,]

ggplot(af7, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"))+
    geom_line(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"))+
    geom_line(aes_string(x = "win_mid", y = "zt23_mean", color = "org"))+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    #ylim(-0.05, 0.05)+
    ggtitle(pop)+
    scale_x_continuous(breaks = seq(0, max(df$win_mid), by = 1000000),labels = comma)+
    theme(axis.text.x = element_text(angle=45))
ggsave("Output/AF/chr7/Chr7_16m-30m_afShifts.pdf", width = 8, height = 5.5)

#finer resolution plot
freqs_win2 <- read.csv(paste0("Output/freqs/PW_shifts_5kb_win_500b_step.csv"), row.names = 1)
freqs_win2$chr_num <- factor(freqs_win2$chr, levels = c(1:26))

af7_1 <- freqs_win2[freqs_win2$chr_num ==7 & freqs_win2$win_start>=18000000& freqs_win2$win_start<28000000,]
af7_2 <- freqs_win2[freqs_win2$chr_num ==7 & freqs_win2$win_start>=18000000& freqs_win2$win_start<=19000000,]

ggplot(af7_2, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"))+
    geom_line(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"))+
    geom_line(aes_string(x = "win_mid", y = "zt23_mean", color = "org"))+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ggtitle(pop)+
    scale_x_continuous(breaks = seq(0, max(df$win_mid), by = 1000000),labels = comma)+
    theme(axis.text.x = element_text(angle=45))
ggsave("Output/AF/chr7//Chr7_18m-16m_afShifts_5k_window.pdf", width = 8, height = 5.5)


ggplot(af7_1, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"))+
    geom_line(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"))+
    geom_line(aes_string(x = "win_mid", y = "zt23_mean", color = "org"))+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ggtitle(pop)+
    scale_x_continuous(breaks = seq(0, max(df$win_mid), by = 1000000),labels = comma)+
    theme(axis.text.x = element_text(angle=45))


#where is the loci?    
#1996-2007 change over 0.2
af7_d<-af7[af7$zt12_mean>0.3,]
af7_d$win_mid
#[1] 18715000 18720000

#find annotations for the region
# write a bed file

srt<-af7_d$start[1]-50000
fin<-af7_d$end[2]+50000
sink("Output/freqs/PWS/chr7.bed")
cat(paste("chr7", "18500000","19000000", sep="\t"))
sink(NULL)


#PWS17

# Create a new vcf file containing the loci in p07_loc.
## at terminal Output/freqs/PWS/chr7.bed
#vcftools --gzvcf Data/vcfs/ph_chr7.vcf.gz --bed Output/freqs/PWS/chr7.bed --out Output/AF/chr7/ch7_bigchang --recode --keep-INFO-all

tabix -p vcf Data/vcfs/ph_chr7.vcf.gz
bcftools view Data/vcfs/ph_chr7.vcf.gz -r chr7:18650000-18790000 -Oz > Output/AF/chr7/chr7_bigchang_loci.vcf.gz

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/AF/chr7/chr7_bigchang_loci.vcf -stats ~/Projects/PacHerring/Output/AF/chr7/chr7_loci > ~/Projects/PacHerring/Output/AF/chr7/Anno.chr_bigchange_loci.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS07_PBS_outlier.vcf > PWS07_annotation


# for larger region of chr7
bcftools view Data/vcfs/ph_chr7.vcf.gz -r chr7:18000000-28000000 -Oz > Output/AF/chr7/chr7_bigchang_loci2.vcf.gz
bgzip -d Output/AF/chr7/chr7_bigchang_loci2.vcf.gz
#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/AF/chr7/chr7_bigchang_loci2.vcf -stats ~/Projects/PacHerring/Output/AF/chr7/chr7_loci2 > ~/Projects/PacHerring/Output/AF/chr7/Anno.chr_bigchange_loci2.vcf
#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/AF/chr7/Anno.chr_bigchange_loci2.vcf > Output/AF/chr7/PWS_chr7_18-28m_annotation

## Annotation summary
#read the annotation info
df2<-read.table("Output/AF/Chr7/chr7_18-28m_annotation", header = F)
#remove the duplicates
annotations2<-data.frame()
for (i in 1: nrow(df2)){
    anns<-unlist(strsplit(df2$V4[i], "\\|"))
    anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
    annotations2<-rbind(annotations2, anns)
}      

colnames(annotations2)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
Ano2<-cbind(df2[,1:3], annotations2)
colnames(Ano2)[1:3]<-c("chr","pos","AF")
write.csv(Ano2, "Output/AF/Chr7/chr7_18-28m_extractedAnnotation_info.csv")            

genes<-unique(c(Ano2$Gene_name, Ano2$Gene_name2))
#437 genes involved in this region

genes2<-genes[!grepl("ENSCHAG0", genes)] # annotated genes
#232 annotated genes

p<-"ENSCHAG00000018984"
nchar(p)

ids<-unique(c(Ano2$Gene_ID, Ano2$Gene_ID2))
ids<-substr(ids, 1,18)
ids<-unique(ids)
ids<-ids[!is.na(ids)] #329 unique gene IDs


sink("Output/AF/Chr7/Region18_28m_geneID.txt")
for (i in 1: length(ids)){
    cat(paste0(ids[i],"; "))
}
sink(NULL)



## Look at the allele frequency of each SNP
##################
library(vcfR)

## Look at the maf005 file for freq change
pwspop<-c("PWS91","PWS96","PWS07","PWS17")
year<-c(1991,1996,2007,2017)

AF<-data.frame()
for (i in 1:4){
    vcf <- read.vcfR(paste0("Output/AF/Chr7/", pwspop[i],"_chr7_region1.vcf.gz"))
    
    #extract the allele freq
    af<-data.frame(maf(vcf, element=2))
    
    positions<-as.integer(gsub("chr7_", '', rownames(af))) #3966
    
    af$pos<-positions
    af$year<-year[i]
    
    AF<-rbind(AF, af[,c("pos","Frequency","year")])
}



subAF<-AF[AF$pos>18650000 & AF$pos<18790000,]
posi<-unique(subAF$pos) #24 loci

#plot 5 at a time 
vec<-seq(1,24,5)
plots<-list()
for (i in 1:length(vec)){
    df<-subAF[subAF$pos %in% posi[vec[i]:(vec[i]+4)],]
    plots[[i]]<-ggplot(df, aes(x=year, y=Frequency, color=factor(pos), group=factor(pos)))+
        geom_point()+
        geom_line ()+
        #ylim(0,(max(df$Frequency)+0.01))+
        theme_classic()+
        labs(color="Chr7 pos")+
        ylim(0,0.5)
    
}
pdf(paste0("Output/AF/chr7/AF_change_maf05_scaled.pdf"), width = 14, height = 7)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()

## annotation information

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/AF/chr7/Anno.chr_bigchange_loci.vcf > Output/AF/chr7/chr7_region1_annotation

#read the annotation info
df<-read.table("Output/AF/Chr7/chr7_region1_annotation", header = F)

annotations<-data.frame()
for (i in 1: nrow(df)){
    anns<-unlist(strsplit(df$V4[i], "\\|"))
    anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
    annotations<-rbind(annotations, anns)
}      

colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
Ano<-cbind(df[,1:3], annotations)
colnames(Ano)[1:3]<-c("chr","pos","AF")
write.csv(Ano, "Output/AF/Chr7/chr7_region1_extractedAnnotation_info.csv")            


###
## Look at the maf00 file for freq change
pwspop<-c("PWS91","PWS96","PWS07","PWS17")
year<-c(1991,1996,2007,2017)

AF<-data.frame()
for (i in 1:4){
    vcf <- read.vcfR(paste0("Output/AF/Chr7/maf00/", pwspop[i],"_chr7_region1.vcf.gz"))
    
    #extract the allele freq
    af<-data.frame(maf(vcf, element=2))
    
    positions<-as.integer(gsub("chr7_", '', rownames(af))) #3966
    
    af$pos<-positions
    af$year<-year[i]
    
    AF<-rbind(AF, af[,c("pos","Frequency","year")])
}

subAF<-AF[AF$pos>18650000 & AF$pos<18790000,]
posi<-unique(subAF$pos) #1396 loci

#plot 10 at a time 
vec<-seq(1,length(posi),20)
plots2<-list()
for (i in 1:length(vec)){
    df<-subAF[subAF$pos %in% posi[vec[i]:(vec[i]+19)],]
    plots2[[i]]<-ggplot(df, aes(x=year, y=Frequency, color=factor(pos), group=factor(pos)))+
        geom_point()+
        geom_line ()+
        ylim(0,(max(df$Frequency)+0.01))+
        theme_classic()+
        theme(legend.title = element_blank())+
        guides(color=guide_legend(ncol=2))
}

pvec<-seq(1, length(plots2), 18)
for (i in 1: length(pvec)){
    pdf(paste0("Output/AF/chr7/maf00/AF_change",i,".pdf"), width = 24, height = 16)
    if (i==4) do.call(grid.arrange, c(plots2[pvec[i]:(pvec[i]+15)], ncol=3))
    else do.call(grid.arrange, c(plots2[pvec[i]:(pvec[i]+17)], ncol=3))
    dev.off()
}





#### Plot which individuals have freq changes in this region
library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)

pws_vcf <- read.vcfR("Data/vcfs/PWS_SS/PWS_SS_chr7.vcf.gz")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws_info<-pop_info[pop_info$pop=="PWS",]
pws_popmap<-pws_info[c("Sample","Population.Year")]

#meta <- read.csv('Data/vcfs/sample_metadata.txt', comment.char = '#', sep = "\t")
#my_popmap <- meta[c("Sample","Population_Year")]

chr7_plot <- genotype_plot(vcf_object  =  pws_vcf,   
                          popmap = pws_popmap,                           
                          snp_label_size = 1000000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=F,                    
                          invariant_filter = TRUE)                      

pdf(paste0("Output/AF/Chr7/chr7_SNPs_PWS.pdf"),width=10,height=8)
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



