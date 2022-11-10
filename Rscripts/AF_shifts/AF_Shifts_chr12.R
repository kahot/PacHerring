library(ggplot2)
library(tidyverse)
#library(reticulate)
library(reshape2)
#library(plyranges)
#library(seqinr)
library(ggpubr)
library(gridExtra)
library(scales)
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

#all_pops<-read.csv("Output/freqs/freq_shifts_all_populations_100kb_10kb.csv", row.names = 1)

#################################
# allele frequency shifts

pop<-"PW"
freqs_win <- read.csv(paste0("Output/freqs/PW_shifts_50kb_win_5kb_step.csv"), row.names = 1)
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

freq<-melt(freqs_win[,c("chr_num","win_mid","zt01_mean" ,"zt12_mean","zt12_mean")], id.vars=c("chr_num", "win_mid"))

##########################
## Identify the regions with big/opposite allele freq changes over time
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

#chr12
df <- freqs_win[freqs_win$chr_num ==12,]
#plot with genome position to identify the region of interest
ggplot(df, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/3)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/3)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+
    ggtitle("PWS Chr12" )+
    scale_x_continuous(breaks = seq(0, max(df$win_mid), by = 1000000),labels = comma)+
    theme(axis.text.x = element_text(angle=90))
#ggsave("Output/AF/chr7/PWS_chr7_smooth0.33.pdf", width=6, height=4)

#between 17 to 27M
af12 <- freqs_win[freqs_win$chr_num ==12 & freqs_win$start>=16000000& freqs_win$start<33000000,]

ggplot(af12, aes_string(x = "win_mid", y = "zt01_mean")) + 
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
ggsave("Output/AF/chr12/Chr12_16m-31m_afShifts.pdf", width = 8, height = 5.5)

#finer resolution plot
freqs_win2 <- read.csv(paste0("Output/freqs/PW_shifts_5kb_win_500b_step.csv"), row.names = 1)
freqs_win2$chr_num <- factor(freqs_win2$chr, levels = c(1:26))

#af7_1 <- freqs_win2[freqs_win2$chr_num ==12 & freqs_win2$win_start>=18000000& freqs_win2$win_start<28000000,]
af12_2 <- freqs_win2[freqs_win2$chr_num ==12 & freqs_win2$win_start>=19000000& freqs_win2$win_start<=20000000,]

ggplot(af12_2, aes_string(x = "win_mid", y = "zt01_mean")) + 
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
ggsave("Output/AF/chr12/Chr12_18m-16m_afShifts_5k_window.pdf", width = 8, height = 5.5)



#where is the loci?    
#1996-2007 change over 0.2
af12_d<-af12[abs(af12$zt12_mean)>0.2|abs(af12$zt01_mean)>0.2|abs(af12$zt23_mean)>0.2,]
af12_d$win_mid
# [1] 16490000 16495000 19030000 19215000 19220000 19225000 19230000 19235000 19240000 19245000
#[11] 19250000 19255000 19260000 19835000 19840000 19845000 22935000 22945000 25535000 25540000
#[21] 25545000 25550000 25555000 25895000 25900000 25905000 25910000 26550000 26745000 26750000
#[31] 26755000 26760000 26765000 26770000 26775000 26780000 26785000 26790000 26885000 26890000
#[41] 26895000 27030000 27035000 27040000 27045000 27050000 27055000 27060000 27065000 27070000
#[51] 27630000 27635000 29065000 29070000 29125000 29130000 29135000 29165000 29170000 29175000
#[61] 29180000 29185000 29370000 29375000 29435000 29440000 29445000 29450000 29455000 29460000
#[71] 29465000 29470000 29475000 29480000 29500000 29870000 29930000 29935000



#### Plot which individuals have freq changes in this region
vcf <- read.vcfR("Data/vcfs/ph_chr12.vcf.gz")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
popmap<-pop_info[c("Sample","Population.Year")]
pops<-unique(popmap$Population.Year)

#somehow loop doesn't create the correct plot output file 
#for (i in 2:length(pops)){
i=1
    print(pops[i])
    p_map<-popmap[popmap$Population.Year==pops[i],]
    chr_plot <- genotype_plot(vcf_object  =  vcf,   
                               popmap = p_map,                           
                               start=14000000,
                               end=30000000,
                               snp_label_size = 1000000,                     
                               colour_scheme=c(yel,"red","blue"),
                               plot_allele_frequency=F,                    
                               invariant_filter = F)  
    
    #Change the na color to lighter gray
    gt_plot<-chr_plot$genotypes   
    updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                       breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
    chr_plot$genotypes<-updated    
    
    
    pdf(paste0("Output/AF/chr12/", pops[i], "_chr12.pdf"),width=12,height=10)
    combine_genotype_plot(chr_plot)
    dev.off()
    
    gt<-chr_plot$genotypes[["data"]]
    gt$Sample<-p_map$Sample
    write.csv(gt, paste0("Output/AF/Chr12/", pops[i], "_gtinfo.csv"))
#}

    
    
ca17<-read.csv("Output/AF/Chr12/CA17_gtinfo.csv", row.names = 1)

#select snps in inversion regions 
ca2<-ca17[ca17$snp>17500000&ca17$snp<25800000,]
#number of snps 
length(unique(ca2$snp))
#2260

ca_gt<-data.frame(table(ca2$GT, ca2$Sample))
ind<-unique(ca_gt$Var2)
ca_alt<-data.frame(Sample=ind)
for (i in 1: length(ind)){
    df<-ca_gt[ca_gt$Var2==ind[i],]
    ca_alt$Alt_homo[i]<-df$Freq[3]/sum(df$Freq)
}

plot(ca_alt$Alt_homo)
#extract the individual ids that have higher alt_homo alleles.
ca_alho<-ca_alt[ca_alt$Alt_homo>0.2,]
#44 individuals
write.csv(ca_alho,"Output/AF/chr12/CA_samples_with_inversions.csv")


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



## Find the loci with HomAlt >0.5
loci0.5<-snp[snp$percent >=0.5 & snp$Var1==2,]


#create a scripts to extract the region from bam files

# see slurm_scipts.R extract_inversion_chr12


#cutoff =0.3
length(ca_alt$Sample[ca_alt$Alt_homo>0.3])
#43 individuals

#near brewkpoint
ggplot(ca2[ca2$snp>=17700000&ca2$snp<=17900000,], aes(x=snp, y=index))+
    geom_tile(aes(fill=factor(GT)))+
    scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                      breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))


### create a script to extract consensus for these individuals
sink("Output/Slurmscripts/extract_inversion2.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_inversion \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_inversion.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load samtools \n")  
cat("module load bcftools \n")  
cat("\n\n")


for (i in 1: nrow(ca_alho)){
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/',ca_alho[i],'.bam "chr12:17750000-25780000" >  /home/ktist/ph/data/bam/chr12/', ca_alho[i],'_chr12_inversion.bam  \n')) 
    cat(paste0('bcftools mpileup -Ou -f /home/ktist/ph/data/bam/c.harengus.fasta /home/ktist/ph/data/bam/chr12/',ca_alho[i],'_chr12_inversion.bam | bcftools call -Ou -mv | bcftools norm -f /home/ktist/ph/data/bam/c.harengus.fasta -Oz -o /home/ktist/ph/data/bam/chr12/',ca_alho[i],'_chr12_inversion.vcf.gz \n'))

    cat(paste0('tabix /home/ktist/ph/data/bam/chr12/0546_CA17.g.vcf.gz \n'))
    cat(paste0('bcftools consensus -f /home/ktist/ph/data/bam/chr12/chr12_17750000-25780000.fa /home/ktist/ph/data/bam/chr12/',ca_alho[i],'_chr12_inversion.vcf.gz > /home/ktist/ph/data/bam/chr12/',ca_alho[i],'_chr12_inversion.fa \n\n'))
}
  

### Read the allele freq of Atlantic Herring data
keypops<-read.table("Data/chr12/Han_invert_pops.txt", header = F)
#Afreq<-read.table("Data/chr12/chr12.60Neff.freq", sep="\t")
#head<-read.table("Data/chr12/Neff.small", sep="\t", header=T)
#colnames(Afreq)<-colnames(head)
#keep<-which(colnames(Afreq) %in% keypops$V1)
#write.csv(Afreq, "Data/chr12/hr12.60Neff.freq.csv")

Afreq<-read.csv("Data/chr12/chr12.han_60Neff.freq.csv", row.names = 1)
keypops<-read.table("Data/chr12/Han_invert_pops.txt", header = F)
keep<-which(colnames(Afreq) %in% keypops$V1)
Afreq1<-Afreq[, c(1,2,keep)]

#select positions in inversion region
Afreq2<-Afreq1[Afreq1$POS>=17750000&Afreq1$POS<=25780000,]
Afreq2

#any overlapping positions?
#create freq file for subset (w/ inversions) of CA17 pop


ca<-read.table("Data/chr12/CA17_inversion_freq.frq",sep="\t", row.names=NULL)
colnames(ca)[1:4]<-colnames(ca)[2:5]
colnames(ca)[5:6]<-c("f1","f2")
ca$pos<-paste0(ca$CHROM,":",ca$POS)
ca$nuc1<-substr(ca$f1,1,1)
ca$nuc2<-substr(ca$f2,1,1)
ca$freq1<-as.numeric(substr(ca$f1,3,10))
ca$freq2<-as.numeric(substr(ca$f2,3,10))

ca<-ca[ca$CHROM=="chr12",]
write.csv(ca,"Output/chr12/CA17_inversion_chr12_freq.csv")
    
#select the inversion region
ca2<-ca[ca$POS>=17750000&ca$POS<=25780000,]
ca2<-ca2[,c("POS","freq2")]
colnames(ca2)[2]<-"CA17"
m1<-merge(Afreq2, ca2, by="POS")
m2<-merge(Afreq2, ca2, by="POS", all.y = T)

m1$Atl.mean<-rowMeans(m1[,3:10], na.rm=T)
m1$diff<-m1$Atl.mean-m1$CA17    

#2019 sites overlapping from Han (2020)
library(ggplot2)
ggplot()+
    geom_point(data=m1, aes(x=POS, y=Atl.mean), color="dodgerblue", size=1, alpha=0.6)+
    geom_point(data=m1, aes(x=POS, y=CA17), color="indianred1", size=1, alpha=0.6)
    geom_hline(yintercept=0.629, color="gray60")

#Plot the difference:
#Same inverted sites should be around ~0.4
#Negative sites -only distinct in CA pop
    
m1$overlap<-'Some overlap'
m1$overlap[m1$diff<0.25&m1$diff>=0&m1$Atl.mean>0.75]<-"Overlapped"
m1$overlap[m1$diff>=0.25&m1$Atl.mean>0.75&m1$CA17<0.2]<-"Atlantic loci"
m1$overlap[m1$diff<0]<-"CA loci"

m1$position<-as.character(m1$POS)
library(scales)

ggplot()+
    geom_point(data=m1,aes(x=POS, y=diff, color=overlap),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)
ggsave("Output/chr12/DifferenceinAF_AtlanticvsCA_inversion.pdf", width = 9, height = 6)
    
table(m1$overlap)
#Atlantic loci       CA loci    Overlapped  Some overlap 
#          482           428           944           165 



# look at maf00 file
ca0<-read.table("Data/chr12/CA17_inversion_chr12_maf00.freq",sep="\t", row.names=NULL)
colnames(ca0)<-colnames(ca)[1:6]

ca0$pos<-paste0(ca0$CHROM,":",ca0$POS)
ca0$nuc1<-substr(ca0$f1,1,1)
ca0$nuc2<-substr(ca0$f2,1,1)
ca0$freq1<-as.numeric(substr(ca0$f1,3,10))
ca0$freq2<-as.numeric(substr(ca0$f2,3,10))

#select the inversion region
ca00<-ca0[ca0$POS>=17750000&ca0$POS<=25780000,]
ca00<-ca00[,c("POS","freq2")]
colnames(ca00)[2]<-"CA17"
m01<-merge(Afreq2, ca00, by="POS") #6282 overlapping loci

m01$Atl.mean<-rowMeans(m01[,3:10], na.rm=T)
m01$diff<-m01$Atl.mean-m01$CA17    

ggplot()+
    geom_point(data=m01, aes(x=POS, y=Atl.mean), color="dodgerblue", size=0.3, alpha=0.6)+
    geom_point(data=m01, aes(x=POS, y=CA17), color="indianred1", size=0.3, alpha=0.6)
geom_hline(yintercept=0.629, color="gray60")

#Plot the difference:
#Same inverted sites should be around ~0.4
#Negative sites -only distinct in CA pop

m01$overlap<-'Some overlap'
m01$overlap[m01$diff<0.25&m01$diff>=0&m01$Atl.mean>0.75]<-"Overlapped"
m01$overlap[m01$diff>=0.25&m01$Atl.mean>0.75&m01$CA17<0.2]<-"Atlantic loci"
m01$overlap[m01$diff<0]<-"CA loci"

ggplot()+
    geom_point(data=m01,aes(x=POS, y=diff, color=overlap),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)
ggsave("Output/chr12/DifferenceinAF.maf00_AtlanticvsCA_inversion.pdf", width = 9, height = 6)

table(m01$overlap)
#Atlantic loci       CA loci    Overlapped  Some overlap 
#         4060           433           957           832

#more fixed differences in Atlantic populations 

summ<-data.frame(table(m01$overlap))
prop.table(summ$Freq)
# 0.64629099 0.06892709 0.15234002 0.13244190
write.csv(m01,"Output/chr12/chr12_inversion_overlappedLoci_withAtlanticHerring.csv")

## created a new vcf with no max depth (min 100) for CA inversion subpopulation
# created freq
# look at maf00 file
ca0<-read.table("Data/chr12/CA17_inversion_freq_filtered_MD100.frq",sep="\t", row.names=NULL)
colnames(ca0)<-c("CHROM","POS","N_ALLELES","N_CHR","f1","f2")

ca0$pos<-paste0(ca0$CHROM,":",ca0$POS)
ca0$nuc1<-substr(ca0$f1,1,1)
ca0$nuc2<-substr(ca0$f2,1,1)
ca0$freq1<-as.numeric(substr(ca0$f1,3,10))
ca0$freq2<-as.numeric(substr(ca0$f2,3,10))

#select the inversion region
ca00<-ca0[ca0$POS>=17750000&ca0$POS<=25780000,]
ca00<-ca00[,c("POS","freq2")]
colnames(ca00)[2]<-"CA17"
m01<-merge(Afreq2, ca00, by="POS") #9112 overlapping loci

m01$Atl.mean<-rowMeans(m01[,3:10], na.rm=T)
m01$diff<-m01$Atl.mean-m01$CA17    

ggplot()+
    geom_point(data=m01, aes(x=POS, y=Atl.mean), color="dodgerblue", size=0.3, alpha=0.6)+
    geom_point(data=m01, aes(x=POS, y=CA17), color="indianred1", size=0.3, alpha=0.6)


#Plot the difference:
#Same inverted sites should be around ~0.4
#Negative sites -only distinct in CA pop

m01$overlap<-'Some overlap'
m01$overlap[m01$diff<0.25&m01$diff>=0&m01$Atl.mean>0.75]<-"Overlapped"
m01$overlap[m01$diff>=0.25&m01$Atl.mean>0.75&m01$CA17<0.2]<-"Atlantic loci"
m01$overlap[m01$diff<0]<-"CA loci"

ggplot()+
    geom_point(data=m01,aes(x=POS, y=diff, color=overlap),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)
ggsave("Output/chr12/DifferenceinAF.maf00_AtlanticvsCA_inversion.pdf", width = 9, height = 6)

summ<-data.frame(table(m01$overlap))
prop.table(summ$Freq)

#Atlantic loci       CA loci    Overlapped  Some overlap 
#        4568          2136          1259          1149 
#   0.5013169     0.2344162      0.1381694     0.1260975



#Compare the rest of the chr12
ca02<-read.table("Data/chr12/CA17_inversion_freq_filtered_MD100_2.frq",sep="\t", row.names=NULL)
colnames(ca02)<-c("CHROM","POS","N_ALLELES","N_CHR","f1","f2")

ca02$pos<-paste0(ca02$CHROM,":",ca02$POS)
ca02$nuc1<-substr(ca02$f1,1,1)
ca02$nuc2<-substr(ca02$f2,1,1)
ca02$freq1<-as.numeric(substr(ca02$f1,3,10))
ca02$freq2<-as.numeric(substr(ca02$f2,3,10))

ca3<-ca0[ca0$POS<17600000,]
ca3<-rbind(ca02, ca3)

ca3<-ca3[,c("POS","freq2")]
colnames(ca3)[2]<-"CA17"

Afreq3<-Afreq[, c(1,2,keep)]
Afreq3<-Afreq3[Afreq3$POS<17600000,]

m3<-merge(Afreq3, ca3, by="POS") #14010 overlapping loci
m3$Atl.mean<-rowMeans(m3[,3:10], na.rm=T)
m3$diff<-m3$Atl.mean-m3$CA17    

m3$overlap<-'Some overlap'
m3$overlap[m3$diff<0.25 &m3$diff>=0&m3$Atl.mean>0.75]<-"Overlapped"
m3$overlap[m3$diff>=0.25&m3$Atl.mean>0.75&m3$CA17<0.2]<-"Atlantic loci"
m3$overlap[m3$diff<0]<-"CA loci"

ggplot()+
    geom_point(data=m3,aes(x=POS, y=diff, color=overlap),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)

prop.table(table(m3$overlap))

#Atlantic loci       CA loci    Overlapped  Some overlap 
#         6833          2176          2965          2036 
#    0.4877231     0.1553176     0.2116345     0.1453248

#Categorize loci between CA and Atlantic

#Calcualte Fst between populations, nuc diversity 

#noninversion region of chr12
m3$group<-'?'
m3$group[m3$diff<0.25 &m3$diff>(-0.25)&m3$Atl.mean<0.5&m3$CA17<0.5]<-"Both ref"
m3$group[m3$diff<0.25 &m3$diff>(-0.25)&m3$Atl.mean>=0.5&m3$CA17>=0.5]<-"Both alt"
m3$group[m3$diff>=0.25&m3$Atl.mean>=0.5&m3$CA17<0.5]<-"Atlantic alt"
m3$group[m3$diff<0&m3$Atl.mean<0.5&m3$CA17>=0.5]<-"CA alt"
table(m3$group)
ggplot()+
    geom_point(data=m01,aes(x=POS, y=diff, color=group),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange","pink"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)+
    ggtitle('chr12: non-inversion region')
ggsave("Output/chr12/AFcomparison_Non-inversion__Atlanticvs.CA.pdf", width=8.5, height = 6)


m01$group<-'?'
m01$group[m01$diff<0.25&m01$diff>(-0.25)&m01$Atl.mean<0.5&m01$CA17<0.5]<-"Both ref"
m01$group[m01$diff<0.25&m01$diff>(-0.25)&m01$Atl.mean>=0.5&m01$CA17>=0.5]<-"Both alt"
m01$group[m01$diff>=0.25&m01$Atl.mean>=0.5&m01$CA17<0.5]<-"Atlantic alt"
m01$group[m01$diff<=(-0.25)&m01$Atl.mean<0.5&m01$CA17>=0.5]<-"CA alt"
table(m01$group)

ggplot()+
    geom_point(data=m01,aes(x=POS, y=diff, color=group),size=0.5, alpha=0.6)+
    scale_color_manual(values=c("gray50", "lightblue","red","orange","pink"))+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle=50, hjust=1))+ylab("Allele freq difference")+
    xlab("chr12 position")+
    scale_x_continuous(labels = comma)
ggsave("Output/chr12/AFcomparison_Inversion_Atlanticvs.CA.pdf", width=8.5, height = 6)

sum<-data.frame(table(m01$group))
colnames(sum)<-c("group","inversion")
sum2<-data.frame(table(m3$group))
colnames(sum2)<-c("group","non-inversion")
sum<-merge(sum, sum2, by="group")
ratio<-sum[,2:3]
rownames(ratio)<-sum[,1]
#run chi-square test
chisq.test(ratio)

m01$Atl.maf=apply(m01["Atl.mean"], 1, function(x) ifelse(x>0.5, 1-x, x))
m01$CA.maf=apply(m01["CA17"], 1, function(x) ifelse(x>0.5, 1-x, x))
m3$Atl.maf=apply(m3["Atl.mean"], 1, function(x) ifelse(x>0.5, 1-x, x))
m3$CA.maf=apply(m3["CA17"], 1, function(x) ifelse(x>0.5, 1-x, x))

mean(m01$Atl.maf, na.rm=T)
mean(m01$CA.maf, na.rm=T)
mean(m3$Atl.maf, na.rm=T)
mean(m3$CA.maf, na.rm=T)

maf<-data.frame(loci=rep(c("Inversion", "Non-inveresion"), times=2), Population=rep(c("Atlantic","CA"), each=2))
maf$maf[1]<-mean(m01$Atl.maf, na.rm=T)
maf$maf[2]<-mean(m3$Atl.maf, na.rm=T)
maf$maf[3]<-mean(m01$CA.maf, na.rm=T)
maf$maf[4]<-mean(m3$CA.maf, na.rm=T)

ggplot(maf, aes(x=loci, y=maf, fill=Population))+
    geom_bar(stat="identity", position=position_dodge())+
    theme_classic()+
    xlab('')+ylab("MAF")+
    scale_fill_manual(values=c(grb,org))
ggsave("Output/chr12/meanMAF_comparison.pdf", width = 6, height = 4)
