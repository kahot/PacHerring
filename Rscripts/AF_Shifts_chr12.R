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
theme_bw()+
    