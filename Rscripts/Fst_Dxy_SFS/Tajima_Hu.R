## Plot Pi and theta estimated from downsampled SFS using bam files (41 individuals, 100 max depth)  -> will chage  
library(gridExtra)
library(windowscanr)
source("Rscripts/BaseScripts.R")

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)

cols<-qualitative_hcl(5, palette="Dark3")



######  Tajima's D and Fay's H comparison between SFS from Joe and SFS using maf00 and all samples and all sites   ########

# Fay & Wu's H
#Distinguish between a DNA sequence evolving randomly ("neutrally") and one evolving under positive selection. 
#This test is an advancement over Tajima's D. Fay and Wu's H is frequently used to identify sequences which have experienced selective sweeps in their evolutionary history.
#A significantly positive Fay and Wu's H indicates a deficit of moderate- and high-frequency derived single nucleotide polymorphisms (SNPs) relative to equilibrium expectations, 
#  whereas a significant negative Fay and Wu's H indicates an excess of high-frequency derived SNPs

pops<-c("PWS96","PWS07","PWS17")
pops2<-c("PWS91","PWS96","PWS07","PWS17")

names<-c("Downsampled","AllSamples", "maf00")

#1. unfolded -for comparison & Fay's H only 
old<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/Theta/unfolded/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    #theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","fayh" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    old<-rbind(old, df1)
}

new<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromBam/unfolded/',pops[i],'_50kwin_10kstep.pestPG'))
    df2<-theta2[,c("Chr","WinCenter","fayh" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    new<-rbind(new, df2)
}

maf0<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromVCF/',pops[i],'_maf00.thetas50kWindow.gz.pestPG'))
    df3<-theta2[,c("Chr","WinCenter","fayh" )]
    df3$pop<-pops[i]    
    df3$data<-names[3]
    maf0<-rbind(maf0, df3)
}

fayH<-rbind(old, new, maf0)    
aggregate(fayH$fayh, by=list(fayH$pop,fayH$data), mean, na.rm=T) 
#1. unfolded
#   Group.1     Group.2            x
#1    PWS07  AllSamples  0.016625489
#2    PWS17  AllSamples  0.004191033
#3    PWS96  AllSamples -0.003257172
#4    PWS07 Downsampled  0.240189072
#5    PWS17 Downsampled  0.264787545
#6    PWS91 Downsampled  0.256330603
#7    PWS96 Downsampled  0.259873412
#8    PWS07       maf00  0.096787339
#9    PWS17       maf00  0.167483837
#10   PWS96       maf00  0.103264660

fayH<-fayH[fayH$pop!="PWS91",]
fayH$pop<-factor(fayH$pop, levels=pops)
ggplot(fayH, aes(x=pop, y=fayh, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab("Fay's H")+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/FayH_estimates_comparison_bam.vcf.pdf", width = 6, height = 4.5)

ggplot(fayH, aes(x=data, y=fayh, color=pop))+
    facet_wrap(~Chr, ncol=6)+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab("Fay's H")+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/FaysH_estimates_comparison_byChrom_groupedbymethod.pdf", width = 18, height = 13)




#2. folded -for comparison of Dajima's D
old2<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    #theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","Tajima" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    old2<-rbind(old2, df1)
}

new2<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromBam/folded/',pops[i],'_50kwin_10kstep.pestPG'))
    df2<-theta2[,c("Chr","WinCenter","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    new2<-rbind(new2, df2)
}

maf02<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromVCF/folded/folded_',pops[i],'_maf00.thetas50kWindow.gz.pestPG'))
    df3<-theta2[,c("Chr","WinCenter","Tajima" )]
    df3$pop<-pops[i]    
    df3$data<-names[3]
    maf02<-rbind(maf02, df3)
}

tajima<-rbind(old2, new2, maf02)    
aggregate(tajima$Tajima, by=list(tajima$pop,tajima$data), mean, na.rm=T) 

#  Group.1     Group.2          x
#1   PWS07  AllSamples -1.9217199
#2   PWS17  AllSamples -1.8698115
#3   PWS96  AllSamples -1.9611638
#4   PWS07 Downsampled -1.3966960
#5   PWS17 Downsampled -1.1368737
#6   PWS96 Downsampled -1.1838244
#7   PWS07       maf00 -0.9577346
#8   PWS17       maf00 -1.5043204
#9   PWS96       maf00 -1.4703454



# Tajima's D for folded
tajima$pop<-factor(tajima$pop, levels=pops)
ggplot(tajima, aes(x=pop, y=Tajima, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab("Tajima's D")+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/TajimaD_estimates_comparison_bam.vcf.pdf", width = 6, height = 4.5)

#by chromosome
ggplot(tajima, aes(x=pop, y=Tajima, color=data))+
    facet_wrap(~Chr, ncol=6)+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab("Tajima's D")+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/TajimaD_estimates_comparison_byChrom.pdf", width = 6, height = 4.5)


ggplot(tajima, aes(x=data, y=Tajima, color=pop))+
    facet_wrap(~Chr, ncol=6)+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab("Tajima's D")+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/TajimaD_estimates_comparison_byChrom_groupedbymethod.pdf", width = 18, height = 13)



