## Plot Pi and theta estimated from downsampled SFS using bam files (41 individuals, 10M max depth) 
library(gridExtra)
library(windowscanr)
source("Rscripts/BaseScripts.R")

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)

Plots<-list()
for (i in 1:length(pops)){
    theta<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta$pi<-theta$tP/theta$nSites
    theta$ch<-as.integer(gsub("chr", "", theta$Chr))
    Plots[[i]]<-ggplot(theta, aes(x = WinCenter, y = pi)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(pi))+ xlab("")+ 
        ggtitle(pops[i])+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = "steelblue")+
        facet_wrap(~ch, ncol = 9)
    
}

png("Output/Pi/Pi_all_pops.png", height = 15, width = 20, res=150, units = "in")
do.call(grid.arrange, c(Plots, ncol=4))
dev.off()

#Line plots
Plots<-list()
thetas<-data.frame()
for (i in 1:length(pops)){
    theta<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta$pi<-theta$tP/theta$nSites
    theta$ch<-as.integer(gsub("chr", "", theta$Chr))
    theta$pop<-pops[i]
    thetas<-rbind(thetas, theta[,c(2,3,15,15,16,17)])
    Plots[[i]]<-ggplot(theta, aes(x = WinCenter, y = pi)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+ylim(0,0.018)+
        theme(axis.text.x=element_blank())+
        ylab(expression(pi))+ xlab("")+ 
        ggtitle(pops[i])+
        geom_line(color = "steelblue", size=0.1)+
        facet_wrap(~ch, ncol = 9)
    
}
write.csv(thetas, "Output/Pi/Pi_windows_all_pop.csv")

#thetas<-read.csv("Output/Pi/Pi_windows_all_pop.csv", row.names = 1)

png("Output/Pi/Pi_all_pops_line.png", height = 15, width = 20, res=150, units = "in")
do.call(grid.arrange, c(Plots, ncol=4))
dev.off()

#PWS
png("Output/Pi/Pi_PWS_line.png", height = 8, width = 18, res=150, units = "in")
grid.arrange(Plots[[5]],Plots[[6]],Plots[[3]],Plots[[4]], ncol=2)
dev.off()

#TB
png("Output/Pi/Pi_TB_line.png", height = 8, width = 12, res=150, units = "in")
grid.arrange(Plots[[12]],Plots[[13]],Plots[[10]],Plots[[11]], ncol=2)
dev.off()

#SS
png("Output/Pi/Pi_SS_line.png", height = 4, width = 18, res=150, units = "in")
grid.arrange(Plots[[9]],Plots[[7]],Plots[[8]], ncol=3)
dev.off()

#Y2017
png("Output/Pi/Pi_2017_line.png", height = 8, width = 18, res=150, units = "in")
grid.arrange(Plots[[11]],Plots[[4]],Plots[[8]],Plots[[1]],Plots[[14]],Plots[[2]], ncol=3)
dev.off()


#Calculate means for each chromosome
pi<-data.frame(chr=1:26)
for (i in 1:length(pops)){
    theta<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta$pi<-theta$tP/theta$nSites
    theta$ch<-as.integer(gsub("chr", "", theta$Chr))
    meanPi<-aggregate(theta$pi, by=list(theta$ch), mean, na.rm=T)
    pi[, paste0(pops[i])]<-meanPi[,2]
    
}
write.csv(pi, "Output/Pi/mean_pi_byChrom.csv")

library(reshape2)
library(colorspace)
cols<-qualitative_hcl(5, palette="Dark3")
pim<-melt(pi, id.vars="chr", value.name = "pi")
pim$chr<-factor(pim$chr, levels=1:26)
ggplot(pim, aes(x=chr, y=pi, color=variable))+
    geom_point(position=position_dodge(width = 0.5))+
    theme_bw()+
    geom_path(aes(x=chr, y=pi, color=variable, group=factor(variable)), position=position_dodge(width=0.5))+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())
ggsave("Output/Pi/Pi_mean_bychrom_bypop.pdf", width=10, height=5)

#Each group
pipws<-pim[grep("PWS",pim$variable),]
pipws$variable<-factor(pipws$variable, levels=paste0("PWS",c(91, 96, "07", 17)))
ggplot(pipws, aes(x=chr, y=pi, color=variable))+
    geom_point(position=position_dodge(width = 0.5))+
    theme_bw()+scale_color_manual(values=cols[c(1,2,3,5)])+
    geom_path(aes(x=chr, y=pi, color=variable, group=factor(variable)), position=position_dodge(width=0.5))+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())
ggsave("Output/Pi/Pi_mean_bychrom_PWS.pdf", width=10, height=5)

pitb<-pim[grep("TB",pim$variable),]
pitb$variable<-factor(pitb$variable, levels=paste0("TB",c(91, 96, "06", 17)))
ggplot(pitb, aes(x=chr, y=pi, color=variable))+
    geom_point(position=position_dodge(width = 0.5))+
    theme_bw()+scale_color_manual(values=cols[c(1,2,3,5)])+
    geom_path(aes(x=chr, y=pi, color=variable, group=factor(variable)), position=position_dodge(width=0.5))+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())
ggsave("Output/Pi/Pi_mean_bychromTB.pdf", width=10, height=5)

piss<-pim[grep("SS",pim$variable),]
piss$variable<-factor(piss$variable, levels=paste0("SS",c(91,96, "06", 17)))

ggplot(piss, aes(x=chr, y=pi, color=variable))+
    geom_point(position=position_dodge(width = 0.5))+
    theme_bw()+scale_color_manual(values=cols[c(2,3,5)])+
    geom_path(aes(x=chr, y=pi, color=variable, group=factor(variable)), position=position_dodge(width=0.5))+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())
ggsave("Output/Pi/Pi_mean_bychrom_SS.pdf", width=10, height=5)


pi17<-pim[grep("17",pim$variable),]
pi17$variable<-factor(pi17$variable, levels=paste0(c("TB","PWS","SS","BC","WA","CA"),17))
hcl_palettes(plot = TRUE)
cols2<-sequential_hcl(7, palette="Heat2")
cols3<-diverging_hcl(6, palette="Blue-Red2")

ggplot(pi17, aes(x=chr, y=pi, color=variable))+
    geom_point(position=position_dodge(width = 0.5))+
    theme_bw()+scale_color_manual(values=cols3)+
    geom_path(aes(x=chr, y=pi, color=variable, group=factor(variable)), position=position_dodge(width=0.5))+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())
ggsave("Output/Pi/Pi_mean_bychrom_2017.pdf", width=10, height=5)


## plot in violin 
#pop mean
Pi.mean<-data.frame(aggregate(thetas[,"pi"], by=list(thetas$pop), mean, na.rm=T))
colnames(Pi.mean)<-c("pop","mean")

thetas$Group<-substr(thetas$pop, 1, 2)
thetas$Group<-factor(thetas$Group, levels=c("TB","PW","SS","BC","WA","CA"))
thetas$pop<-factor(thetas$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))
Pi.mean$pop<-factor(Pi.mean$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))

ggplot() +
    geom_violin(data=thetas, aes(x=pop, y=pi, fill=Group)) + xlab("\npopulation/year") + ylab(expression(pi))+
    scale_fill_manual(name = "Population",values=cols3,guide=FALSE)+
    geom_point(data=Pi.mean, aes(x=pop,y=mean ) )+
    theme_minimal()

ggplot() +
    geom_boxplot(data=thetas, aes(x=pop, y=pi, fill=Group), outlier.alpha = 0.3, outlier.size = 0.5, outlier.color = "gray") +
    xlab("\npopulation/year") + ylab(expression(pi))+
    scale_fill_manual(name = "Population",values=cols3,guide=FALSE)+
    geom_point(data=Pi.mean, aes(x=pop,y=mean ) )+
    theme_minimal()
ggsave("Output/Pi/Mean_pi_perPop.pdf", width = 8, height = 5)

#zoom in to pi<0.0075
ggplot() +
    geom_boxplot(data=thetas, aes(x=pop, y=pi, fill=Group), outlier.alpha = 0.5, outlier.size = 0.5, outlier.color = "gray") +
    xlab("\npopulation/year") + ylab(expression(pi))+
    scale_fill_manual(name = "Population",values=cols3,guide=FALSE)+
    geom_point(data=Pi.mean, aes(x=pop,y=mean ) )+
    theme_minimal()+ylim(0, 0.006)
ggsave("Output/Pi/Mean_pi_perPop_zommed.pdf", width = 8, height = 4)


### Theta and pi estimates from a new SFS file (MD7000, no caps in # of sites):

#Mean Pi by chromosome from Joe's downsampled theta estimates
pi<-read.csv("Output/Pi/mean_pi_byChrom.csv", row.names = 1)

#New theta estimates
theta<-read.delim('Data/new_vcf/angsd/fromBam/folded/PWS07.thetas.idx.pestPG')
theta$pi<-theta$tP/theta$nSites
theta$ch<-as.integer(gsub("chr", "", theta$Chr))
theta$ch<-factor(theta$ch, levels=1:26)
theta<-theta[order(theta$ch),]

pi$PWS07_new<-theta$pi

pi2<-pi[,c("chr", "PWS07","PWS07_new")]
pi2m<-melt(pi2, id.vars="chr")

ggplot(pi2m,  aes(x=chr, y=value, color=variable))+
    geom_point()+
     geom_path()+
    theme_bw()+
    ylab(expression(pi))+xlab("Chromosome")+theme(legend.title = element_blank())+
    scale_color_manual(values=c(blu,red),labels=c("PWS07_downsampled", "PWS07_new"))
ggsave("Output/Pi/Pi_PWS07_comparison_downsampled.vs.allsamples.pdf", width = 8, height = 4)




######  Pi comparison between SFS from Joe and SFS using all samples and all sites   ########

pops<-c("PWS91","PWS96","PWS07","PWS17")
names<-c("Downsampled","AllSamples")
#1. unfolded -for comparison & Fay's H only 
old<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/Theta/unfolded/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","pi","Tajima" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    old<-rbind(old, df1)
}

new<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromBam/unfolded/',pops[i],'_50kwin_10kstep.pestPG'))
    theta2$pi<-theta2$tP/theta2$nSites
    df2<-theta2[,c("Chr","WinCenter","pi","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    new<-rbind(new, df2)
}
theta<-rbind(old, new)    

# Add pi estiamtes from maf00 files
vcf<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromVCF/',pops[i],'_maf00.thetas50kWindow.gz.pestPG'))
    theta2$pi<-theta2$tP/theta2$nSites
    df2<-theta2[,c("Chr","WinCenter","pi","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-"vcf.maf00"
    vcf<-rbind(vcf, df2)
}

theta<-rbind(theta, vcf)    
aggregate(theta$pi, by=list(theta$pop,theta$data), mean, na.rm=T) 
#   Group.1     Group.2           x
#1    PWS07  AllSamples 0.004203786
#2    PWS17  AllSamples 0.003938429
#3    PWS91  AllSamples 0.004129560
#4    PWS96  AllSamples 0.004183521
#5    PWS07 Downsampled 0.005527874
#6    PWS17 Downsampled 0.004981327
#7    PWS91 Downsampled 0.005205000
#8    PWS96 Downsampled 0.005164741
#9    PWS07   vcf.maf00 0.025836189
#10   PWS17   vcf.maf00 0.018439158
#11   PWS91   vcf.maf00 0.019099129
#12   PWS96   vcf.maf00 0.021110088

#mean pi and Tajima's D for folded
theta$pop<-factor(theta$pop, levels=pops)

ggplot(theta, aes(x=pop, y=pi, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab(expression(pi))+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5,3.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/Pi_estimates_comparison_bam.vcf.pdf", width = 6, height = 4.5)

#from bam only
ggplot(theta[theta$data!="vcf.maf00",], aes(x=pop, y=pi, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab(expression(pi))+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5,3.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red,gre))
ggsave("Output/Pi/Pi_estimates_comparison_downsampled.vs.all.pdf", width = 6, height = 4)



## 2. folded SFS
old<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","pi","Tajima" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    old<-rbind(old, df1)
}

new<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromBam/folded/',pops[i],'_50kwin_10kstep.pestPG'))
    theta2$pi<-theta2$tP/theta2$nSites
    df2<-theta2[,c("Chr","WinCenter","pi","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    new<-rbind(new, df2)
}
    
theta<-rbind(old, new)    
    
aggregate(theta$pi, by=list(theta$pop,theta$data), mean, na.rm=T) 
#2. folded
#Group.1     Group.2           x
#1.  PWS07  AllSamples 0.004196675
#2   PWS17  AllSamples 0.003933540
#3   PWS96  AllSamples 0.004176939
#4   PWS07 Downsampled 0.003726764
#5   PWS17 Downsampled 0.003289054
#6   PWS96 Downsampled 0.003428644

# relative differences between years
sum<-data.frame(aggregate(theta$pi, by=list(theta$pop,theta$data), mean, na.rm=T) )

# pws07 vs. pws17
sum$x[1]/sum$x[2] #1.067376
# pws96/pws17
sum$x[3]/sum$x[2] #1.062231

sum$x[4]/sum$x[5] #1.109719
# pws96/pws17
sum$x[6]/sum$x[5] #1.03682


aggregate(theta$pi, by=list(theta$pop,theta$data), median, na.rm=T)    

aggregate(theta$Tajima, by=list(theta$pop,theta$data), mean, na.rm=T)    
# Group.1     Group.2          x
#1   PWS07  AllSamples -1.9175842
#2   PWS17  AllSamples -1.8654831
#3   PWS07 Downsampled -0.6132072
#4   PWS17 Downsampled -0.2510502

#folded
#1   PWS07  AllSamples -1.921720
#2   PWS17  AllSamples -1.869812
#3   PWS96  AllSamples -1.961164
#4   PWS07 Downsampled -1.396696
#5   PWS17 Downsampled -1.136874
#6   PWS96 Downsampled -1.183824



#mean pi and Tajima's D for folded
theta$pop<-factor(theta$pop, levels=pops)
ggplot(theta, aes(x=pop, y=pi, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab(expression(pi))+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red))




### Folded vs. unfolded comparison ####

######  Pi comparison between SFS from Joe   ########

pops<-c("PWS91","PWS96","PWS07","PWS17")
names<-c("Folded","Unfolded")
fold<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","pi","Tajima" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    fold<-rbind(fold, df1)
}

unfold<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/Theta/unfolded/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG'))
    theta2$pi<-theta2$tP/theta2$nSites
    df2<-theta2[,c("Chr","WinCenter","pi","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    unfold<-rbind(unfold, df2)
}

theta<-rbind(fold, unfold)    

aggregate(theta$pi, by=list(theta$pop,theta$data), mean, na.rm=T)    
aggregate(theta$pi, by=list(theta$pop,theta$data), median, na.rm=T)    

aggregate(theta$Tajima, by=list(theta$pop,theta$data), mean, na.rm=T)    


theta$pop<-factor(theta$pop, levels=pops)
ggplot(theta, aes(x=pop, y=pi, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab(expression(pi))+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5, 2.5,3.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red))
ggsave("Output/Pi/Pi_compare_downsampled.foldedvs.unfoldedPWS.pdf", width = 5.5, height = 4)



###### 2.  Pi comparison between fold vs. unfolded SFS from using all samples and all sites   ########

pops<-c("PWS96","PWS07","PWS17")
names<-c("Folded","Unfolded")
fold<-data.frame()
for (i in 1:length(pops)){
    theta1<-read.delim(paste0('Data/new_vcf/angsd/fromBam/folded/',pops[i],'_50kwin_10kstep.pestPG'))
    theta1$pi<-theta1$tP/theta1$nSites
    df1<-theta1[,c("Chr","WinCenter","pi","Tajima" )]
    df1$pop<-pops[i]
    df1$data<-names[1]
    fold<-rbind(fold, df1)
}

unfold<-data.frame()
for (i in 1:length(pops)){
    theta2<-read.delim(paste0('Data/new_vcf/angsd/fromBam/unfolded/',pops[i],'_50kwin_10kstep.pestPG'))
    theta2$pi<-theta2$tP/theta2$nSites
    df2<-theta2[,c("Chr","WinCenter","pi","Tajima" )]
    df2$pop<-pops[i]    
    df2$data<-names[2]
    unfold<-rbind(unfold, df2)
}

theta<-rbind(fold, unfold)    

aggregate(theta$pi, by=list(theta$pop,theta$data), mean, na.rm=T) 
#  Group.1  Group.2           x
#1   PWS07   Folded 0.004196675
#2   PWS17   Folded 0.003933540
#3   PWS96   Folded 0.004176939
#4   PWS07 Unfolded 0.004203786
#5   PWS17 Unfolded 0.003938429
#6   PWS96 Unfolded 0.004183521
#

aggregate(theta$Tajima, by=list(theta$pop,theta$data), mean, na.rm=T)    
# Group.1  Group.2         x
#1   PWS07   Folded -1.921720
#2   PWS17   Folded -1.869812
#3   PWS07 Unfolded -1.917584
#4   PWS17 Unfolded -1.865483

theta$pop<-factor(theta$pop, levels=pops)
ggplot(theta, aes(x=pop, y=pi, color=data))+
    geom_boxplot(position=position_dodge(width = 0.8), outlier.alpha = 0.6,outlier.size = 0.7,width=0.6)+
    geom_point(stat = "summary", fun = "mean",position=position_dodge(width = 0.8))+
    ylab(expression(pi))+xlab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = c(1.5,2.5), color="gray", size=0.5)+
    scale_color_manual(values=c(blu,red))
ggsave("Output/Pi/Pi_compare_foldvs.unfold.allSitesSamples.pws.pdf", width = 6.5, height = 4)

