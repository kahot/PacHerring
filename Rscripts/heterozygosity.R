library(ggplot2)
library(reshape2)
library(colorspace)
source("Rscripts/Pcorrection.R")
source("Rscripts/BaseScripts.R")

#Heterozygosity estimate from ANGSD (not helpful)
#######
files<-list.files("~/Projects/PacHerring/Data/new_vcf/angsd/fromVCF/estml/", pattern="_est.ml")

het<-data.frame(pop=gsub("_maf00_est.ml",'',files))

for (i in 1:length(files)){
    a<-scan(paste0("Data/new_vcf/angsd/fromVCF/estml/", files[i]))
    het$Het[i]<-a[2]/sum(a)
}

#    pop          Het
#1   BC17 0.000000e+00
#2   CA17 0.000000e+00
#3  PWS07 0.000000e+00
#4  PWS17 0.000000e+00
#5  PWS91 3.999496e-14
#6  PWS96 0.000000e+00
#7   SS06 5.582628e-13
#8   SS17 0.000000e+00
#9   SS96 0.000000e+00
#10  TB06 1.112637e-13
#11  TB17 4.213958e-14
#12  TB91 0.000000e+00
#13  TB96 3.201209e-11
#14  WA17 0.000000e+00


##########

####### 
# Calculate Ho and He from bcftools stats output files

sfiles<-list.files("Data/het/", pattern="output")
pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)

Het_sum<-data.frame(pop=pops)

for (i in 1: length(sfiles)){
    df<-read.table(paste0("Data/het/", sfiles[i]), sep="\t", header=F)
    df<-df[,c(3:10, 14:15)]
    colnames(df)<-c("Sample","nRefHom","nNonRefHom","nHets", "nTransitions", "nTransversions","nIndels","average depth","nMissing","window_no")
    df$p<-(2*df$nRefHom+df$nHets)/(rowSums(df[,c("nRefHom","nNonRefHom","nHets")])*2)
    df$q<-(2*df$nNonRefHom+df$nHets)/(rowSums(df[,c("nRefHom","nNonRefHom","nHets")])*2)
    
    df$He<-2*df$p*df$q
    df$Ho<-df$nHets/rowSums(df[,c("nRefHom","nNonRefHom","nHets")])
    Ho<-aggregate(df[,"Ho"], by=list(df$window_no), mean )
    He<-aggregate(df[,"He"], by=list(df$window_no), mean )
    het<-cbind(Ho, He$x)
    colnames(het)<-c("window_id","Ho","He")
    write.csv(het,paste0("Output/Stats_window/Hetero_",pops[i],".csv"))
    
    Het_sum$Ho_mean[i]<-  mean(het$Ho, na.rm=T)
    Het_sum$Ho_median[i]<-median(het$Ho,na.rm=T)
    Het_sum$He_mean[i]<-  mean(het$He,na.rm=T)
    Het_sum$He_median[i]<-median(het$He,na.rm=T)
    print(i)
}

write.csv(Het_sum, "Output/Stats_window/Heterozygosity_pop_summary.csv")

Het_sum<-read.csv("Output/Stats_window/Heterozygosity_pop_summary.csv", row.names = 1)


Hetm<-melt(Het_sum[,c(1,2,4)], id.vars="pop")
Hetm2<-melt(Het_sum[,c(1,3,5)], id.vars="pop")


Hetm$pop<-factor(Hetm$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))
Hetm2$pop<-factor(Hetm2$pop,levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))

ggplot(Hetm,aes(x=pop, y=value, color=variable))+
    geom_point(size=3)+
    theme_classic()+
    xlab('')+ylab("Heterozygosity")+
    theme(legend.title = element_blank())+
    ggtitle("mean")
ggsave("Output/Stats_window/Hetero_mean_allPopulations.png", width = 8, height = 5, dpi=300 )

ggplot(Hetm2,aes(x=pop, y=value, color=variable))+
    geom_point(size=3)+
    theme_classic()+theme(legend.title = element_blank())+
    xlab('')+ylab("Heterozygosity")+
    ggtitle("median")
ggsave("Output/Stats_window/Hetero_median_allPopulations.png", width = 8, height = 5, dpi=300 )


hfiles<-list.files("Output/Stats_window/", pattern="^Hetero_.+.csv")

Hetero_all<-data.frame()                  
for (i in 1:length(hfiles)){
    df<-read.csv(paste0("Output/Stats_window/", hfiles[i]), row.names = 1, stringsAsFactors = F)
    df[,1]<-pops[i]
    Hetero_all<-rbind(Hetero_all, df)
}

colnames(Hetero_all)[1]<-"pop"
Heterom<-melt(Hetero_all,id.vars="pop")
Heterom$pop<-factor(Heterom$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))
Heterom$value<-as.numeric(Heterom$value)
Hetero_all$location<-substr(Hetero_all$pop, 1,3)
Hetero_all$location<-gsub("1|0|9",'',Hetero_all$location)
Hetero_all$location<-factor(Hetero_all$location, levels=c("TB","PWS",'SS','BC','WA','CA'))
Hetero_all$pop<-factor(Hetero_all$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))


ggplot(data=Hetero_all, aes(x=pop, y=Ho, color=location, fill=location))+
    geom_boxplot(outlier.alpha = 0.2,  alpha=0.6)+
    geom_point(stat="summary", fun=mean, size=2)+
    theme_classic()+xlab('')+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=paste0(cols, "99"))+
    theme(legend.title = element_blank())
ggsave("Output/Stats_window/Ho_byPops.png",height =4, width = 8, dpi=300 )

ggplot(data=Hetero_all, aes(x=pop, y=Ho, color=location, fill=location))+
    geom_boxplot(outlier.alpha = 0.2,  alpha=0.6, outlier.size = 0.6)+
    geom_point(stat="summary", fun=mean, size=2)+
    theme_classic()+xlab('')+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=paste0(cols,"33"))+
    theme(legend.title = element_blank())+
    ylim(0,0.008)
ggsave("Output/Stats_window/Ho_byPops_ylim2.png",height =3.4, width = 8.2, dpi=300 )

    
ggplot(data=Hetero_all, aes(x=pop, y=He,, color=location, fill=location))+
    geom_boxplot(outlier.alpha = 0.2)+
    geom_point(stat="summary", fun=mean, size=2)+
    theme_classic()+xlab('')+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=paste0(cols, "99"))+
    theme(legend.title = element_blank())
ggsave("Output/Stats_window/He_byPops.png",height =4, width = 8, dpi=300)


ggplot(Heterom, aes(x=pop, y=value, color=variable, fill=variable))+
    geom_boxplot(position=position_dodge(), outlier.alpha = 0.2, outlier.size = 0.7)+
    scale_color_manual(values=c(blu, gre))+
    scale_fill_manual(values=paste0(c(blu, gre),"99"))+xlab('')+ylab("Heterozygosity")+
    theme_classic()+
    theme(legend.title = element_blank())+
    ylim(0,0.04)
ggsave("Output/Stats_window/HeHo_byPops.png",height =4, width = 10, dpi=300)


ggplot(Hetero_all, aes(x=He, color=pop))+
    geom_histogram(position=position_dodge())

ggplot(Hetero_all, aes(x=He, color=pop))+
    geom_freqpoly()

ggplot(Hetero_all, aes(x=He, color=pop))+
    geom_density()+
    theme_classic()

# Population by populations

#PWS
p<-Hetero_all[Hetero_all$location=="PWS",1:3]
p$pop<-factor(p$pop, levels=c("PWS91","PWS96","PWS07","PWS17"))
ggplot(p, aes(x=He, color=pop))+
    geom_freqpoly(bins=150, size=0.6)+
    theme_classic()

ggplot(p, aes(x=He, color=pop))+
    geom_density()+
    theme_classic()

ggplot(p, aes(x=Ho, color=pop))+
    geom_density()+
    theme_classic()+theme(legend.title = element_blank())
ggplot(p, aes(x=Ho, color=pop))+
    geom_freqpoly(bins=150, size=0.6)+
    theme_classic()

pm<-melt(p,id.vars="pop")

library(scales)
show_col(hue_pal()(4))
show_col(hue_pal()(5))
cols4<-c("#00B0F6","#00BF7D","#E76BF3", "#F8766D")

ggplot(pm, aes(x=value, color=pop, linetype=variable))+
    geom_freqpoly(bins=150, size=0.6)+
    scale_color_manual(values=cols4)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    xlab("H")+xlim(0,0.045)+
    ggtitle("PWS")
ggsave("Output/Stats_window/Het/PWS_Ho_He_distribution.png", height = 5, width = 7, dpi=300)


s<-Heterom[grep("SS", Heterom$pop),]
ggplot(s, aes(x=value, color=pop, linetype=variable))+
    geom_freqpoly(bins=150, size=0.6)+
    scale_color_manual(values=cols4[2:4])+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    xlab("H")+xlim(0,0.045)+
    ggtitle("SS")
ggsave("Output/Stats_window/Het/SS_Ho_He_distribution.png", height = 5, width = 7, dpi=300)

tb<-Heterom[grep("TB", Heterom$pop),]
ggplot(tb, aes(x=value, color=pop, linetype=variable))+
    geom_freqpoly(bins=150, size=0.6)+
    scale_color_manual(values=cols4)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    xlab("H")+xlim(0,0.045)+
    ggtitle("TB")
ggsave("Output/Stats_window/Het/TB_Ho_He_distribution.png", height = 5, width = 7, dpi=300)

##Ho only

bins=150, size=0.6


ggplot(p, aes(x=Ho, color=pop))+
    geom_freqpoly(bins=50, size=0.6)+
    scale_color_manual(values=cols4)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    ggtitle("PWS")+xlim(0,0.02)
ggsave("Output/Stats_window/Het/PWS_Ho_distribution.png", height = 5, width = 7, dpi=300)

ss<-Hetero_all[Hetero_all$loc=="SS",]
ggplot(ss, aes(x=Ho, color=pop))+
    geom_freqpoly(bins=50, size=0.6)+
    scale_color_manual(values=cols4[2:4])+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    ggtitle("SS")+xlim(0,0.02)
ggsave("Output/Stats_window/Het/SS_Ho_distribution.png", height = 5, width = 7, dpi=300)

t<-Hetero_all[Hetero_all$loc=="TB",]
ggplot(t, aes(x=Ho, color=pop))+
    geom_freqpoly(bins=50, size=0.6)+
    scale_color_manual(values=cols4)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    ggtitle("TB")+xlim(0,0.02)
ggsave("Output/Stats_window/Het/TB_Ho_distribution.png", height = 5, width = 7, dpi=300)


y17<-Hetero_all[grep("17",Hetero_all$pop), ]
ggplot(y17, aes(x=Ho, color=pop))+
    geom_freqpoly()+
    #scale_color_manual(values=cols)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    ggtitle("Y2017")+xlim(0,0.02)
ggsave("Output/Stats_window/Het/Year2017_Ho_distribution.png", height = 5, width = 7, dpi=300)

ggplot(y17, aes(x=He, color=pop))+
    geom_freqpoly()+
    #scale_color_manual(values=cols)+
    theme_classic()+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 3) ) )+
    ggtitle("Y2017")+xlim(0,0.04)
ggsave("Output/Stats_window/Het/Year2017_He_distribution.png", height = 5, width = 7, dpi=300)




#run statistical test

#PWS 91 vs. 96 etc.
pops<-unique(Het_sum$pop)
comb<-combn(pops,2)
comb<-t(comb)
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

wil.res<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    v1<-Hetero_all$Ho[Hetero_all$pop==comb[i,1]]
    v2<-Hetero_all$Ho[Hetero_all$pop==comb[i,2]]
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res$rawP[i]<-re1[[3]]
    wil.res$mean1[i]<-mean(v1, na.rm=T)
    wil.res$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
}

wil.res<-Pcorrection(wil.res)
write.csv(wil.res, paste0("Output/Stats_window/Het/Wilcox_results_Ho_all.csv"))


pwsRes<-wil.res[grep("PWS",wil.res$Test),]
# PWS07 is higher than all other years. 
#PWS07 > PWS96 > PWS91 > PWS17
# Losing heterozygosity in 2017


ssRes<-wil.res[grep("SS",wil.res$Test),]
#SS96 > SS17 > SS06  

tbRes<-wil.res[grep("TB",wil.res$Test),]
# TB91 > TB06 > TB96 > TB17


wil.res2<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    v1<-Hetero_all$He[Hetero_all$pop==comb[i,1]]
    v2<-Hetero_all$He[Hetero_all$pop==comb[i,2]]
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res2$rawP[i]<-re1[[3]]
    wil.res2$mean1[i]<-mean(v1, na.rm=T)
    wil.res2$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res2$higher[i]<-ifelse((wil.res2$mean1[i]-wil.res2$mean2[i])>0, comb[i,1],comb[i,2])
}

wil.res2<-Pcorrection(wil.res2)
write.csv(wil.res2, paste0("Output/Stats_window/Het/Wilcox_results_He_all.csv"))







### Nuclotide diversity
pifiles<-list.files("Output/Pi/", pattern="Pi_pixy_per50kWindow.csv")
Pi<-data.frame()
for (i in 1: length(pifiles)){
    pop.name<-gsub("_Pi_pixy_per50kWindow.csv", "", pifiles[i])
    df<-read.csv(paste0("Output/Pi/", pifiles[i]), row.names = 1)
    df$pop<-pop.name
    Pi<-rbind(Pi, df[,c("pop","Chr","WinCenter","pi")])
}

Pi$pop<-factor(Pi$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))
Pi$Chr<-factor(Pi$Chr, levels=paste0("chr", 1:26))
Pi$location<-substr(Pi$pop,1,2)
avePi<-aggregate(Pi[,c("pi")], by=list(Pi$pop), mean)
colnames(avePi)<-c("pop","mean")
avePi$pop<-factor(avePi$pop, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))

ggplot()+
    geom_boxplot(data=Pi, aes(x=pop, y=pi, color=location, fill=location), position=position_dodge(), outlier.alpha = 0.2, alpha=0.6)+
    geom_point(data=avePi, aes(x=pop, y=mean))+
    xlab('')+ylab("Pi (100kb window)")+
    theme_classic()+
    theme(legend.title = element_blank())
#ggsave("Output/Stats_window/pi/Pi_byPops.pdf",height =4, width = 8 )


ggplot(Pi, aes(x=pop, y=pi, color=Chr, fill=Chr))+
    geom_bar(stat="identity")+
    xlab('')+ylab("Pi (100kb window)")+
    theme_classic()+
    theme(legend.title = element_blank())+
    guides(fill = guide_legend(override.aes = list(size = 2) ) )+
    theme(axis.text.x = element_text(angle=45, hjust=1))
#ggsave("Output/Stats_window/pi/Pi_byPops_byChr.pdf",height =4, width = 8 )


wil.res3<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    v1<-Pi$pi[Pi$pop==comb[i,1]]
    v2<-Pi$pi[Pi$pop==comb[i,2]]
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res3$rawP[i]<-re1[[3]]
    wil.res3$mean1[i]<-mean(v1, na.rm=T)
    wil.res3$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res3$higher[i]<-ifelse((wil.res3$mean1[i]-wil.res3$mean2[i])>0, comb[i,1],comb[i,2])
}

wil.res3<-Pcorrection(wil.res3)
write.csv(wil.res3, paste0("Output/Pi/Wilcox_results_pi_all.csv"))


#### Theta from Joe

Theta<-data.frame(pop=pops)
AllPi<-data.frame()
for (i in 1:length(pops)){
    
    thetas1 <-read.delim(paste('Data/Theta/',pops[i],'_minQ20_minMQ30_50kb_win_10kb_step.pestPG',sep = ""))
    #thetas1$pi<- (thetas1$tP*thetas1$nSites)/10000
    thetas1$pi<- thetas1$tP/thetas1$nSites
    thetas1<-thetas1[,c("Chr","pi")]
    thetas1$pop<-pops[i]
    AllPi<-rbind(AllPi, thetas1)
    Theta$mean[i]<-mean(thetas1$pi, na.rm=T)
    Theta$median[i]<-median(thetas1$pi, na.rm=T)
}


ggplot()+
    geom_point(data=AllPi, aes(x=pop, y=pi), color="lightblue", position=position_jitter(width = .1), size=0.3)+
    geom_point(data=Theta, aes(x=pop, y=mean))+
    theme_bw()+ylab("Pi")+xlab('')+
    theme(panel.grid.major.x = element_blank())

AllPi$loc<-substr(AllPi$pop, 1,2)
Theta$loc<-substr(Theta$pop, 1,2)

AllPi$pop<-factor(AllPi$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))
Theta$pop<-factor(Theta$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))

ggplot()+
    geom_boxplot(data=AllPi, aes(x=pop, y=pi, color=loc, fill=loc), outlier.alpha = 0.2, alpha=0.6)+
    geom_point(data=Theta, aes(x=pop, y=mean))+
    theme_classic()+ylab("Pi")+xlab('')+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Stats_window/pi/Pi_from_Theta_joe.pdf", width = 8, height = 5 )


