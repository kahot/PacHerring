
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(vcfR)

bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'

# allele frequency changes over time
pops<-c("PW","SS","TB")

files<-list.files("Data/freqs/persite/")

window_size <- 100000
step_size <- 10000

#look at the mean etc.
summary=data.frame(pop=pops)

for (i in 1:length(files)){
    print(files[i])
    df<-read.table(paste0("Data/freqs/persite/", files[i]))
    if (i==2) summary$t0_mean[i]<-NA
    else summary$t0_mean[i]<-mean(df$t0_AF,na.rm=T)
    summary$t1_mean[i]<-mean(df$t1_AF,na.rm=T)
    summary$t2_mean[i]<-mean(df$t2_AF,na.rm=T)
    summary$t3_mean[i]<-mean(df$t3_AF,na.rm=T)
    if (i==2) summary$t0_median[i]<-NA
    else summary$t0_median[i]<-median(df$t0_AF,na.rm=T)
    summary$t1_median[i]<-median(df$t1_AF,na.rm=T)
    summary$t2_median[i]<-median(df$t2_AF,na.rm=T)
    summary$t3_median[i]<-median(df$t3_AF,na.rm=T)
}

write.csv(summary,"Output/AF/AF_change_overtime_mean_median.csv")

summ<-melt(summary,id.vars="pop")
summ$yr<-substr(summ$variable, 1,2)    
summ$stat<-substr(summ$variable, 4,9) 

ggplot(data=summ, aes(x=yr, y=value, fill=stat))+
    geom_bar(stat="identity",position=position_dodge())+
    facet_wrap(~pop)+
    scale_fill_manual(values=c(grb,blu))+
    theme_bw()+
    ylab("None-ref AF")+
    xlab('')+
    theme(legend.title = element_blank())
ggsave("Output/AF/MeanAF_overTime.pdf", width = 8, height=4.5)

AFchange<-data.frame(pop=pops)
for (i in 1:length(files)){
    print(files[i])
    df<-read.table(paste0("Data/freqs/persite/", files[i]))
    length(df$pos)
    if (i==2) AFchange$nLoci_AltMaj_t0[i]<-NA
    else AFchange$nLoci_AltMaj_t0[i]<-nrow(df[df$t0_AF>0.5,])
    AFchange$nLoci_AltMaj_t1[i]<-nrow(df[df$t1_AF>0.5,])
    AFchange$nLoci_AltMaj_t2[i]<-nrow(df[df$t2_AF>0.5,])
    AFchange$nLoci_AltMaj_t3[i]<-nrow(df[df$t3_AF>0.5,])
    
    if (i==2) AFchange$percentAltMaj_t0[i]<-NA
    else AFchange$percentAltMaj_t0[i]<-AFchange$nLoci_AltMaj_t0[i]/length(df$pos)*100
    AFchange$percentAltMaj_t1[i]<-AFchange$nLoci_AltMaj_t1[i]/length(df$pos)*100
    AFchange$percentAltMaj_t2[i]<-AFchange$nLoci_AltMaj_t2[i]/length(df$pos)*100
    AFchange$percentAltMaj_t3[i]<-AFchange$nLoci_AltMaj_t3[i]/length(df$pos)*100
} 
write.csv(AFchange,"Output/AF/nonRefAllele_percent_summary.csv")
    
afPer<-melt(AFchange[,c(1,6:9)], id.vars="pop")

afPer$yr<-gsub("percentAltMaj_",'', afPer$variable)    


ggplot(afPer, aes(x=yr, y=value))+facet_wrap(~pop)+
    geom_bar(stat="identity", position=position_dodge(), fill=grb, width =0.7)+
    ylab("% nonRef allele >50%")+
    
    xlab('')+
    theme(legend.title = element_blank())
ggsave("Output/AF/nonRefAllele_percent_overTime.pdf", width = 6, height = 4.5)
    

#which sites? are they same over years? 
#for (i in 1:length(files)){
    print(files[i])
    df<-read.table(paste0("Data/freqs/persite/", files[i]))
    
    if (i==2) AFchange$nLoci_AltMaj_t0[i]<-NA
    #else t0<-df[df$t0_AF>0.5,]
    t1<-df[df$t1_AF>0.5,]
    t2<-df[df$t2_AF>0.5,]
    t3<-df[df$t3_AF>0.5,]

    t0.t1<-intersect(t0$pos,t1$pos)#1366
    t1.t2<-intersect(t1$pos,t2$pos)#1448
    t2.t3<-intersect(t2$pos,t3$pos) #1353
    t0.t3<-intersect(t0$pos,t3$pos) #1279
    
    t012<-intersect(t0.t1, t2$pos) #1211
    t0123<-intersect(t0.t1,t2.t3) #1071
    
    Overlap<-df[df$pos %in% t0123,1:6]
    
    Overlap$delta<-apply(Overlap[,3:6], 1, function(x) max(x[3:6], na.rm=T)-min(x[3:6], na.rm=T))
    
    plot(Overlap$delta)
    
    overlap<-df[df$pos %in% t0123,2:6]
    overt<-t(overlap)
    colnames(overt)<-paste0("pos.",overt[1,])
    overt<-overt[-1,]
    overt<-data.frame(overt)
    overt$Time<-as.integer(substr(rownames(overt),2,2))
    overm<-melt(overt, id.vars="Time")
    
    poss<-unique(overm$variable)
    ov1<-overm[overm$variable %in% poss[1:20],]
    ggplot(ov1, aes(x=Time, y=value, color=variable))+
        geom_line()+ylim(0.49,1)+
        theme(legend.position = "none") 
    

    ggplot(overm, aes(x=Time, y=value, color=variable))+
        geom_line()+ylim(0.49,1)+
        theme(legend.position = "none") 
    
    
### find the sites that are fixed to the alternate and use that as a reference:
    
freqs<-list.files("Output/AF/", pattern="_freq.frq")
Alt<-list()
for (i in 1: length (freqs)){
    frq<-read.table(paste0("Output/AF/",freqs[i]), sep="\t", row.names=NULL)
    colnames(frq)[1:4]<-colnames(frq)[2:5]
    colnames(frq)[5:6]<-c("f1","f2")
    frq$pos<-paste0(frq$CHROM,":",frq$POS)
    
    frq$nuc1<-substr(frq$f1,1,1)
    frq$nuc2<-substr(frq$f2,1,1)
    frq$freq1<-as.numeric(substr(frq$f1,3,10))
    frq$freq2<-as.numeric(substr(frq$f2,3,10))
    write.csv(frq, paste0("Output/AF/AF_", gsub("_freq.frq","",freqs[i]),".csv"))
    
    df<-frq[frq$freq2>0.5,c(1:2,7:11)]
    Alt[[i]]<-df
    names(Alt)[i]<- gsub("_freq.frq","",freqs[i])
}


#Find fixed to alternate positions
#exclude CA and TB
common_loci <- Reduce(intersect, list(Alt[[1]]$pos,Alt[[3]]$pos,Alt[[4]]$pos,Alt[[5]]$pos,Alt[[6]]$pos,Alt[[7]]$pos,Alt[[8]]$pos,Alt[[9]]$pos,Alt[[14]]$pos))
#874 loci

#Add TB
common_tb<-Reduce(intersect, list(Alt[[10]]$pos, Alt[[11]]$pos,Alt[[12]]$pos,Alt[[13]]$pos))
#2058
common_loci2<-Reduce(intersect, list(common_loci,common_tb))
#401 loci only

#Add CA
common_loci3<-Reduce(intersect, list(common_loci,Alt[[2]]$pos))  #824

# the 874 sites should be reversed when calculating allele freq changes


df<-Alt[[1]]
df<-df[df$pos %in% common_loci,]
names(Alt)


#PWS
pw<-read.table(paste0("Data/freqs/persite/PWS_shifts_persite.txt"))
pw$loci<-paste0(pw$chr,":",pw$pos)

#The divide the files into 2

pw1<-pw[!(pw$loci%in% common_loci),]
pw2<-pw[(pw$loci%in% common_loci),]
pw2master<-pw2

pw2$t0_AF<-1-pw2$t0_AF
pw2$t1_AF<-1-pw2$t1_AF
pw2$t2_AF<-1-pw2$t2_AF
pw2$t3_AF<-1-pw2$t3_AF
pw2$t0_transformed_freq<-asin(sqrt(pw2$t0_AF))
pw2$t1_transformed_freq<-asin(sqrt(pw2$t1_AF))
pw2$t2_transformed_freq<-asin(sqrt(pw2$t2_AF))
pw2$t3_transformed_freq<-asin(sqrt(pw2$t3_AF))
pw2$zt01<-pw2$t1_transformed_freq-pw2$t0_transformed_freq
pw2$zt02<-pw2$t2_transformed_freq-pw2$t0_transformed_freq
pw2$zt03<-pw2$t3_transformed_freq-pw2$t0_transformed_freq
pw2$zt12<-pw2$t2_transformed_freq-pw2$t1_transformed_freq
pw2$zt13<-pw2$t3_transformed_freq-pw2$t1_transformed_freq
pw2$zt23<-pw2$t3_transformed_freq-pw2$t2_transformed_freq


pw_new<-rbind(pw1,pw2)
colMeans(pw_new[,11:16])
#         zt01          zt02          zt03          zt12          zt13          zt23 
#0.0088443887  0.0204686041 -0.0008572306  0.0116242153 -0.0097016193 -0.0213258346 
colMeans(pw_new[,11:16])
#zt01         zt02         zt03         zt12         zt13         zt23 
#0.008768065  0.020383915 -0.000792345  0.011615851 -0.009560410 -0.021176260 

#Did not differ much

write.csv(pw_new, "Output/AF/PWS_shifts_persite_new.txt")
 


## Look at the AF near CYP1A
# CYP1A is at chr6 2,922,179 - 2,927,184
cyp<-pw_new[pw_new$chr=="chr6",]
cyp1<-cyp[cyp$pos>= 27500000 & cyp$pos <=3250000,]

cyp1<-cyp1[,2:6]
colnames(cyp1)[2:5]<-c(1991,1996,2007,2017)
cyp1m<-melt(cyp1, id.vars="pos")

ggplot(cyp1m, aes(x=variable, y=value, color=factor(pos)))+
    geom_point()+geom_line(aes(x=variable, y=value,group=factor(pos)))+
    ylim(0,0.5)+xlab('')+ylab('AF')


cyp2<-cyp[cyp$pos>= 2900000 & cyp$pos <=3000000,]

cyp2<-cyp2[,2:6]
colnames(cyp2)[2:5]<-c(1991,1996,2007,2017)
cyp2m<-melt(cyp2, id.vars="pos")

ggplot(cyp2m, aes(x=variable, y=value, color=factor(pos)))+
    geom_point()+geom_line(aes(x=variable, y=value,group=factor(pos)))+
    ylim(0,0.5)+xlab('')+ylab('AF')

ch6<-common_loci[grep("chr6",common_loci)] 
ch6


#Look at if copy number variation exist near CYP1A like the killifish
vcf <- read.vcfR('Data/vcfs/ph_chr6.vcf.gz')
vcf<-read.vcfR("Data/new_vfc/PH_MD7000_maf01_chr6.vcf.gz")
#extract the read depth information 
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)

# 1. look at the PWS populations
dp_pws<-dp[,grep("PWS",colnames(dp))]

write.csv(dp_pws, "Output/AF/PWS_chr6_read_depth.csv")

dp_ca<-dp[,grep("CA",colnames(dp))]

#extract the position as numbers from rownames
positions<-as.integer(gsub("chr6_", '', rownames(dp_pws)))
n<-which(positions>=2000000&positions<=3500000)

#near CYP1A
dp_pws2<-data.frame(dp_pws[n,])
dp_pws2$pos<-positions[n]

dp_pwsm<-melt(dp_pws2, id.vars="pos")
ggplot(dp_pwsm, aes(x=pos, y=value))+
    geom_point()

#Change the sample ids as character
dp_pwsm$variable<-as.character(dp_pwsm$variable)

#select only PWS91
p91<-dp_pwsm[grep("PWS91", dp_pwsm$variable),]
samples<-unique(p91$variable) #individual IDs

#plot each individual
plots<-list()
for (i in 1: length(samples)){
    plots[[i]]<-ggplot(p91[p91$variable==samples[i],], aes(x=pos, y=value))+
        geom_point(size=0.5)+ylim(0,max(p91$value, na.rm=T))+
        ggtitle(gsub("X","",samples[i]))
}
#pdf("Output/AF/dp7000_PWS91_copyNo_nearCYP1A.pdf", width = 15, height = 40)
pdf("Output/AF/chr6/dp7000_PWS91_ch6.20-35m.pdf", width = 15, height = 40)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()


#select PWS96
p96<-dp_pwsm[grep("PWS96", dp_pwsm$variable),]
samples<-unique(p96$variable) #individual IDs

#plot each individual
plots<-list()
for (i in 1: 30){
    plots[[i]]<-ggplot(p96[p96$variable==samples[i],], aes(x=pos, y=value))+
        geom_point(size=0.5)+ylim(0,max(p96$value, na.rm=T))+
        ggtitle(gsub("X","",samples[i]))
}
pdf("Output/AF/dp7000_PWS96_copyNo_nearCYP1A.pdf", width = 15, height = 26)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()

#select PWS07
p07<-dp_pwsm[grep("PWS07", dp_pwsm$variable),]
samples<-unique(p07$variable) #individual IDs

#plot each individual
plots<-list()
for (i in 1: 20){
    plots[[i]]<-ggplot(p07[p07$variable==samples[i],], aes(x=pos, y=value))+
        geom_point(size=0.5)+ylim(0,20)+
        ggtitle(gsub("X","",samples[i]))
}
pdf("Output/AF/dp7000_PWS07_copyNo_nearCYP1A_1.pdf", width = 15, height = 15.5)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()



#select PWS17
p17<-dp_pwsm[grep("PWS17", dp_pwsm$variable),]
samples<-unique(p17$variable) #individual IDs

#plot each individual
plots<-list()
for (i in 1: 30){
    plots[[i]]<-ggplot(p17[p17$variable==samples[i],], aes(x=pos, y=value))+
        geom_point(size=0.5)+ylim(0,20)+
        ggtitle(gsub("X","",samples[i]))
}
pdf("Output/AF/dp7000_PWS17_copyNo_nearCYP1A_1.pdf", width = 15, height = 26)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()


## CA pop
#near CYP1A
dp_ca2<-data.frame(dp_ca[n,])
dp_ca2$pos<-positions[n]

dp_cam<-melt(dp_ca2, id.vars="pos")

#Change the sample ids as character
dp_cam$variable<-as.character(dp_cam$variable)

samples<-unique(dp_cam$variable) #individual IDs
#plot each individual
plots<-list()
for (i in 1: 30){
    plots[[i]]<-ggplot(dp_cam[dp_cam$variable==samples[i],], aes(x=pos, y=value))+
        geom_point(size=0.5)+ylim(0,max(dp_cam$value, na.rm = T))+
        ggtitle(gsub("X","",samples[i]))
}
pdf("Output/AF/dp7000_CA17_copyNo_nearCYP1A_2.pdf", width = 15, height = 26)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()




##################
## Look at the maf00 file
pwspop<-c("PWS91","PWS96","PWS07","PWS17")
for (p in 1:4){
    vcf <- read.vcfR(paste0("Output/AF/Chr6/", pwspop[p],"_chr6"))
    
    #extract the read depth information 
    dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
    
    positions<-as.integer(gsub("chr6_", '', rownames(dp))) #3966
    
    dp<-data.frame(dp)
    dp$pos<-positions
    
    dpm<-melt(dp, id.vars="pos")
    samples<-unique(dpm$variable) #individual IDs
    
    #plot each individual
    plots<-list()
    for (i in 1: length(samples)){
        plots[[i]]<-ggplot(dpm[dpm$variable==samples[i],], aes(x=pos, y=value))+
            geom_point(size=0.5)+ylim(0,20)+
            ggtitle(gsub("X","",samples[i]))
    }
    pdf(paste0("Output/AF/chr6/",pwspop[p],"_copyNo_nearCYP1A.pdf"), width = 15, height = length(samples)*0.7)
    do.call(grid.arrange, c(plots, ncol=3))
    dev.off()
}
# No copy number variations observed

### check the allele freq change


pwspop<-c("PWS91","PWS96","PWS07","PWS17")
year<-c(1991,1996,2007,2017)
AF<-data.frame()
for (i in 1:4){
    vcf <- read.vcfR(paste0("Output/AF/Chr6/", pwspop[i],"_chr6"))
    
    #extract the allele freq
    af<-data.frame(maf(vcf, element=2))
     
    positions<-as.integer(gsub("chr6_", '', rownames(af))) #3966
    
    af$pos<-positions
    af$year<-year[i]
    
    AF<-rbind(AF, af[,c("pos","Frequency","year")])
}

subAF<-AF[AF$pos>2920000 & AF$pos<2930000,]
posi<-unique(subAF$pos)

vec<-seq(1,125,10)
plots<-list()
for (i in 1:length(vec)){
    df<-subAF[subAF$pos %in% posi[vec[i]:(vec[i]+9)],]
    plots[[i]]<-ggplot(df, aes(x=year, y=Frequency, color=factor(pos), group=factor(pos)))+
        geom_point()+
        geom_line ()+
        ylim(0,(max(df$Frequency)+0.01))+
        theme_classic()
    
}
pdf(paste0("Output/AF/chr6/AF_change_maf00.pdf"), width = 15, height = 18)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()

plots[[1]]


ggplot(df, aes(x=year, y=Frequency, color=factor(pos), group=factor(pos)))+
    geom_point(position=position_jitter(width=0.1,height=0.1,seed = 123))+
    geom_path(position=position_jitter(width=0.1,height=0.1, seed = 123))+
    ylim(0,(max(df$Frequency)+0.01))+
    theme_classic()



df[df$pos==2922057,]
    #plot each individual
    plots<-list()
    for (i in 1: length(samples)){
        plots[[i]]<-ggplot(dpm[dpm$variable==samples[i],], aes(x=pos, y=value))+
            geom_point(size=0.5)+ylim(0,20)+
            ggtitle(gsub("X","",samples[i]))
    }
    pdf(paste0("Output/AF/chr6/",pwspop[p],"_copyNo_nearCYP1A.pdf"), width = 15, height = length(samples)*0.7)
    do.call(grid.arrange, c(plots, ncol=3))
    dev.off()
}
























    
dp_bc<-dp[,grep("BC",colnames(dp))]
dp_bc<-data.frame(dp_bc)  
po<-
    positions<-gsub("chr6_", '', rownames(dp_pws))
positions<-as.integer(positions)
n<-which(positions>=2500000&positions<=3500000)

dp_bc$pos<-positions
bcm<-melt(dp_bc, id.vars="pos")
samples<-unique(bcm$variable)
i=1
ggplot(bcm[bcm$variable==samples[i],], aes(x=pos, y=value))+
    geom_point()+ylim(0,20)

dp_bc$mean<-rowMeans(dp_bc[,1:64], na.rm=T)

ggplot(dp_bc[],aes(x=pos, y=mean))+
    geom_point(size=0.1)

p17<-dp_pwsm[grep("PWS17", dp_pwsm$variable),]
ggplot(p17, aes(x=pos, y=value))+
    geom_point()+ylim(0,20)


par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")
