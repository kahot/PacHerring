#Selection results from PCAngsd
source("Rscripts/BaseScripts.R")
library(stringr)
library(gridExtra)
library(RcppCNPy) # Numpy library for R


cols<-c("#0072b2","#cc79a7","#009e73","#d55e00","#56b4e9","#e69f00","#f0e442")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected")]
colnames(pop_info)[3]<-"year"


#######
### Selection ###
library(RcppCNPy) # Numpy library for R

## function for QQplot
qqchi<-function(x,...){
    lambda<-round(median(x)/qchisq(0.5,1),2)
    qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
    legend("topleft",paste("lambda=",lambda))
}

### read in seleciton statistics (chi2 distributed)
# Each column reflect the selection statistics along a tested PC (they are χ²-distributed with 1 degree of freedom.)
s<-npyLoad("Data/PCAangsd/selection/Y2017_selection.selection.npy")

## make QQ plot to QC the test statistics
qqchi(s)

ncol(s)

## read positions 
p<-read.table("Data/PCAangsd/selection/Y2017_selection.sites",colC=c("factor","integer"),sep=":")
names(p)<-c("chr","pos")


# 1. more than 1 axis:
#convert test statistic to p-value
pval<-1-pchisq(s,1)
## make manhatten plot
plot(-log10(pval),col=p$chr,xlab="Chromosomes",main="Manhattan plot", pch=".")

# 2. if more than 1 axis (ncol(s)>1)
p$pval1<-pval[,1]
p$pval2<-pval[,2]
p$loc<-1:nrow(p)
p$pval1.log<--log10(p$pval1)
p$pval2.log<--log10(p$pval2)

## make manhatten plot
pdf("Output/PCA/selection/Year2017pruned_selection.plot.pdf", width = 8, height = 4)
plot(-log10(p$pval1),col=p$chr,xlab="Chromosomes",main="Year 2017", pch=".")
dev.off()


evens<-paste0("chr",seq(2,26, by=2))

#count the number of sites per chromosomes
poss<-data.frame(chr=paste0("chr",1:26))
k=1
for (i in 1:26){
    df<-p[p$chr==paste0("chr",i),]
    poss$start[i]<-k
    poss$end[i]<-k+nrow(df)-1
    k=k+nrow(df)
}

poss$x<-poss$start+(poss$end-poss$start)/2

p$color<-"steelblue"
p$color[p$chr %in% evens]<-"lightblue"
ggplot(data=p, aes(x=loc, y=pval1.log, color=color))+
    geom_point(size=0.1)+
    scale_color_manual(values=c("lightblue","steelblue"), guide='none')+
    scale_x_continuous(name="Chromosome position", breaks=poss$x, labels=1:26)+
    theme_classic()+ylab("-log10(p-value)")+
    ggtitle("PC 1")
ggsave("Output/PCA/selection/pcangsd_selection_2017pops.rpuned_pc1.png", width = 12, height = 6, dpi=300)    

ggplot(data=p, aes(x=loc, y=pval2.log, color=color))+
    geom_point(size=0.1)+
    scale_color_manual(values=c("lightblue","steelblue"), guide='none')+
    scale_x_continuous(name="Chromosome position", breaks=poss$x, labels=1:26)+
    theme_classic()+ylab("-log10(p-value)")+
    ggtitle("PC 2")
ggsave("Output/PCA/selection/pcangsd_selection_2017pops.pruned_pc2.png", width = 12, height = 6, dpi=300)    


## PWS together
s<-npyLoad("Data/PCAangsd/selection/PWS_selection.selection.npy")

## make QQ plot to QC the test statistics
qqchi(s)

# convert test statistic to p-value
pval<-1-pchisq(s,1)

## read positions 
p<-read.table("Data/PCAangsd/selection/PWS_selection.sites",colC=c("factor","integer"),sep=":")
names(p)<-c("chr","pos")

## make manhatten plot
pdf("Output/PCA/selection/PWS_selection_scan.pdf", width = 10, height = 5)
plot(-log10(pval),col=p$chr,xlab="Chromosomes",main="Manhattan plot", pch=".")
dev.off()

p$pval1<-pval[,1]
p$loc<-1:nrow(p)
p$pval1.log<--log10(p$pval1)

evens<-paste0("chr",seq(2,26, by=2))

#count the number of sites per chromosomes

poss<-data.frame(chr=paste0("chr",1:26))
k=1
for (i in 1:26){
    df<-p[p$chr==paste0("chr",i),]
    poss$start[i]<-k
    poss$end[i]<-k+nrow(df)-1
    k=k+nrow(df)
}

poss$x<-poss$start+(poss$end-poss$start)/2

p$color<-"steelblue"
p$color[p$chr %in% evens]<-"lightblue"
ggplot(data=p, aes(x=loc, y=pval1.log, color=color))+
    geom_point(size=0.1)+
    scale_color_manual(values=c("lightblue","steelblue"), guide='none')+
    scale_x_continuous(name="Chromosome position", breaks=poss$x, labels=1:26)+
    theme_classic()+ylab("-log10(p-value)")+
    ggtitle("PC 1")
ggsave("Output/PCA/selection/PWS_pcangsd_selection_pc1.png", width = 12, height = 6, dpi=150)    

ggplot(data=p, aes(x=loc, y=pval2.log, color=color))+
    geom_point(size=0.1)+
    scale_color_manual(values=c("lightblue","steelblue"), guide='none')+
    scale_x_continuous(name="Chromosome position", breaks=poss$x, labels=1:26)+
    theme_classic()+ylab("-log10(p-value)")+
    ggtitle("PC 2")
ggsave("Output/PCA/pcangsd_selection_2017pops_pc2.pdf", width = 12, height = 6)    


########  PWS pairwise run ###
pwapops<-c("PWS91" ,"PWS96","PWS07", "PWS17")
comb<-combn(pwapops, 2)
comb<-t(comb)
#comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

##TB
tbpops<-c("TB91" ,"TB96","TB06", "TB17")
comb<-combn(tbpops, 2)
comb<-t(comb)

comb1<-comb[c(1,4,6,3),]
evens<-paste0("chr",seq(2,26, by=2))

for (i in 1:nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    #s<-npyLoad(paste0("Data/PCAangsd/selection/",pop1,".",pop2,"_selection.selection.npy"))
    s<-npyLoad(paste0("Data/PCAangsd/selection/pruned_",pop1,".",pop2,"_selection.selection.npy"))
    
    #how many PC axes were evaluated?
    print(ncol(s))
    
    ## make QQ plot to QC the test statistics
    #qqchi(s)
    
    # convert test statistic to p-value
    pval<-1-pchisq(s,1)
    
    #p<-read.table(paste0("Data/PCAangsd/selection/",pop1,".",pop2,"_selection.sites"),colC=c("factor","integer"),sep=":")
    p<-read.table(paste0("Data/PCAangsd/selection/pruned_",pop1,".",pop2,"_selection.sites"),colC=c("factor","integer"),sep=":")
    names(p)<-c("chr","pos")
    
    ## make manhatten plot
    #pdf(paste0("Output/PCA/",pop,".pcangsd_selection_manhattan.plot.pdf"), width = 12, height = 6)
    #plot(-log10(pval),col=p$chr,xlab="Chromosomes",main="Manhattan plot", pch=".")
    #dev.off()
    
    p$pval<-pval
    p$loc<-1:nrow(p)
    p$pval.log<--log10(p$pval)
    
    #write.csv(p[p$pval.log>4,], paste0("Output/PCA/selection/Selected_highP_sites_", pop1,".",pop2,".csv"))
    write.csv(p[p$pval.log>4,], paste0("Output/PCA/selection/Selected_highP_sites_Pruned", pop1,".",pop2,".csv"))
    
    #count the number of sites per chromosomes
    poss<-data.frame(chr=paste0("chr",1:26))
    k=1
    for (i in 1:26){
        df<-p[p$chr==paste0("chr",i),]
        poss$start[i]<-k
        poss$end[i]<-k+nrow(df)-1
        k=k+nrow(df)
    }
    
    poss$x<-poss$start+(poss$end-poss$start)/2
    
    p$color<-"1"
    p$color[p$chr %in% evens]<-"2"
    ggplot(data=p, aes(x=loc, y=pval.log, color=color))+
        geom_point(size=0.15)+
        scale_color_manual(values=c("steelblue","lightblue"), guide='none')+
        scale_x_continuous(name="Chromosome", breaks=poss$x, labels=c(1:23,'',25,26))+
        theme_classic()+ylab("-log10(p-value)")+
        ggtitle(paste0(pop1," vs. ",pop2))
    #ggsave(paste0("Output/PCA/selection/pcangsd_selection_",pop1,"_",pop2,"_plot.png"), width = 6, height = 3, dpi = 300)    
    #ggsave(paste0("Output/PCA/selection/pcangsd_selection_",pop1,"_",pop2,"_plot.png"), width = 6, height = 3, dpi = 300)    
    
}


#compare pruned vs. non-pruned results for pws91 vs. pws17
pop1="PWS91"
pop2="PWS17"
s<-npyLoad(paste0("Data/PCAangsd/selection/",pop1,".",pop2,"_selection.selection.npy"))
sp<-npyLoad(paste0("Data/PCAangsd/selection/pruned_",pop1,".",pop2,"_selection.selection.npy"))

pval<-1-pchisq(s,1)
pval.p<-1-pchisq(sp,1)

p<-read.table(paste0("Data/PCAangsd/selection/",pop1,".",pop2,"_selection.sites"),colC=c("factor","integer"),sep=":")
pp<-read.table(paste0("Data/PCAangsd/selection/pruned_",pop1,".",pop2,"_selection.sites"),colC=c("factor","integer"),sep=":")
names(p)<-c("chr","pos")
names(pp)<-c("chr","pos")

p$pval<-pval
p$loc<-1:nrow(p)
p$pval.log<--log10(p$pval)

pp$pval<-pval.p
pp$loc<-1:nrow(pp)
pp$pval.log<--log10(pp$pval)


# how many sites over -log10(P)>4
nrow(p[p$pval.log>=4.5,]) #11 sites
nrow(p[p$pval.log>=4,]) #29 sites
nrow(pp[pp$pval.log>=4.5,]) #11
nrow(pp[pp$pval.log>=4,]) #30 sites


highPsites1<-p[p$pval.log>=4,]
highPsites2<-pp[pp$pval.log>=4,]
highPsites1$id<-paste0(highPsites1$chr,"_",highPsites1$pos)
highPsites2$id<-paste0(highPsites2$chr,"_",highPsites2$pos)

intersect(highPsites1$id, highPsites2$id)
# [1] "chr1_18838039"  "chr1_28002097"  "chr2_31558154"  "chr3_2074360"   "chr3_11717670"  "chr3_19942012" 
#[7] "chr3_22120097"  "chr6_30766004"  "chr7_2937190"   "chr8_23042229"  "chr10_5648124"  "chr10_29357400"
#[13] "chr11_12214430" "chr11_12371202" "chr11_25671857" "chr12_1328605"  "chr13_1134639"  "chr13_24316172"
#[19] "chr13_25016173" "chr15_14258937" "chr15_19053152" "chr16_16795729" "chr17_16655115" "chr19_22975760"
#[25] "chr21_3340668"  "chr22_20774031"

