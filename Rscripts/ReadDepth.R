# Check read depth across genomes to see if there are any regions with higher coverage that stands out

library(windowscanr)
library(ggplot2)
library(gridExtra)

snps <- read.table("~/Projects/pac_herring_joe/slurm_scripts/fastq_to_vcf/combine_gvcfs/chr1_1_out.INFO", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

mean(snps$DP) #1825
median(snps$DP) #1736
sd(snps$DP) #1471
median(snps$DP)+2*sd(snps$DP) #4705.837

snps<-snps[,c("CHROM","POS","DP"),]
colnames(snps)<-c("V1","V2","V3")
write.table(snps, "Data/new_vcf/Raw_AN/chr1.1.info", row.names=FALSE, col.names = FALSE, quote = FALSE)



AN.summary<-data.frame(chr=c(paste0("chr",1:26, ".1"), paste0("chr",1:26, ".2")))
AN.summary$mean<-''

plots1<-list()
plots2<-list()
for (i in 2: 26){
    df1<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".1.info"), header = FALSE, stringsAsFactors = FALSE)
    AN.summary$mean[AN.summary$chr==paste0("chr",i,".1")]<-mean(df1$V3,na.rm=T )    
    AN.summary$median[AN.summary$chr==paste0("chr",i,".1")]<-median(df1$V3,na.rm=T ) 
    AN.summary$sd[AN.summary$chr==paste0("chr",i,".1")]<-sd(df1$V3,na.rm=T) 
    df2<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".2.info"), header = FALSE, stringsAsFactors = FALSE)
    AN.summary$mean[AN.summary$chr==paste0("chr",i,".2")]<-mean(df2$V3,na.rm=T )    
    AN.summary$median[AN.summary$chr==paste0("chr",i,".2")]<-median(df2$V3,na.rm=T ) 
    AN.summary$sd[AN.summary$chr==paste0("chr",i,".2")]<-sd(df2$V3,na.rm=T ) 


    win1<-winScan(x = df1, position = "V2", values = "V3", 
                  win_size = 100000, win_step = 10000,
                  funs = c("mean"),cores=10)
    win2<-winScan(x = df2, position = "V2", values = "V3", 
                  win_size = 100000, win_step = 10000,
                  funs = c("mean"),cores=10)
    win2<-win2[!is.nan(win2$V3_mean),]
    #
    plots1[[i]]<-ggplot(win1, aes(x=win_mid, y=V3_mean))+
        geom_point(color="steelblue", size=0.2)+
        ylab("read deapth")+xlab("")+
        ggtitle(paste0("chr",i,".1"))
    plots2[[i]]<-ggplot(win2, aes(x=win_mid, y=V3_mean))+
        geom_point(color="steelblue", size=0.2)+
        ylab("read deapth")+xlab("")+
        ggtitle(paste0("chr",i,".2"))
    print(i)
}


AN.summary$max_cutoff<-AN.summary$median+2*AN.summary$sd


#dir.create("Output/Depth/")
write.csv(AN.summary,"Output/Depth/Read_depth_allsamples_median.csv")


#plot each chromosome together

#rearragne the two lists
Plots<-list()
k=1
for (i in 1:26){
    Plots[[k]]<-plots1[[i]]
    k=k+1
    Plots[[k]]<-plots2[[i]]
    k=k+1
}
pdf("Output/Depth/readdepth_allSamples_100kscan.pdf", width=20, height=25)
do.call(grid.arrange, c(Plots, ncol=4))
dev.off()

mean(AN.summary$max_cutoff)
#[1] 6851.77


#look for high coverage areas 

high_depth_sites<-data.frame()
for (i in 1: 26){
    df1<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".1.info"), header = FALSE, stringsAsFactors = FALSE)
    df1<-df1[order(df1$V3, decreasing = T),]
    print(df1$V3[1])
    df1_high<-df1[df1$V3>10000,]
    high_depth_sites<-rbind(high_depth_sites, df1_high)
    
    df2<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".2.info"), header = FALSE, stringsAsFactors = FALSE)
    df2<-df2[order(df2$V3, decreasing = T),]
    print(df2$V3[1])
    df2_high<-df2[df2$V3>10000,]
    high_depth_sites<-rbind(high_depth_sites, df2_high)
}

colnames(high_depth_sites)<-c("chr","pos","ReadDepth")
write.csv(high_depth_sites, "Output/Depth/HighDepthSites.csv")

#group them by regions
library(dplyr)

gen.size<-read.table("Data/vcfs/chr_sizes.txt")
colnames(gen.size)<-c("chr", "size")
groups<-list()
for (i in 1:26){
    ch<-paste0("chr",i)
    df<-high_depth_sites[high_depth_sites$chr==ch,]
    df<-df%>% mutate(region=cut(pos, breaks=seq(1, gen.size$size[i]+1000000, 100000)))
    length(unique(df$region))
    sum<-data.frame(table(df$region))
    sum<-sum[sum$Freq>0,]
    
    
    pos<-high_depth_sites$pos[i]
    groups[i]
    
}


gt<-read.csv("Output/chr7/DP7000/CA17_gtinfo.csv", row.names = 1)


### Plot depth distribution
hist_log<-list()
histgrams<-list()
for (i in 1:26){
    df1<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".1.info"), header = FALSE, stringsAsFactors = FALSE)
    df2<-read.table(paste0("Data/new_vcf/raw_DP/chr",i,".2.info"), header = FALSE, stringsAsFactors = FALSE)
    df<-rbind(df1, df2)
    df$log<-log10(df$V3)
    hist_log[[i]]<-ggplot(df, aes(x=log))+
        geom_histogram()+xlab("log10(read depth)")+
        theme_bw()+ggtitle(paste0("chr",i))+
        geom_vline(xintercept=c(2.778151, 3.845098), color="blue",size=.3 )
    
    histgrams[[i]]<-ggplot(df, aes(x=V3))+
        geom_histogram()+xlab("read depth (max 10,000")+
        theme_bw()+ggtitle(paste0("chr",i))+xlim(0,10000)+
        geom_vline(xintercept=c(600, 7000), color="blue",size=.3 )
    
    
}

pdf("Output/Depth/Read_depth_hist_perChrom_log.pdf", height = 10, width = 18)
do.call(grid.arrange, c(hist_log, ncol=7))
dev.off()

pdf("Output/Depth/Read_depth_hist_perChrom.pdf", height = 10, width = 18)
do.call(grid.arrange, c(histgrams, ncol=7))
dev.off()

### 