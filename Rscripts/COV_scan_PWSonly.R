#Find the regions with a high temporal covariance 
source("Rscripts/BaseScripts.R")


winsize<-c("100k","250k")
evens<-paste0("chr",seq(2,26, by=2))
cov.list<-list()

for (i in 1: length(winsize)){
    cov12<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/PWSonly_maf05_cov12_1996-1991_2006-1996_md7000_",winsize[i],"window.csv"), header = F)
    cov23<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/PWSonly_maf05_cov23_2017-2006_2006-1996_md7000_",winsize[i],"window.csv"), header = F)
    cov13<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/PWSonly_maf05_cov13_2017-2006_1996-1991_md7000_",winsize[i],"window.csv"), header = F)
    iv<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/PWSonly_MD7000_intervals_",winsize[i],"window.csv"), row.names = 1)
                 
    covs<-cbind(iv, cov12, cov23,cov13)
    colnames(covs)[4:6]<-c("cov12","cov23","cov13")
    covs$index=1:nrow(covs)
    
    covs$color<-"col1"
    covs$color[covs$chrom %in% evens]<-"col2"
    
    covs$cov12[is.nan(covs$cov12)]<-NA
    covs$cov12[is.infinite(covs$cov12)]<-NA
    
    cov.list[[i]]<-covs
    names(cov.list)[i]<-winsize[i]
    covsm<-melt(covs[,c("index","color","cov12","cov23","cov13")], id.vars = c("index", "color"))
    ymax<-max(covsm$value, na.rm=T)
    y<-min(covsm$value, na.rm=T)
    ymin<-ifelse (y<=-0.1,-0.1, y) 
    ggplot(covsm, aes(x=index, y=value, color=color))+
        facet_wrap(~variable, nrow=3)+
        geom_point(size=1, alpha=0.5)+
        theme_classic()+
        ylim(ymin,ymax)+
        scale_color_manual(values=c("gray70","steelblue"), guide="none")+
        ylab("Covariance")+xlab('Chromosome')+
        theme(axis.text.x = element_blank())+
        ggtitle(paste0("PWS"," ", winsize[i]," window"))
    ggsave(paste0("Output/COV/PWSonly_tempCovs_acrossGenome_",winsize[i], "Window.pdf"), width = 8, height = 8)    
}


#find how outliers overlap between different windows
cov12<-data.frame()
cov23<-data.frame()
cov13<-data.frame()

for (i in 1:length(winsize)){
    covs<-cov.list[[i]]
    covs<-covs[order(covs$cov12, decreasing=T),]
    n<-ceiling(nrow(covs)*0.01) #top1% region
    covs12_top<-covs[1:n,c(1:4)]
    covs12_top<-covs12_top[order(covs12_top$chrom, covs12_top$start),]
    covs12_top$window<-winsize[i]
    cov12<-rbind(cov12, covs12_top)
    
    covs<-covs[order(covs$cov13, decreasing=T),]
    n<-ceiling(nrow(covs)*0.01) #top1% region
    covs13_top<-covs[1:n,c(1:3,6)]
    covs13_top<-covs13_top[order(covs13_top$chrom, covs13_top$start),]
    covs13_top$window<-winsize[i]
    cov13<-rbind(cov13, covs13_top)
    
    covs<-covs[order(covs$cov23, decreasing=T),]
    n<-ceiling(nrow(covs)*0.01) #top1% region
    covs23_top<-covs[1:n,c(1:3,5)]
    covs23_top<-covs23_top[order(covs23_top$chrom, covs23_top$start),]
    covs23_top$window<-winsize[i]
    cov23<-rbind(cov23, covs23_top)
    
}

write.csv(cov12, "Output/COV/PWSonly_top1percent_outlier_regions.cov12.csv",row.names = F)
write.csv(cov23, "Output/COV/PWSonly_top1percent_outlier_regions.cov23.csv",row.names = F)
write.csv(cov13, "Output/COV/PWSonly_top1percent_outlier_regions.cov13.csv",row.names = F)


#1. Find the overlapping reigons between cov12, cov23 and over13

cov12<-read.csv("Output/COV/PWSonly_top1percent_outlier_cov12_overlap.csv")
cov23<-read.csv("Output/COV/PWSonly_top1percent_outlier_cov23_overlap.csv")
cov13<-read.csv("Output/COV/PWSonly_top1percent_outlier_cov13_overlap.csv")


# Overlapping regions between 100k and 250k windows
#dataframe to store overlapping regions
overlaps<-data.frame()
#covariacne time periods
covs<-c("cov12","cov23","cov13")
#response
ov<-c("Y","M")
for (i in 1:3){
    df<-read.csv(paste0("Output/COV/PWSonly_top1percent_outlier_",covs[i],"_overlap.csv"))
    df$Cov<-covs[i]
    colnames(df)[which(colnames(df)==covs[i])]<-"cov"
    overlaps<-rbind(overlaps, df[df$overlap %in% ov,])
}

write.csv(overlaps, "Output/COV/PWSonly_outliers_overlaps_100k.250k.csv")


#Create bed files for top 1% region (100k) and overlapping regions

#Bedfile
write.table(overlaps[overlaps$Cov=="cov12",1:3], "Output/COV/PWSonly_outliers_overlap_cov12.bed", quote = F, row.names = F, col.names = F,sep = "\t")
write.table(overlaps[overlaps$Cov=="cov13",1:3], "Output/COV/PWSonly_outliers_overlap_cov13.bed", quote = F, row.names = F, col.names = F,sep = "\t")
write.table(overlaps[overlaps$Cov=="cov23",1:3], "Output/COV/PWSonly_outliers_overlap_cov23.bed", quote = F, row.names = F, col.names = F,sep = "\t")
# remove the overlapping regions from the bed file manually. 


#Run snpEff



covs2<-covs[order(covs$cov12, decreasing=T),]
#top 1% (100k)
covs2_top<-covs2[1:72,c(1:4)]
covs2_top<-covs2_top[order(covs2_top$chrom, covs2_top$start),]
#create a bed file to find genes in these regions
write.table(covs2_top[,1:3], "Output/COV/pws_tempcov12_outliers.bed", quote = F, row.names = F, col.names = F,sep = "\t")

#top 1% 1m
covs2_top<-covs2[1:7,c(1:4)]
covs2_top<-covs2_top[order(covs2_top$chrom, covs2_top$start),]
write.table(format(covs2_top[,1:3],scientific=FALSE), "Output/COV/pws_tempcov12_1M_outliers.bed", quote = F, row.names = F, col.names = F,sep = "\t")


#top 1% 250k
covs2_top<-covs2[1:28,c(1:4)]
covs2_top<-covs2_top[order(covs2_top$chrom, covs2_top$start),]
write.table(format(covs2_top[,1:3],scientific=FALSE), "Output/COV/pws_tempcov12_250K_outliers.bed", quote = F, row.names = F, col.names = F,sep = "\t")



#over 0.01
cov2_over01<-covs2[covs2$cov12>=0.01,]
#6 regions
cov2_over01[,1:3]
#   chrom   start     end
#576  chr6 3.0e+06 4.0e+06
#541  chr4 3.1e+07 3.2e+07
#162 chr14 1.0e+07 1.1e+07
#72  chr11 9.0e+06 1.0e+07
#199 chr15 1.8e+07 1.9e+07
#121 chr12 2.8e+07 2.9e+07


# 2nd time period
covs$cov23[is.nan(covs$cov23)]<-NA
covs$cov23[is.infinite(covs$cov23)]<-NA
ggplot(covs, aes(x=index, y=cov23, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())+ylim(-0.34,0.4)
ggsave("Output/COV/PWS_tempCovs23_acrossGenome.pdf", width = 6, height = 3)    


covs3<-covs[order(covs$cov23, decreasing=T),]
#top 1%
covs3_top<-covs3[1:72,c(1:5)]
covs3_top<-covs3_top[order(covs3_top$chrom, covs3_top$start),]
write.table(covs3_top[,1:3], "Output/COV/pws_tempcov23_outliers.bed", quote = F, row.names = F, col.names = F,sep = "\t")

#1m
covs3_top<-covs3[1:7,c(1:5)]
covs3_top<-covs3_top[order(covs3_top$chrom, covs3_top$start),]
write.table(format(covs3_top[,1:3],scientific=FALSE), "Output/COV/pws_tempcov23_1M_outliers.bed", quote = F, row.names = F, col.names = F,sep = "\t")


