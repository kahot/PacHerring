#import temporal covariance values generated from cvtkpy
library(ggplot2)

##### MD7000 MAF0.01 100k windows
pops<-c("PWS","TB","SS")
covs<-data.frame()
Variance<-data.frame()

for (p in 1: length(pops)){
    #covariance output file
    cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf01_temp_cov_matrix_",pops[p],"_100k.csv"))
    
    cov<-cov[,-1]
    
    ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf01_",pops[p],"_Cov_CIs_bootstrap5000.csv"))
    ci<-ci[,-1]
    
    #reshape the matrix
    mat1<-cov[1:3,]
    mat2<-cov[4:6,]
    
    covdf<-data.frame()
    k=1
    for (i in 1:nrow(mat1)){
        for (j in 1:ncol(mat1)){
            covdf[k,1]<-mat2[i,j]
            covdf[k,2]<-mat1[i,j]
            k=k+1
        }
    }
    colnames(covdf)<-c("label","value")
    covdf$value<-as.numeric(covdf$value)
    covar<-covdf[grep("cov",covdf$label),]
    vars<-covdf[grep("var",covdf$label),]
    
    #remove the redundant values
    #assign the starting time period and covering period values
    if (pops[p]!="SS") covar<-covar[!duplicated(covar[, 2]),] 
    if (pops[p]=="SS") covar<-covar[c(1,2,4),]
    
    #assign the starting time period and covering period values
    covar$year<-c(1,2,2)
    covar$series<-c("1991","1991","1996")
    
    vars$year<-c(1,2,2)
    vars$series<-c("1991","1991","1996")
    
    #assign population name
    covar$location<-pops[p]
    vars$location<-pops[p]
    
    #attach ci info
    covar$ci_l<-c(ci[1,2], ci[1,3],ci[2,3])
    covar$ci_u<-c(ci[4,2], ci[4,3],ci[5,3])
    
    #combine in to one matrix
    covs<-rbind(covs, covar)
    Variance<-rbind(Variance, vars)
}

covs$ci_l<-as.numeric(covs$ci_l)
covs$ci_u<-as.numeric(covs$ci_u)

ggplot(data=covs, aes(x=year, y=value, color=location, shape=series, group=interaction(location, series)))+
    geom_point(size=3, position=position_dodge(width = 0.1,preserve ="total"))+
    #geom_errorbar(data=covs, aes(x=year, y=value, ymin=ci_l, ymax=ci_u), width=.2, size=.2, position=position_dodge(width = 0.1,preserve ="total"))+
    geom_line(data=covs, aes(x=year, y=value,color=location, group=interaction(location, series)), position=position_dodge(width = 0.1,preserve ="total"))+
    ylab("Covariance")+xlab('')+theme_classic()+
    theme(axis.text.x = element_blank(),legend.title = element_blank())+
    geom_hline(yintercept = 0,color="gray70", size=0.3)+
    geom_errorbar(aes(ymin=ci_l, ymax=ci_u), width=.2, size=.2, position=position_dodge(width = 0.1,preserve ="total"))+
    scale_shape_manual(values=c(16,17),labels=c("1991-","1996-"))+
    scale_x_continuous(breaks = c(1,2))
ggsave("Output/COV/MD7000_maf01_Cov_overtime_CI_100k.window.pdf",width = 4.7, height = 3)

## create manhattan plot for maf01 (100k)
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_maf01_100kwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_maf01_100kwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_maf01_100kwindow.csv", header = F)

iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf01_intervals_100Kwindow.csv", row.names = 1)

covs<-cbind(iv, cov12, cov23,cov13)
colnames(covs)[4:6]<-c("cov12","cov23","cov13")
covs$index=1:nrow(covs)

evens<-paste0("chr",seq(2,26, by=2))
covs$color<-"col1"
covs$color[covs$chrom %in% evens]<-"col2"

covs$cov12[is.nan(covs$cov12)]<-NA
covs$cov12[is.infinite(covs$cov12)]<-NA
write.csv(covs,"Output/COV/PWS_tempCovs_md7000_maf01_100k.csv")

ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+ylim(-0.1,0.1)+xlab('Chromosome')+
    theme(axis.text.x = element_blank())
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_MD7000_maf01.pdf", width = 6, height = 3)    


### 10kb-window


## create manhattan plot for maf01 (100k)
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_maf01_10kwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_maf01_10kwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_maf01_10kwindow.csv", header = F)

iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf01_intervals_10Kwindow.csv", row.names = 1)

covs<-cbind(iv, cov12, cov23,cov13)
colnames(covs)[4:6]<-c("cov12","cov23","cov13")
covs$ch<-gsub("chr","", covs$chrom)
covs$ch<-as.integer(covs$ch)
covs<-covs[order(covs$ch),]

evens<-paste0("chr",seq(2,26, by=2))
covs$color<-"col1"
covs$color[covs$chrom %in% evens]<-"col2"

covs$cov12[is.nan(covs$cov12)]<-NA
covs$cov12[is.infinite(covs$cov12)]<-NA
covs$cov23[is.nan(covs$cov23)]<-NA
covs$cov23[is.infinite(covs$cov23)]<-NA
covs$cov13[is.nan(covs$cov13)]<-NA
covs$cov13[is.infinite(covs$cov13)]<-NA

covs$index=1:nrow(covs)

write.csv(covs,"Output/COV/PWS_tempCovs_md7000_maf01_10k.csv")

ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_MD7000_maf01_10k.pdf", width = 6, height = 3)    

#zoomed
ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())+
    ylim(-0.2,0.2)
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_MD7000_maf01_10k_zoomedIn.pdf", width = 6, height = 3)    

ggplot(covs, aes(x=index, y=cov23, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())+ylim(-1.5,1)
ggsave("Output/COV/PWS_tempCovs23_acrossGenome_MD7000_maf01_10k_zoomed.pdf", width = 6, height = 3)    


ggplot(covs, aes(x=index, y=cov13, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank()) #+ylim(-1.5,1)
ggsave("Output/COV/PWS_tempCovs13_acrossGenome_MD7000_maf01_10k.pdf", width = 6, height = 3)    



# Where are the largest cov regions?

covs2<-covs[order(covs$cov12, decreasing=T),]
#top 1%
covs2_top<-covs2[1:70,c(1:4)]
covs2_top<-covs2_top[order(covs2_top$chrom, covs2_top$start),]
#create a bed file to find genes in these regions
write.table(covs2_top[,1:3], "Output/COV/pws_tempcov12_outliersMD7000_maf01.bed", quote = F, row.names = F, col.names = F,sep = "\t")


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

## Run snpEff
# Create a new vcf file containing the loci in p07_loc.
## at terminal
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_tempcov12_outliers.bed --out Output/COV/annotation/PWS_cov12_outlier --recode --keep-INFO-all
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_tempcov23_outliers.bed --out Output/COV/annotation/PWS_cov23_outlier --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12_outlier.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12 > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov12_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23_outlier.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23 > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov23_outlier.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov12_outlier.vcf > Output/COV/annotation/PWS_cov12_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov23_outlier.vcf > Output/COV/annotation/PWS_cov23_annotation


#snpEff results
compa<-c("cov12","cov23")

for (i in 1:2){
    df<-read.table(paste0("Output/COV/annotation/PWS_",compa[i],"_annotation"), header = F)
    annotations<-data.frame()
    for (j in 1: nrow(df)){
        anns<-unlist(strsplit(df$V4[j], "\\|"))
        anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
        annotations<-rbind(annotations, anns)
    }     
    
    colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
    Ano<-cbind(df[,1:3], annotations)
    colnames(Ano)[1:3]<-c("chr","pos","AF")
    #remove the duplicated annotations for deeper digging
    remove<-!duplicated(annotations)
    Ano2<-Ano[remove,]
    write.csv(Ano2, paste0("Output/COV/annotation/PWS_",compa[i],"_outlier_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids1<-geneids[nchar(geneids)<=18]
    geneids1<-unique(geneids1)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_outlier_geneid_list1.txt"))
    cat(paste0(geneids1,"; "))
    sink(NULL)
    
    #split the intergenic ids into two
    geneids2<-geneids[nchar(geneids)>18]
    ids2<-unlist(str_split(geneids2, "-", 2))
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-c(geneids, ids2)
    geneids<-unique(geneids)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_outlier_geneid_list2.txt"))
    cat(paste0(geneids,"; "))
    sink(NULL)
    
    #gene names
    genenames<-c(Ano2$Gene_name, Ano2$Gene_name2)
    genenames<-unique(genenames)
    genenames2<-genenames[-grep("ENSCHAG",genenames)]
    genenames2<-genenames2[-grep("si\\:",genenames2)]
    long<-genenames2[grep("\\-",genenames2)]
    longids<-unlist(str_split(long, "-", 2))
    genenames2<-genenames2[-grep("\\-",genenames2)]
    genenames2<-c(genenames2, longids)
    genenames2<-unique(genenames2)
    
    write.table(genenames2, paste0("Output/COV/annotation/PWS_",compa[i],"_outlier_genenames_list_test.txt"), quote=F, row.names = F, col.names = F)
}





##### Tiled covs    
pw<-read.csv("~/Projects/Pacherring_Vincent/notebooks/PWS_cov12_1996-1991_2006-1996.csv")
pw[pw=="NaN"]<-NA
pw$pos<-1:nrow(pw)
plot(pw$pos, pw$nan, pch=".")    
pw1<-pw[!is.na(pw$nan),]
colnames(pw1)[1]<-"cov"
hist(pw1$cov, xlim=c(-0.04,0.04), breaks=30)

ggplot(pw1, aes(x=cov))+
    geom_histogram(data=pw1,  color="gray60", alpha = 0.5, binwidth =0.005 ) +
    xlim(-0.2,0.2)+
    theme_classic()

