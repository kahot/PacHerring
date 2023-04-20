# Why Fst values differ between files obtained using 3DSFS

fst1 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_CA17_PWS17_50kWindow")
colnames(fst1)[5:7]<-c("Fst_BC.CA","Fst_BC.PWS","Fst_CA.PWS")

fst2 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_CA17_WA17_50kWindow")
colnames(fst2)[5:7]<-c("Fst_BC.CA","Fst_BC.WA","Fst_CA.WA")
mean(fst1$Fst_BC.CA) #0.02984197
mean(fst1$Fst_BC.CA) #0.02984197
#BC.vs.CA is the same

fst3 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_SS17_TB17_50kWindow")
colnames(fst3)[5:7]<-c("Fst_BC.SS","Fst_BC.TB","Fst_SS.TB")

fst4 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_CA17_PWS17_SS17_50kWindow")
colnames(fst4)[5:7]<-c("Fst_CA.PWS","Fst_CA.SS","Fst_PWS.SS")

mean(fst4$Fst_CA.PWS) #0.03114319
mean(fst1$Fst_CA.PWS) #0.03138774

fst5 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_CA17_PWS17_TB17_50kWindow")
colnames(fst5)[5:7]<-c("Fst_CA.PWS","Fst_CA.TB","Fst_PWS.TB")

mean(fst5$Fst_CA.PWS) #0.03124173

fst6 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_PWS17_SS17_WA17_50kWindow")
colnames(fst6)[5:7]<-c("Fst_PWS.SS","Fst_PWS.WA","Fst_SS.WA")

mean(fst6$Fst_PWS.SS) #0.008251199
mean(fst4$Fst_PWS.SS) #0.008366741

fst7 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_SS17_TB17_WA17_50kWindow")
colnames(fst7)[5:7]<-c("Fst_SS.TB","Fst_SS.WA","Fst_TB.WA")


merged67<-merge(fst6[,c(2,3,4,5,6,7)], fst7[,c(2,3,4,5,6,7)],by=c("chr","midPos"))


mean(fst6$Fst_SS.WA) # 0.009572484
mean(fst7$Fst_SS.WA) #0.009179432

mean(fst7$Fst_SS.TB) # 0.1006478
mean(fst3$Fst_SS.TB) #0.1285945

mean(fst7$PBS2) #0.01149243
mean(fst3$PBS2) #0.1455083


fst<-read.csv("Output/SFS/Fst_window_year2017_combined.csv")
mean(fst$Fst_SS.TB) #0.1286036
mean(fst$Fst_SS.WA) # 0.009572478

#Per-site Fst comparison

pfst1<-read.delim("Data/new_vcf/angsd/fromVCF/fst_folded__BC17_SS17_TB17.persite.txt")
mean(pfst1$X0.000218)
persite1<-read.table("Data/new_vcf/angsd/fromVCF/fst_BC_SS_TB.txt", header=F)
persite2<-read.table("Data/new_vcf/angsd/fromVCF/fst_SS_TB_WA.txt", header=F)
#SS vs. TB

mean(persite1$V5)
mean(persite2$V3)

#1. 0.0001610366


#Calculate 2D Fst for all combinations

fst_2d<-read.delim("Data/new_vcf/angsd/fromVCF/2D/fst_SS17_TB17_50kWindow_maf00")
mean(fst_2d[,4]) #0.1369242

#TB pops

tb1<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/3D/fst_TB91_TB96_TB06_maf00_50kWindow"))
tb2<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/3D/fst_TB91_TB96_TB17_maf00_50kWindow"))
tb3<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/3D/fst_TB91_TB06_TB17_maf00_50kWindow"))

#91 vs. 96
mean(tb1$Fst01) # 0.007100834
mean(tb2$Fst01) # 0.007075912

#91. vs 17
mean(tb2$Fst02) #0.007508833
mean(tb3$Fst02) #0.007564142



pairs<-c()
fst<-data.frame()
for (i in 1:3){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    pop3<-comb[i,3]
    
    df<-read.delim(paste0("../Data/new_vcf/angsd/fromVCF/3D/fst_",pop1,"_",pop2,"_",pop3,"_maf00_50kWindow"))
    
    pair<-paste0(pop1,".vs.",pop2)
    if (!(pair %in% pairs)){
        df2<-df[,c("chr","midPos","Fst01")]
        df2$pop<-pair
        colnames(df2)[colnames(df2)=="Fst01"]<-"Fst"
        df2$ch=as.integer(gsub("chr","", df2$chr))
        df2<-df2[order(df2$ch),]
        df2$loc<-1:nrow(df2)
        fst<-rbind(fst, df2)
        pairs<-c(pairs, pair)
    }
    
    pair<-paste0(pop1,".vs.",pop3)
    if (!(pair %in% pairs)){
        df2<-df[,c("chr","midPos","Fst02")]
        df2$pop<-pair
        colnames(df2)[colnames(df2)=="Fst02"]<-"Fst"
        df2$ch=as.integer(gsub("chr","", df2$chr))
        df2<-df2[order(df2$ch),]
        df2$loc<-1:nrow(df2)
        fst<-rbind(fst, df2)
        pairs<-c(pairs, pair)
    }
    
    pair<-paste0(pop2,".vs.",pop3)
    if (!(pair %in% pairs)){
        df2<-df[,c("chr","midPos","Fst12")]
        df2$pop<-pair
        colnames(df2)[colnames(df2)=="Fst12"]<-"Fst"
        df2$ch=as.integer(gsub("chr","", df2$chr))
        df2<-df2[order(df2$ch),]
        df2$loc<-1:nrow(df2)
        fst<-rbind(fst, df2)
        pairs<-c(pairs, pair)
    }
}

evens<-paste0("chr",seq(2,26, by=2))
#Plot Fst values across Genome
fst$color<-"col1"
fst$color[fst$chr %in% evens]<-"col2"
fst$pop<-factor(fst$pop, levels=unique(fst$pop))


