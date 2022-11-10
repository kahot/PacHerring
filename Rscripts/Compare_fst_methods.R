#Comparison of Fst values calculated from different SFS 

#fst.fo<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS91_50kWindow")


#Fst using maf05 from Joe's vcf

Fsum<-data.frame(method=c("Joe_maf05","maf05","maf00"))

fst.joe<-read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt")
Fsum$pw07.vs.17[1]<-mean(fst.joe$Fst12,na.rm=T) #0.01080878
Fsum$pw91.vs.17[1]<-mean(fst.joe$Fst02,na.rm=T) #0.007911872
fst.joe<-read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt")
Fsum$pw91.vs.96[1]<-mean(fst.joe$Fst01,na.rm=T) #0.007098329



#Fst using maf00 of new_vcf
fst00<-read.delim("Data/new_vcf/angsd/fromVCF/fst_PWS07_PWS17_50kWindow_maf00")
Fsum$pw07.vs.17[3]<-mean(fst00$Nsites,na.rm=T)
#0.01220287
#Fst using maf05 of new_vcf
fst05<-read.delim("Data/new_vcf/angsd/fromVCF/maf05/fst_PWS07_PWS17_50kWindow_maf05")
Fsum$pw07.vs.17[2]<-mean(fst05$Nsites,na.rm=T)
#0.009756816

#91 vs. 17
fst00.2<-read.delim("Data/new_vcf/angsd/fromVCF/fst_PWS91_PWS17_50kWindow_maf00")
Fsum$pw91.vs.17[3]<-mean(fst00.2$Nsites,na.rm=T)
#0.008094017
fst05.2<-read.delim("Data/new_vcf/angsd/fromVCF/maf05/fst_PWS91_PWS17_50kWindow_maf05")
Fsum$pw91.vs.17[2]<-mean(fst05.2$Nsites,na.rm=T)
# 0.00627597

#91 vs. 96
fst00.3<-read.delim("Data/new_vcf/angsd/fromVCF/fst_PWS91_PWS16_50kWindow_maf00")
Fsum$pw91.vs.96[3]<-mean(fst00.3Nsites,na.rm=T)

fst05.3<-read.delim("Data/new_vcf/angsd/fromVCF/maf05/fst_PWS91_PWS96_50kWindow_maf05")
Fsum$pw91.vs.96[2]<-mean(fst05.3$Nsites,na.rm=T)
# 0.006126065

files5<-list.files("Data/new_vcf/angsd/fromVCF/maf05/",pattern="50kWindow_maf05")
files0<-list.files("Data/new_vcf/angsd/fromVCF/",pattern="50kWindow_maf00")

summary<-data.frame(method=c("maf00","maf05"))
for (i in 1:6){
    fst00<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/", files0[i]))
    fst05<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/maf05/",files5[i]))
    compare<-substr(files0[i], 5, 15)
    summary<-cbind(summary, c(mean(fst00$Nsites,na.rm=T), mean(fst05$Nsites,na.rm=T)))
    colnames(summary)[i+1]<-compare
}

#Add Joe's estimates
fst.joe<-read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt")
summary<-rbind(summary,c("Joe maf05", mean(fst.joe$Fst12,na.rm=T), mean(fst.joe$Fst01,na.rm=T),mean(fst.joe$Fst02,na.rm=T), NA,NA,NA))
fst.joe2<-read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt")
summary$PWS91_PWS96[3]<-mean(fst.joe2$Fst01,na.rm=T)
summary$PWS96_PWS07[3]<-mean(fst.joe2$Fst12,na.rm=T)
fst.joe3<-read.delim("Data/fst_pbs/maf05/Fst_pbs_50kb_win_10kb_step_folded_PWS96_PWS07_PWS17.txt")
summary$PWS96_PWS17[3]<-mean(fst.joe3$Fst02,na.rm=T)



source("Rscripts/BaseScripts.R")

fs<-melt(summary, id.vars="method", value.name = "Fst")
fs$Fst<-as.numeric(fs$Fst)
ggplot(fs,aes(x=variable, y=Fst, color=method))+ theme_bw()+
    geom_point(size=3)+xlab("")+theme(axis.text.x = element_text(angle=90))
ggsave("Output/Fst/Fst_estimates_angsd_comparison.pdf", width = 6, height = 4)


### Compare pi 