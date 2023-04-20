source("Rscripts/BaseScripts.R")



#pi
pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)
pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","TB17")
sss<-c("SS96","SS06","SS17")
y17<-c("TB17","PWS17","SS17","BC17","WA17","CA17")
comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-data.frame(t(comb4))
comb<-data.frame(rbind(comb1, comb2, comb3))
comb$Y2017<-"N"
comb4$Y2017<-"Y"
comb<-rbind(comb, comb4)

meanPi<-data.frame(pop.yr=pops)
for (i in 1:length(pops)){
    df<-read.csv(paste0("Output/Pi/",pops[i], "_Pi_pixy_per50kWindow.csv"), row.names =1 )
    meanPi$mean[i]<-mean(df$pi, na.rm=T)
    meanPi$pop[i]<-gsub("\\d.+","",pops[i])
}

pidiff<-data.frame(pop1=comb[,1], pop2=comb[,2])
for (i in 1:nrow(pidiff)){
    pidiff$pi_diff[i]<-abs(meanPi$mean[meanPi$pop.yr==comb[i,1]]-meanPi$mean[meanPi$pop.yr==comb[i,2]])
    pidiff$Y2017[i]<-comb[i,3]
}

#Fst files
fst1<-read.csv("Output/Fst/Mean_Fst_overTime_3pops_allcomp.csv", row.names = 1)
colnames(fst1)<-c("pop1","pop2","Fst","comp")
fst_mat<-read.csv("Output/SFS/Fst_matrix_2017_all.csv")
melted_cormat <- melt(fst_mat, na.rm = TRUE)
melted_cormat[melted_cormat==0]<-NA
f17<-melted_cormat[!is.na(melted_cormat$value),] 
colnames(f17)<-c("pop1","pop2","Fst")
f17$comp<-paste0(f17$pop1,"_",f17$pop2)

fst<-rbind(fst1, f17)

#merge the two files
pidiff$comp<-paste0(pidiff$pop1,"_",pidiff$pop2)

pidiff<-merge(pidiff, fst[, c("comp","Fst")], by="comp")

#Time-series
cor.test(pidiff$pi_diff[pidiff$Y2017=="N"], pidiff$Fst[pidiff$Y2017=="N"], method="spearman") 
#S = 476, p-value = 0.5934
#    rho 
#0.15 

#Y2017 between pops
cor.test(pidiff$pi_diff[pidiff$Y2017=="Y"], pidiff$Fst[pidiff$Y2017=="Y"], method="spearman") 
#S = 216, p-value = 0.01708
#    rho 
#0.6142857 

#all together
cor.test(pidiff$pi_diff, pidiff$Fst, method="spearman")
#S = 2194, p-value = 0.00429
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.5119021 

ggplot(pidiff, aes(x=pi_diff, y=Fst, color=Y2017))+
    geom_point(size=3, alpha=0.8)+
    theme_bw()+xlab("Diffrence in Pi")+
    scale_color_manual(values=c("steelblue",red))+
    annotate('text', x=0.00001, y=0.19, label="rho=0.51** (p=0.004)",hjust=0)
ggsave("Output/QC/Fst_pi.diff_correlation.png", width = 4.8, height = 4, dpi=300)


ggplot(pidiff[pidiff$Y2017=="N",], aes(x=pi_diff, y=Fst))+
    geom_point(size=3, color="steelblue")+
    theme_bw()+xlab("Diffrence in Pi")+
    annotate('text', x=0.00001, y=0.012, label="rho=0.15 (p=0.6)",hjust=0)
ggsave("Output/QC/Fst_pi.diff_correlation_no2017Cross.png", width = 4.5, height = 4.2, dpi=300)


ggplot(pidiff[pidiff$Y2017=="Y",], aes(x=pi_diff, y=Fst))+
    geom_point(size=3, color="steelblue")+
    theme_bw()+xlab("Diffrence in Pi")+
    annotate('text', x=0.00001, y=0.19, label="rho=0.61* (P=0.0017)",hjust=0)
ggsave("Output/QC/Fst_pi.diff_correlation_2017Cross.png", width = 4.5, height = 4.2, dpi=300)



#### Dxy (from calcDxy.R)
lens<-read.table("Data/new_vcf/chr_sizes.bed")
L<-sum(lens$V3)#725670187 bases
71914.4048090072/L



pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","TB17")
sss<-c("SS96","SS06","SS17")
y17<-c("TB17","PWS17","SS17","BC17","WA17","CA17")
comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-data.frame(t(comb4))
comb<-data.frame(rbind(comb1, comb2, comb3))
comb$Y2017<-"N"
comb4$Y2017<-"Y"
comb<-rbind(comb, comb4)

dxy<-data.frame(pop1=comb[,1], pop2=comb[,2])

for (i in 1: nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    df<-read.table(paste0("Scripts/calculateDxy/Dxy_persite_",pop1,"_",pop2,".txt"), header = T)
    dxy$mean_Dxy[i]<-sum(df$dxy, na.rm=T)/L
}

write.csv(dxy, "Output/Dxy/Dxy_from.calcDxy.csv", row.names = F)

dxy$comp<-paste0(dxy$pop1,"_",dxy$pop2)

pidiff<-merge(pidiff, dxy[,c("comp","mean_Dxy")], by="comp")

#Time-series
cor.test(pidiff$pi_diff[pidiff$Y2017=="N"], pidiff$mean_Dxy[pidiff$Y2017=="N"], method="spearman") 
#S = 726, p-value = 0.2827
#    rho 
#-0.2964286 

#Y2017 between pops
cor.test(pidiff$pi_diff[pidiff$Y2017=="Y"], pidiff$mean_Dxy[pidiff$Y2017=="Y"], method="spearman") 
#S = 404, p-value = 0.3138
#    rho 
#0.2785714 

#all together
cor.test(pidiff$pi_diff, pidiff$mean_Dxy, method="spearman")
#S = 3692, p-value = 0.3434
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.1786429 

ggplot(pidiff, aes(x=pi_diff, y=mean_Dxy, color=Y2017))+
    geom_point(size=3, alpha=0.8)+
    theme_bw()+xlab("Diffrence in Pi")+
    scale_color_manual(values=c("steelblue",red))+ylab("Dxy")
ggsave("Output/QC/Dxy_pi.diff_correlation.png", width = 4.8, height = 4, dpi=300)




#### Dxy (calculated manually from mafs)

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)
pwss<-pops[grep("PWS",pops)]
tbs<-pops[grep("TB",pops)]
sss<-pops[grep("SS",pops)]
y17<-pops[grep("17",pops)]
comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-data.frame(t(comb4))
comb<-data.frame(rbind(comb1, comb2, comb3))
comb$Y2017<-"N"
comb4$Y2017<-"Y"
comb<-rbind(comb, comb4)

dxy<-data.frame(pop1=comb[,1], pop2=comb[,2])
for (i in 1: nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    dxy$Dxy[i]<-mean(rol_win$dxy, na.rm=T)
}



#### Read depth vs. Pi comparison ####

## look at the correlation between read coverage and pi

pops<-c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17")
year<-c(1991,1996,2006,2017,1991,1996,2007,2017,1996,2006,2017,2017,2017,2017)
meanPi<-data.frame(pop.yr=pops, year=year)
for (i in 1:length(pops)){
    df<-read.csv(paste0("Output/Pi/",pops[i], "_Pi_pixy_per50kWindow.csv"), row.names =1 )
    meanPi$mean[i]<-mean(df$pi, na.rm=T)
    meanPi$pop[i]<-gsub("\\d.+","",pops[i])
}

### mean number of reads per sample (using samtools to get the number of mapped reads)
byPopYr<-read.csv("Output/QC/Bam_metrics_summary.csv", row.names = 1)
colnames(byPopYr)<-c("pop.yr","depth")
data<-merge(byPopYr, meanPi, by="pop.yr")

write.csv(data, "Output/read.depth_pi.csv", row.names = F)

cor.test(data$depth, data$mean, method = "spearman")
#Spearman's rank correlation rho
#
#data:  data$depth and data$mean
#S = 114, p-value = 0.002992
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7494505 

ho<-read.csv("Output/Stats_window/Heterozygosity_pop_summary.csv", row.names = 1)
data<-merge(data, ho[,c("pop","Ho_mean")], by.x = "pop.yr", by.y = "pop")

#cor.test(data$depth, data$Ho_mean, method = "spearman")
#data:  data$depth and data$Ho_mean
#S = 36, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.9208791 
data$pop<-factor(data$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))
ggplot(data, aes(x=depth, y=mean, color=pop))+
    geom_point(size=3)+
    xlab("Mean read depth")+
    ylab(expression(paste("Mean ", pi)))+
    theme_bw()+
    annotate('text', x=1.28, y=0.00325, label="rho = 0.75**", size=3)+
    scale_color_manual(values=cols)+ theme(legend.title = element_blank())
ggsave("Output/QC/Cor_plot_depth.vs.Pi.png", width = 5.2, height = 4, dpi=300)

ggplot(data, aes(x=depth, y=Ho_mean,color=pop))+
    geom_point(size=3)+
    xlab("Mean read depth")+
    ylab("Mean Ho")+
    theme_bw()+  scale_color_manual(values=cols)+
    annotate('text', x=1.28, y=0.0041, label="rho = 0.92***", size=3)+
    theme(legend.title = element_blank())
ggsave("Output/QC/Cor_plot_depth.vs.Ho.png", width = 5.2, height = 4, dpi=300)


# Read depth vs. Pi from ANGSD with downsampled data

thetas<-read.csv("Output/Pi/Pi_windows_all_pop.csv", row.names = 1)
Pi.mean<-data.frame(aggregate(thetas[,"pi"], by=list(thetas$pop), mean, na.rm=T))
colnames(Pi.mean)<-c("pop","Pi_downsampled")

data<-merge(data, Pi.mean,  by.x = "pop.yr", by.y = "pop")


cor.test(data$depth, data$Pi_downsampled, method = "spearman")
#data:  data$depth and data$mean
#S = 158, p-value = 0.01371
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.6527473 

data$pop<-factor(data$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))

ggplot(data, aes(x=depth, y=Pi_downsampled, color=pop))+
    geom_point(size=3)+
    xlab("Mean read depth")+
    ylab(expression(paste("Mean ", pi)))+
    theme_bw()+
    annotate('text', x=0.9, y=0.00365, label="rho = 0.65, p=0.014", size=4, hjust=0)+
    scale_color_manual(values=cols)+ theme(legend.title = element_blank())
ggsave("Output/QC/Mean_depth.vs.meanPi_downsampled.png", width = 6, height = 4, dpi=300)


 




