#Estimate effective population size
#Nest/PoolSeq package is used. Ref: doi: 10.1534/genetics.116.191197
library(poolSeq)

#Read the allele freq data for PWS

pws1<-read.table("Data/new_vcf/AF/PWSonly91_maf05_af.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws2<-read.table("Data/new_vcf/AF/PWSonly96_maf05_af.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws3<-read.table("Data/new_vcf/AF/PWSonly07_maf05_af.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws4<-read.table("Data/new_vcf/AF/PWSonly17_maf05_af.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))

#combine AF for all years
pws<-cbind(pws1, pws2[,3],pws3[,3],pws4[,3])
colnames(pws)<-c("chr","pos","F0","F1","F2","F3")

#Read depth information
p_reads<-data.frame()
for (i in 1:26){
    df<-read.csv(paste0("Output/CNV/chr",i,"_330k_depth.csv"), row.names = 1)
    df_p91<-df[,grep("PWS91", colnames(df))]
    df_p96<-df[,grep("PWS96", colnames(df))]
    df_p07<-df[,grep("PWS07", colnames(df))]
    df_p17<-df[,grep("PWS17", colnames(df))]
    reads<-data.frame(F0=rowSums(df_p91, na.rm = T), F1=rowSums(df_p96, na.rm = T),F2=rowSums(df_p07, na.rm = T),F3=rowSums(df_p17, na.rm = T))
    p_reads<-rbind(p_reads,reads)
}
write.csv(p_reads,"Output/Depth/PWS_read_depth.csv")

#Find SNPs with extreme values and uninformative loci and remove them
retain<-checkSNP(pws[,"F0"],pws[,"F3"],p_reads[,"F0"], p_reads[,"F3"])
length(retain[retain==F])

#filtered the snp dataset
pws_filtered<-pws[retain,]
pwsS_filtered<-p_reads[retain,]

#Look at F0 and F1
retain<-checkSNP(pws_filtered[,"F0"],pws_filtered[,"F1"],pwsS_filtered[,"F0"], pwsS_filtered[,"F1"])
length(retain[retain==F]) #0

retain<-checkSNP(pws_filtered[,"F0"],pws_filtered[,"F2"],pwsS_filtered[,"F0"], pwsS_filtered[,"F2"])
length(retain[retain==F]) #0

retain<-checkSNP(pws_filtered[,"F1"],pws_filtered[,"F2"],pwsS_filtered[,"F1"], pwsS_filtered[,"F2"])
length(retain[retain==F]) #0

retain<-checkSNP(pws_filtered[,"F2"],pws_filtered[,"F3"],pwsS_filtered[,"F2"], pwsS_filtered[,"F3"])
length(retain[retain==F]) #0

write.csv(pws_filtered, "Output/COV/PWS_snps_filtered_forNeestimate_vcfAF.csv")


#### 

pws_filtered<-read.csv("Output/COV/PWS_snps_filtered_forNeestimate_angsdAF.csv", row.names = 1)
#Estimate genome-wide Ne
methods<-c("W.planI","W.planII","JR.planI","JR.planII","P.planI","P.planII","P.alt.1step.planII")

pws_Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    pws_Ne$Ne[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                  method=methods[i], Ncensus=1000,poolSize=c(58,56))
    pws_Ne$Ne_10000[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                                   method=methods[i], Ncensus=10000,poolSize=c(58,56))
    
}

write.csv(pws_Ne,"Output/COV/Ne_estimation_PWS91-17_vcfAF.csv")

Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne$Ne01_t1[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F1"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F1"], t=1, 
                         method=methods[i], Ncensus=10000,poolSize=c(58,72))
    Ne$Ne12_t2[i]<-estimateNe(p0=pws_filtered[,"F1"], pt=pws_filtered[,"F2"], cov0=pwsS_filtered[,"F1"], covt=pwsS_filtered[,"F2"], t=1.8, 
                              method=methods[i], Ncensus=10000,poolSize=c(72,46))
    Ne$Ne23_t2[i]<-estimateNe(p0=pws_filtered[,"F2"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F2"], covt=pwsS_filtered[,"F3"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(72,46))
}
write.csv(Ne, "Output/COV/Ne_estimation_PWS_eachTimePeriod_vcfAF.csv")


# Stika Sound
S1 <- read.delim("Data/new_vcf/AF/SS96.mafs")
S2 <- read.delim("Data/new_vcf/AF/SS06.mafs")
S3 <- read.delim("Data/new_vcf/AF/SS17.mafs")

#AF table for all years
ss<-cbind(S1[,c(1,2,6)], S2[,c(6)],S3[,c(6)])
colnames(ss)[3:5]<-c("F1","F2","F3")

#Read depth information
s_reads<-data.frame()
for (i in 1:26){
    df<-read.csv(paste0("Output/CNV/chr",i,"_depth.csv"), row.names = 1)
    df_96<-df[,grep("SS96", colnames(df))]
    df_06<-df[,grep("SS06", colnames(df))]
    df_17<-df[,grep("SS17", colnames(df))]
    reads<-data.frame(F1=rowSums(df_96, na.rm = T),F2=rowSums(df_06, na.rm = T),F3=rowSums(df_17, na.rm = T))
    s_reads<-rbind(s_reads,reads)
}


# TB
T1 <- read.delim("Data/new_vcf/AF/TB91.mafs")
T2 <- read.delim("Data/new_vcf/AF/TB96.mafs")
T3 <- read.delim("Data/new_vcf/AF/TB06.mafs")
T4 <- read.delim("Data/new_vcf/AF/TB17.mafs")

#AF tables for all years
tb<-cbind(T1[,c(1,2,6)], T2[,c(6)],T3[,c(6)],T4[,c(6)])
colnames(tb)[3:6]<-c("F0","F1","F2","F3")

#Read depth information
t_reads<-data.frame()
for (i in 1:26){
    df<-read.csv(paste0("Output/CNV/chr",i,"_depth.csv"), row.names = 1)
    df_91<-df[,grep("TB91", colnames(df))]
    df_96<-df[,grep("TB96", colnames(df))]
    df_06<-df[,grep("TB06", colnames(df))]
    df_17<-df[,grep("TB17", colnames(df))]
    reads<-data.frame(F0=rowSums(df_91, na.rm = T), F1=rowSums(df_96, na.rm = T),F2=rowSums(df_06, na.rm = T),F3=rowSums(df_17, na.rm = T))
    t_reads<-rbind(t_reads,reads)
}

# run checkSNP on SS pop
remove<-checkSNP(ss[,"F1"],ss[,"F3"],s_reads[,"F1"], s_reads[,"F3"])
#filtered the snp dataset
ss_filtered<-ss[remove,]
ssS_filtered<-s_reads[remove,]

ss_Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    ss_Ne$Ne[i]<-estimateNe(p0=ss_filtered[,"F1"], pt=ss_filtered[,"F3"], cov0=ssS_filtered[,"F1"], covt=ssS_filtered[,"F3"], t=2, 
                             method=methods[i], Ncensus=1000,poolSize=c(78,64))
    ss_Ne$Ne_10000[i]<-estimateNe(p0=ss_filtered[,"F1"], pt=ss_filtered[,"F3"], cov0=ssS_filtered[,"F1"], covt=ssS_filtered[,"F3"], t=2, 
                                   method=methods[i], Ncensus=10000,poolSize=c(78,64))
    
}
write.csv(ss_Ne,"Output/COV/Ne_estimation_SS96-17_angsdAF.csv")

Ne_s<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne_s$Ne12_t2[i]<-estimateNe(p0=ss_filtered[,"F1"], pt=ss_filtered[,"F2"], cov0=ssS_filtered[,"F1"], covt=ssS_filtered[,"F2"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(78,64))
    Ne_s$Ne23_t2[i]<-estimateNe(p0=ss_filtered[,"F2"], pt=ss_filtered[,"F3"], cov0=ssS_filtered[,"F2"], covt=ssS_filtered[,"F3"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(78,64))
    
}
write.csv(Ne_s, "Output/COV/Ne_estimation_SS_eachTimePeriod_angsdAF.csv")


#TB
remove<-checkSNP(tb[,"F1"],tb[,"F3"],t_reads[,"F1"], t_reads[,"F3"])

tb_filtered<-tb[remove,]
tbS_filtered<-t_reads[remove,]

tb_Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    tb_Ne$Ne[i]<-estimateNe(p0=tb_filtered[,"F0"], pt=tb_filtered[,"F3"], cov0=tbS_filtered[,"F0"], covt=tbS_filtered[,"F3"], t=5, 
                            method=methods[i], Ncensus=1000,poolSize=c(74,72))
    tb_Ne$Ne_10000[i]<-estimateNe(p0=tb_filtered[,"F0"], pt=tb_filtered[,"F3"], cov0=tbS_filtered[,"F0"], covt=tbS_filtered[,"F3"], t=5, 
                                  method=methods[i], Ncensus=10000,poolSize=c(74,72))
    
}
write.csv(tb_Ne,"Output/COV/Ne_estimation_TB91-17_angsdAF.csv")

Ne_t<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne_t$Ne01_t2[i]<-estimateNe(p0=tb_filtered[,"F0"], pt=tb_filtered[,"F1"], cov0=tbS_filtered[,"F0"], covt=tbS_filtered[,"F1"], t=1, 
                                method=methods[i], Ncensus=10000,poolSize=c(74,73))
    
    Ne_t$Ne12_t2[i]<-estimateNe(p0=tb_filtered[,"F1"], pt=tb_filtered[,"F2"], cov0=tbS_filtered[,"F1"], covt=tbS_filtered[,"F2"], t=1.7, 
                                method=methods[i], Ncensus=10000,poolSize=c(73,52))
    Ne_t$Ne23_t2[i]<-estimateNe(p0=tb_filtered[,"F2"], pt=tb_filtered[,"F3"], cov0=tbS_filtered[,"F2"], covt=tbS_filtered[,"F3"], t=1.7, 
                                method=methods[i], Ncensus=10000,poolSize=c(52,72))
    
}
write.csv(Ne_t, "Output/COV/Ne_estimation_TB_eachTimePeriod_angsdAF.csv")


## Plot Ne estimates for 3 populations

pws<-read.csv("Output/COV/Ne_estimation_PWS91-17_vcfAF.csv", row.names = 1)
#tb<- read.csv("Output/COV/Ne_estimation_SS96-17_angsdAF.csv", row.names = 1)
#ss<- read.csv("Output/COV/Ne_estimation_SS96-17_angsdAF.csv", row.names = 1)

#use PlanII (row 2, 4,6,7)
neest<-pws[c(2,4,6,7),c(1,3)]
neest$pop<-"PWS"

tb<-tb[c(2,4,6,7),c(1,3)]
tb$pop<-"TB"
ss<-ss[c(2,4,6,7),c(1,3)]
ss$pop<-"SS"

neest<-rbind(neest, tb, ss)

neest$methods[neest$methods=="W.planII"] <-"Waples"
neest$methods[neest$methods=="JR.planII"] <-"Jorde&Ryman"
neest$methods[neest$methods=="P.planII"] <-"Jónás"
neest$methods[neest$methods=="P.alt.1step.planII"] <-"Futschik"

library(colorspace) 
colors<-qualitative_hcl(5, palette = "Set2")
hcl_palettes("qualitative", plot = TRUE)

ggplot(neest, aes(x=methods, y=Ne_10000, fill=pop))+
    geom_bar(stat="identity", position=position_dodge(width = 0.7), width=0.6)+
    theme_classic()+ylab("Ne")+xlab("Method")+
    scale_fill_manual(values=colors)+theme(legend.title=element_blank())
ggsave("Output/COV/Ne_estimates_3populations.pdf", width = 6, height = 4)



## plot over years

pws<-read.csv("Output/COV/Ne_estimation_PWS_eachTimePeriod_angsdAF.csv", row.names = 1)
ss<- read.csv("Output/COV/Ne_estimation_SS_eachTimePeriod_angsdAF.csv", row.names = 1)
tb<- read.csv("Output/COV/Ne_estimation_TB_eachTimePeriod_angsdAF.csv", row.names = 1)

#use PlanII (row 2, 4,6,7)
neest<-data.frame(pop=c("PWS","SS","TB"), t1=c(pws[6,2],NA,tb[6,2]),t2=c(pws[6,3],ss[6,2],tb[6,3]),t3=c(pws[6,4],ss[6,3],tb[6,4]))

neestm<-melt(neest, id.vars="pop",value.name ="Ne")

ggplot(neestm, aes(x=variable, y=Ne, color=pop))+
    geom_point()+
    theme_classic()+ylab("Ne")+xlab("Time period")+
    geom_path(aes(x=variable, y=Ne, group=pop,color=pop))+
    scale_color_manual(values=colors)+theme(legend.title=element_blank())
ggsave("Output/COV/Ne_estimates_overtime_3populations.pdf", width = 6, height = 4)





########### 
#Estimate Ne from overlapping pop
#Ref: Waples & Yokota (2007)

# # of loci
a=nrow(pws_filtered)
#F0 to F3

#Without weighting allele freq with age structure
fc<-apply(pws_filtered[,c("F0","F3")], 1, function(x) (x["F3"]-x["F0"])^2/((x["F3"]+x["F0"])/2 -x["F3"]*x["F0"]))
Fc<-sum(fc)/a 
Fc #0.0405317
# the number of generations between the time periods
g=5    

#Ne estimates from mnon-weighted P (allele freq)
g/(2*(Fc-1/58+1/57))  
# 61.22322 (65.26726 for AF from vcf files)  


#age 3-9
pws_age<-read.csv("Data/age_structure_pws.csv")

#Weighted allele frequency
#Pw=sum(vx * Nx * Px)/ sum(vx * Nx)
# vx= reproductive value, Nx = No. of individuals of age x in the population
# Ref: Vx (Reproductive values) from Ware 1985
vx<-data.frame(age=c(3,4,5,6,7,8,9))
vx$v<-c(0.00036, 0.00047, 0.00054, 0.00059, 0.00063, 0.00065,0.00067)

#add vxNx
vx$vN0=vx$v*pws_age$F0
vx$vN1=vx$v*pws_age$F1
vx$vN2=vx$v*pws_age$F2
vx$vN3=vx$v*pws_age$F3

