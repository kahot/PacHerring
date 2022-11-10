#Estimate effective population size
#Nest/PoolSeq package is used. Ref: doi: 10.1534/genetics.116.191197
library(poolSeq)

#Read the allele freq data for PWS

pws1<-read.table("Data/new_vcf/population/PWS91_freq.frq",skip=1, stringsAsFactors = FALSE, sep="\t",col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws2<-read.table("Data/new_vcf/population/PWS96_freq.frq",skip=1, stringsAsFactors = FALSE, sep="\t",col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws3<-read.table("Data/new_vcf/population/PWS07_freq.frq",skip=1, stringsAsFactors = FALSE, sep="\t",col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
pws4<-read.table("Data/new_vcf/population/PWS17_freq.frq",skip=1, stringsAsFactors = FALSE, sep="\t",col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))


for (i in c(1:4)){
    df<-get(paste0("pws",i))
    df$maf<-substr(df$MAF, 3,10)
    df<-df[,c(1,2,7)]
    df$maf<-as.numeric(df$maf)
    assign(paste0("pws.t",i),df)
}

#AF tables for all years
pws<-cbind(pws.t1, pws.t2[,3],pws.t3[,3],pws.t4[,3])
colnames(pws)[3:6]<-c("F0","F1","F2","F3")

#Read depth information
p_reads<-data.frame()
for (i in 1:26){
    df<-read.csv(paste0("Output/CNV/chr",i,"_depth.csv"), row.names = 1)
    df_p91<-df[,grep("PWS91", colnames(df))]
    df_p96<-df[,grep("PWS96", colnames(df))]
    df_p07<-df[,grep("PWS07", colnames(df))]
    df_p17<-df[,grep("PWS17", colnames(df))]
    reads<-data.frame(F0=rowSums(df_p91, na.rm = T), F1=rowSums(df_p96, na.rm = T),F2=rowSums(df_p07, na.rm = T),F3=rowSums(df_p17, na.rm = T))
    p_reads<-rbind(p_reads,reads)
}


#Find SNPs with extreme values and uninformative loci and remove them
remove<-checkSNP(pws[,"F0"],pws[,"F3"],p_reads[,"F0"], p_reads[,"F3"])
#filtered the snp dataset
pws_filtered<-pws[remove,]
pwsS_filtered<-p_reads[remove,]

#Look at F0 and F1
remove<-checkSNP(pws_filtered[,"F0"],pws_filtered[,"F1"],pwsS_filtered[,"F0"], pwsS_filtered[,"F1"])
length(remove[remove==F]) #0

remove<-checkSNP(pws_filtered[,"F0"],pws_filtered[,"F2"],pwsS_filtered[,"F0"], pwsS_filtered[,"F2"])
length(remove[remove==F]) #0

remove<-checkSNP(pws_filtered[,"F1"],pws_filtered[,"F2"],pwsS_filtered[,"F1"], pwsS_filtered[,"F2"])
length(remove[remove==F]) #5376
pws_filtered<-pws_filtered[remove,]
pwsS_filtered<-pwsS_filtered[remove,]

remove<-checkSNP(pws_filtered[,"F2"],pws_filtered[,"F3"],pwsS_filtered[,"F2"], pwsS_filtered[,"F3"])
length(remove[remove==F])
pws_filtered<-pws_filtered[remove,]
pwsS_filtered<-pwsS_filtered[remove,]

write.csv(pws_filtered, "Output/COV/PWS_snps_filtered_forNeestimate.csv")



#estimate Ne using 1mb window size using W.planI (Waples 1989)
pws.ne1<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=1000000, p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
              unit="bp", method="W.planI", Ncensus=1000, poolSize=c(58,56))
mean(pws.ne1$Nw.planI) #121.4605

#Change the window size to 10M
pws.ne10m<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=10000000, p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                       unit="bp", method="W.planI", Ncensus=1000, poolSize=c(58,56))
mean(pws.ne10m$Nw.planI) #117.2566


#Change the window size to 100K
pws.ne100k<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=100000, p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                       unit="bp", method="W.planI", Ncensus=1000, poolSize=c(58,56))
mean(pws.ne100k$Nw.planI) #126.7154

#planII method
pws.ne11<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=1000000, p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                       unit="bp", method="W.planII", poolSize=c(58,56))
mean(pws.ne11$Nw.planII)#128.2156

est<-unlist(as.vector(pws.ne1[,1]))
mean(est) #-229.5022

#between 1991-1996
pws.ne01<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=1000000, p0=pws_filtered[,"F0"], pt=pws_filtered[,"F1"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F1"], t=1, 
                       unit="bp", method="W.planII", poolSize=c(58,72))
mean(pws.ne01$Nw.planII) #29.50389
pws.ne12<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=1000000, p0=pws_filtered[,"F1"], pt=pws_filtered[,"F2"], cov0=pwsS_filtered[,"F1"], covt=pwsS_filtered[,"F2"], t=1.8, 
                        unit="bp", method="W.planI", Ncensus=72, poolSize=c(71,46))
mean(pws.ne12$Nw.planI) #216.7715

pws.ne23<-estimateWndNe(chr=pws_filtered$chr, pos=pws_filtered$pos, wndSize=1000000, p0=pws_filtered[,"F2"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F2"], covt=pwsS_filtered[,"F3"], t=1.67, 
                        unit="bp", method="W.planI", Ncensus=46, poolSize=c(46,56))
mean(pws.ne23$Nw.planI) #103.5082

#### 

#Estimate genome-wide Ne
methods<-c("W.planI","W.planII","JR.planI","JR.planII","P.planI","P.planII","P.alt.1step.planII")

pws_Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    pws_Ne$Ne[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                  method=methods[i], Ncensus=1000,poolSize=c(58,56))
    pws_Ne$Ne_10000[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F3"], t=5, 
                                   method=methods[i], Ncensus=10000,poolSize=c(58,56))
    
}
write.csv(pws_Ne,"Output/COV/Ne_estimation_PWS91-17.csv")

Ne<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne$Ne01_t1[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F1"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F1"], t=1, 
                         method=methods[i], Ncensus=10000,poolSize=c(58,72))
    Ne$Ne01_t5[i]<-estimateNe(p0=pws_filtered[,"F0"], pt=pws_filtered[,"F1"], cov0=pwsS_filtered[,"F0"], covt=pwsS_filtered[,"F1"], t=5, 
                                  method=methods[i], Ncensus=10000,poolSize=c(58,72))
    Ne$Ne12_t2[i]<-estimateNe(p0=pws_filtered[,"F1"], pt=pws_filtered[,"F2"], cov0=pwsS_filtered[,"F1"], covt=pwsS_filtered[,"F2"], t=1.8, 
                              method=methods[i], Ncensus=10000,poolSize=c(72,46))
    Ne$Ne12_t11[i]<-estimateNe(p0=pws_filtered[,"F1"], pt=pws_filtered[,"F2"], cov0=pwsS_filtered[,"F1"], covt=pwsS_filtered[,"F2"], t=11, 
                              method=methods[i], Ncensus=10000,poolSize=c(72,46))
    Ne$Ne23_t2[i]<-estimateNe(p0=pws_filtered[,"F2"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F2"], covt=pwsS_filtered[,"F3"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(72,46))
    Ne$Ne23_t10[i]<-estimateNe(p0=pws_filtered[,"F2"], pt=pws_filtered[,"F3"], cov0=pwsS_filtered[,"F2"], covt=pwsS_filtered[,"F3"], t=10, 
                               method=methods[i], Ncensus=10000,poolSize=c(72,46))
    
}


# Stika Sound
S1 <- read.table("Data/freqs/MD7000/SS96_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
S2 <- read.table("Data/freqs/MD7000/SS06_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
S3 <- read.table("Data/freqs/MD7000/SS17_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))

for (i in c(1:3)){
    df<-get(paste0("S",i))
    df$maf<-substr(df$MAF, 3,10)
    df<-df[,c(1,2,7)]
    df$maf<-as.numeric(df$maf)
    assign(paste0("ss.t",i),df)
}

#AF tables for all years
ss<-cbind(ss.t1, ss.t2[,3],ss.t3[,3])
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
T1 <- read.table("Data/freqs/MD7000/TB91_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
T2 <- read.table("Data/freqs/MD7000/TB96_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
T3 <- read.table("Data/freqs/MD7000/TB06_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
T4 <- read.table("Data/freqs/MD7000/TB17_md7000_maf05_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))


for (i in c(1:4)){
    df<-get(paste0("T",i))
    df$maf<-substr(df$MAF, 3,10)
    df<-df[,c(1,2,7)]
    df$maf<-as.numeric(df$maf)
    assign(paste0("tb.t",i),df)
}

#AF tables for all years
tb<-cbind(tb.t1, tb.t2[,3],tb.t3[,3],tb.t4[,3])
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
write.csv(ss_Ne,"Output/COV/Ne_estimation_SS96-17.csv")

Ne_s<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne_s$Ne12_t2[i]<-estimateNe(p0=ss_filtered[,"F1"], pt=ss_filtered[,"F2"], cov0=ssS_filtered[,"F1"], covt=ssS_filtered[,"F2"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(78,64))
    Ne_s$Ne23_t2[i]<-estimateNe(p0=ss_filtered[,"F2"], pt=ss_filtered[,"F3"], cov0=ssS_filtered[,"F2"], covt=ssS_filtered[,"F3"], t=1.7, 
                              method=methods[i], Ncensus=10000,poolSize=c(78,64))
    
}
write.csv(Ne_s, "Output/COV/Ne_estimation_SS_eachTimePeriod.csv")


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
write.csv(tb_Ne,"Output/COV/Ne_estimation_TB91-17.csv")

Ne_t<-data.frame(methods=methods)
for (i in 1: length(methods)){
    Ne_t$Ne01_t2[i]<-estimateNe(p0=tb_filtered[,"F0"], pt=tb_filtered[,"F1"], cov0=tbS_filtered[,"F0"], covt=tbS_filtered[,"F1"], t=1, 
                                method=methods[i], Ncensus=10000,poolSize=c(74,73))
    
    Ne_t$Ne12_t2[i]<-estimateNe(p0=tb_filtered[,"F1"], pt=tb_filtered[,"F2"], cov0=tbS_filtered[,"F1"], covt=tbS_filtered[,"F2"], t=1.7, 
                                method=methods[i], Ncensus=10000,poolSize=c(73,52))
    Ne_t$Ne23_t2[i]<-estimateNe(p0=tb_filtered[,"F2"], pt=tb_filtered[,"F3"], cov0=tbS_filtered[,"F2"], covt=tbS_filtered[,"F3"], t=1.7, 
                                method=methods[i], Ncensus=10000,poolSize=c(52,72))
    
}
write.csv(Ne_t, "Output/COV/Ne_estimation_TB_eachTimePeriod.csv")



########### 
#Estimate Ne from overlapping pop
#Ref: Waples & Yokota (2007)

# # of loci
a=nrow(pws_filtered)
#F0 to F3

#Without weighting allele freq with age structure
fc<-apply(pws_filtered[,c("F0","F3")], 1, function(x) (x["F3"]-x["F0"])^2/((x["F3"]+x["F0"])/2 -x["F3"]*x["F0"]))
Fc<-sum(fc)/a #0.03800157

# the number of generations between the time periods
g=5    

#Ne estimates from mnon-weighted P (allele freq)
g/(2*(Fc-1/58+1/57))  
# 65.26726  


#age 3-9
pws_age<-read.csv("Data/age_structure_pws.csv")

#Weighted allele frequency
# Pw=sum(vx * Nx * Px)/ sum(vx * Nx)
# vx= reproductive value, Nx = No. of individuals of age x in the population
# Ref: Vx (Reproductive values) from Ware 1985

vx<-data.frame(age=c(3,4,5,6,7,8,9))
#reproductive values (v)
vx$v<-c(0.00036, 0.00047, 0.00054, 0.00059, 0.00063, 0.00065,0.00067)

#add vxNx
vx$vN0=vx$v*pws_age$F0
vx$vN1=vx$v*pws_age$F1
vx$vN2=vx$v*pws_age$F2
vx$vN3=vx$v*pws_age$F3


### calculate weighted P ###

#Need to calculate allele freq for each age class 
G91<-read.csv("Output/PWS91_genotypes.csv",row.names = 1)
G96<-read.csv("Output/PWS96_genotypes.csv",row.names = 1)
G07<-read.csv("Output/PWS07_genotypes.csv",row.names = 1)
G17<-read.csv("Output/PWS17_genotypes.csv",row.names = 1)


#randomly assign the individuals to the age bin 

# subset eaech time point into 7 age groups 
# If only 1 individual, skip 

#Calculate allele freq for each age group
Gs<-c("G91","G96","G07","G17")
AF_list<-list()
for (j in 1:4){
    df<-get(Gs[j])
    AF<-data.frame(matrix(nrow=nrow(df)))
    for (i in 1:7){
        if (pws_age[i,(j+1)]==0|pws_age[i,(j+1)]==1) {
            AF[,i]<-NA
            colnames(AF)[i]<-paste0("age", i+2)}
        else{
            a<-sample(1:ncol(df), pws_age[i,(j+1)], replace=F)
            agedf<-data.frame(df[, a])
            df<-data.frame(df[,-a])
            age_af<-apply(agedf, 1, function(x) {fq<-x[x!=9]
                    af<-sum(fq)/(length(fq)*2)
                    return(af)})
            AF[,i]<-age_af
            colnames(AF)[i]<-paste0("age", i+2)
        }    
    }
    AF_list[[j]]<-AF
    names(AF_list)[j]<-Gs[j]
}
        
#Weighted freq (P) calculation for each locus    
#Pw=sum(vx * Nx * Px)/ sum(vx * Nx)

Pw<-data.frame(pos=rownames(G91))
for (j in 1:4){
    df<-AF_list[[j]]
    n<-j+2
    pw<-apply(df, 1, function(x) sum(x*vx[,n], na.rm=T)/sum(vx[,n],na.rm=T))
    Pw[(j+1)]<-pw
    colnames(Pw)[j+1]<-paste0("F",j-1)
}

write.csv(Pw, "Output/COV/Weighted_allele_freq_PWS.csv")

#run the checkSNP

p_reads2<-p_reads
p_reads2$pos<-rownames(p_reads2)
p_reads2<-p_reads2[p_reads2$pos %in% Pw$pos,]

remove<-checkSNP(Pw[,"F0"],Pw[,"F3"],p_reads2[,"F0"], p_reads2[,"F3"])
length(remove[remove==F]) #10296
Pw_filtered<-Pw[remove,]


#calculate Fc
fc<-apply(Pw_filtered[,c("F0","F3")], 1, function(x) (x["F3"]-x["F0"])^2/((x["F3"]+x["F0"])/2 -x["F3"]*x["F0"]))
Fc<-sum(fc)/nrow(Pw_filtered)

# the number of generations between the time periods
g=5    

#Ne estimates from weighted Pw
g/(2*(Fc-1/58+1/57))  
# 1: 53.90429   
# 2: 53.90449
# 3: 62.95774


#run 100 times
Gs<-c("G91","G96","G07","G17")
# the number of generations between the time periods
g=5    

p_reads2<-p_reads
p_reads2$pos<-rownames(p_reads2)
p_reads2<-p_reads2[p_reads2$pos %in% Pw$pos,]

age_ne<-data.frame(iter=1:100)
for (k in 1:100){
    
    AF_list<-list()
    m=1
    for (j in c(1,4)){
        df<-get(Gs[j])
        AF<-data.frame(matrix(nrow=nrow(df)))
        for (i in 1:7){
            if (pws_age[i,(j+1)]==0|pws_age[i,(j+1)]==1) {
                AF[,i]<-NA
                colnames(AF)[i]<-paste0("age", i+2)}
            else{
                a<-sample(1:ncol(df), pws_age[i,(j+1)], replace=F)
                agedf<-data.frame(df[, a])
                df<-data.frame(df[,-a])
                age_af<-apply(agedf, 1, function(x) {fq<-x[x!=9]
                af<-sum(fq)/(length(fq)*2)
                return(af)})
                AF[,i]<-age_af
                colnames(AF)[i]<-paste0("age", i+2)
            }    
        }
        
        AF_list[[m]]<-AF
        names(AF_list)[m]<-Gs[j]
        m=m+1
    }
    
    #Weighted freq (P) calculation for each locus    
    #Pw=sum(vx * Nx * Px)/ sum(vx * Nx)
    
    Pw<-data.frame(pos=rownames(G91))
    for (j in 1:2){
        df<-AF_list[[j]]
        n<-j+2
        pw<-apply(df, 1, function(x) sum(x*vx[,n], na.rm=T)/sum(vx[,n],na.rm=T))
        Pw[(j+1)]<-pw
        if (j==1) cname="F0"
        if (j==2) cname="F3"
        colnames(Pw)[j+1]<-cname
    }
    
  
    remove<-checkSNP(Pw[,"F0"],Pw[,"F3"],p_reads2[,"F0"], p_reads2[,"F3"])
    Pw_filtered<-Pw[remove,]
    
    #calculate Fc
    fc<-apply(Pw_filtered[,c("F0","F3")], 1, function(x) (x["F3"]-x["F0"])^2/((x["F3"]+x["F0"])/2 -x["F3"]*x["F0"]))
    Fc<-sum(fc)/nrow(Pw_filtered)
    
    #Ne estimates from weighted Pw
    age_ne$ne[k]<- g/(2*(Fc-1/58+1/57)) 
    
}

