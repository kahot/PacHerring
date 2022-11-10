#import temporal covariance values generated from cvtkpy
source("Rscripts/BaseScripts.R")
pops<-c("PWS","TB","SS")
covs<-data.frame()
Variance<-data.frame()

#Use 1mb_window covariances

winsize<-c("1M","250k","100k")
winsize<-c("10k")

for (w in 1: length(winsize)){
    covs<-data.frame()
    for (p in 1: length(pops)){
        #covariance output file
        cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf05_temp_cov_matrix_",pops[p],"_",winsize[w],".csv"))
        cov<-cov[,-1]
        
        ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf05_",pops[p],"_Cov_CIs_bootstrap5000_",winsize[w],"window.csv"))
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
    ggsave(paste0("Output/COV/MD7000_Cov_overtime_CI_",winsize[w],".window_new.pdf"),width = 4.7, height = 3)
}    

## For new MD5000 vcf file
winsize<-c("100k","1M")

for (w in 1: length(winsize)){
    covs<-data.frame()
    for (p in 1: length(pops)){
        #covariance output file
        cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD5000_maf05_temp_cov_matrix_",pops[p],"_",winsize[w],".csv"))
        cov<-cov[,-1]
        
        ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD5000_maf05_",pops[p],"_Cov_CIs_bootstrap5000_",winsize[w],"window.csv"))
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
    ggsave(paste0("Output/COV/MD5000_Cov_overtime_CI_",winsize[w],".window_new.png"),width = 4.7, height = 3, dpi=100)
}    





## For new MD4000 vcf file
winsize<-c("250k","100k","1M")

for (w in 1: length(winsize)){
    covs<-data.frame()
    for (p in 1: length(pops)){
        #covariance output file
        cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD4000_maf05_temp_cov_matrix_",pops[p],"_",winsize[w],".csv"))
        cov<-cov[,-1]
        
        ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD4000_maf05_",pops[p],"_Cov_CIs_bootstrap5000_",winsize[w],"window.csv"))
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
    ggsave(paste0("Output/COV/MD4000_Cov_overtime_CI_",winsize[w],".window_new.png"),width = 4.7, height = 3, dpi=100)
}    



## MD3000 results
winsize<-c("100k","250k","1M")

for (w in 1: length(winsize)){
    covs<-data.frame()
    for (p in 1: length(pops)){
        #covariance output file
        cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD3000_maf05_temp_cov_matrix_",pops[p],"_",winsize[w],".csv"))
        cov<-cov[,-1]
        
        ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD3000_maf05_",pops[p],"_Cov_CIs_bootstrap5000_",winsize[w],"window.csv"))
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
    ggsave(paste0("Output/COV/MD3000_Cov_overtime_CI_",winsize[w],".window_new.png"),width = 4.7, height = 3, dpi=100)
}    








### using angsd AF esimates for temp cov 
covs<-data.frame()
for (p in 1: length(pops)){
    #covariance output file
    cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf05_temp_cov_matrix_",pops[p],"_1m_ANGSD_AF_try2.csv"))
    cov<-cov[,-1]
    
    ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/MD7000/MD7000_maf05_",pops[p],"_Cov_CIs_bootstrap5000_1mwindow_ANGSD_AF_try2.csv"))
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
            k=k+1}}
    colnames(covdf)<-c("label","value")
    covdf$value<-as.numeric(covdf$value)
    covar<-covdf[grep("cov",covdf$label),]
    vars<-covdf[grep("var",covdf$label),]
    
    #remove the redundant values
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
ggsave("Output/COV/MD7000_Cov_overtime_CI_1M.window_ANGSD.MF_new.pdf",width = 4.7, height = 3)





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

pw<-read.csv("~/Projects/Pacherring_Vincent/notebooks/PWS_cov12_1996-1991_2006-1996.csv")
pw[pw=="NaN"]<-NA
pw$pos<-1:nrow(pw)
plot(pw$pos, pw$nan, pch=".")    
pw1<-pw[!is.na(pw$nan),]
colnames(pw1)[1]<-"cov"

#### Variance/Ne plot
library(reshape2)
library(ggplot2)     
library(colorspace)
#cols<-qualitative_hcl(7, palette="Dark3")
cols<-sequential_hcl(8, palette = "Blue Yellow")
hcl_palettes("sequential",n=8, plot = TRUE)

vars<-list.files("Output/COV/vars/", pattern="bootstrap.csv")
labels<-list.files("Output/COV/vars/", pattern="labels.csv")
gentime<-1/6
pops<-c("PWS","SS","TB")

covfiles<-list.files("Output/COV", pattern="temp_cov_matrix")

ave_vars<-data.frame(pop=pops)
var_dfm<-data.frame()
for (i in 1:3){
    vars1000<-read.csv(paste0("Output/COV/vars/", vars[i]))
    #colnames(vars1000)<-c("1991-1996","1996-2006","2006-2017")
    varlabels<-read.csv(paste0("Output/COV/vars/", labels[i]))
    varlabels$nyears=varlabels$X1-varlabels$X2 
    
    var_df<-data.frame(T1=varlabels$nyears[1]*gentime*0.5/vars1000[,1])
    var_df$T2<-varlabels$nyears[2]*gentime*0.5/vars1000[,2]
    var_df$T3<-varlabels$nyears[3]*gentime*0.5/vars1000[,3]
    
    vardfm<-melt(var_df)
    vardfm$pop<-pops[i]
    var_dfm<-rbind(var_dfm,vardfm)
    
    df<-read.csv(paste0("Output/COV/",covfiles[i]))
    df<-df[1:3,-1]
    df<-data.frame(apply(df,1,function(x) x<-as.numeric(x)))
    ave_vars$T1[i]<-varlabels$nyears[1]*gentime*0.5/df[1,1]
    ave_vars$T2[i]<-varlabels$nyears[2]*gentime*0.5/df[2,2]
    ave_vars$T3[i]<-varlabels$nyears[3]*gentime*0.5/df[3,3]
    
}

colnames(var_dfm)[1:2]<-c("year","var")
ave_varm<-melt(ave_vars, id.vars="pop")
colnames(ave_varm)<-c("pop","year","mean")


ave_est<-data.frame(aggregate(var_dfm["var"], by=list(var_dfm$var, var_dfm$pop), mean,na.rm=T))
colnames(ave_est)<-c("year","pop","mean")


ggplot(var_dfm, aes(x=var, color=year, fill=year))+
    geom_histogram( binwidth =5, alpha=0.6)+
    facet_wrap(~pop, ncol=1)+xlim(0,2000)+xlab("Estimated Ne")+
    ylab("")+
    geom_vline(data=ave_varm, aes(xintercept=mean), color="red", size=0.5)+
    geom_vline(data=ave_est, aes(xintercept=mean, color=year), size=0.5)+
    theme_bw()+
    scale_fill_manual(values=paste0(cols[c(1,2,5)],"E6"),labels=c("1991-1996","1996-2006","2006-2017") )+
    scale_color_manual(values=paste0(cols[c(1,2,5)]),guide="none" )+
    theme(legend.title = element_blank())
ggsave("Output/COV/Ne_estimates_distribution.pdf", width = 7, height = 4)


ggplot(var_dfm, aes(x=var, color=year, fill=year))+
    geom_histogram( bins=500, alpha=0.5)+
    facet_grid(rows=vars(pop))+xlim(0,2000)+xlab("Estimated Ne")+
    ylab("")+
    geom_vline(data=ave_varm, aes(xintercept=mean), color="red", size=0.3)+
    geom_vline(data=ave_est, aes(xintercept=mean, color=year), size=0.3, type=2)+
    theme_bw()+
    scale_fill_manual(values=paste0(cols[c(1,2,5)],"E6"),labels=c("1991-1996","1996-2006","2006-2017") )+
    scale_color_manual(values=paste0(cols[c(1,2,5)]),guide="none" )+
    theme(legend.title = element_blank())
ggsave("Output/COV/Ne_estimates_distribution2.pdf", width = 7, height = 4)



### Interpopulation comparisons
#decode the samples to create the right matrix
cv<-read.csv("~/Projects/Pacherring_Vincent/notebooks/gw_covs_interPopulations_MD7000_100k.csv", header = F)
labs<-read.csv("~/Projects/Pacherring_Vincent/notebooks/gw_covs_interPopulations_MD7000_100k_labels.csv")
labs<-labs[,-1]
cvm<-data.frame(label=as.vector(t(labs)), cov=as.vector(t(cv)))

#rearrange based on comparions: covariance between populations within each time period
#1-2 vs. 4-5 vs. 7-8 (1996-2006)
#2-3 vs. 5-6 vs. 8-9 (2006-2017)

Covs<-data.frame(pops=rep(c("PWS.vs.SS", "PWS.vs.TB",  "SS.vs.TB"), times=2),
                 period=c(rep("1996-2006", times=3), rep("2006-2017", times=3)))

Covs$cov<-c(cvm$cov[cvm$label=="cov(PH: 2-1, PH: 5-4)"],cvm$cov[cvm$label=="cov(PH: 2-1, PH: 8-7)"],cvm$cov[cvm$label=="cov(PH: 5-4, PH: 8-7)"], 
            cvm$cov[cvm$label=="cov(PH: 6-5, PH: 3-2)"],cvm$cov[cvm$label=="cov(PH: 9-8, PH: 3-2)"],cvm$cov[cvm$label=="cov(PH: 9-8, PH: 6-5)"])
#C.I.
cis<-read.csv("~/Projects/Pacherring_Vincent/notebooks/Interpop_comparison_MD7000_CIs.csv")
cis<-cis[,-1]
cim<-data.frame(label=as.vector(t(labs)), ci_l=as.vector(t(cis[1:8,])))
cim$ci_h<-as.vector(t(cis[9:16,]))

Covs$ci_l<-as.numeric(c(cim$ci_l[cim$label=="cov(PH: 2-1, PH: 5-4)"],cim$ci_l[cim$label=="cov(PH: 2-1, PH: 8-7)"],cim$ci_l[cim$label=="cov(PH: 5-4, PH: 8-7)"], 
             cim$ci_l[cim$label=="cov(PH: 6-5, PH: 3-2)"],cim$ci_l[cim$label=="cov(PH: 9-8, PH: 3-2)"],cim$ci_l[cim$label=="cov(PH: 9-8, PH: 6-5)"]))

Covs$ci_h<-as.numeric(c(cim$ci_h[cim$label=="cov(PH: 2-1, PH: 5-4)"],cim$ci_h[cim$label=="cov(PH: 2-1, PH: 8-7)"],cim$ci_h[cim$label=="cov(PH: 5-4, PH: 8-7)"], 
             cim$ci_h[cim$label=="cov(PH: 6-5, PH: 3-2)"],cim$ci_h[cim$label=="cov(PH: 9-8, PH: 3-2)"],cim$ci_h[cim$label=="cov(PH: 9-8, PH: 6-5)"]))

ggplot(Covs, aes(x=period, y=cov, fill=pops))+
    geom_bar(stat="identity",position=position_dodge(width = 0.7), width=0.8)+
    ylab("Covariance")+xlab('')+theme_classic()+
    geom_hline(yintercept = 0,color="gray70", size=0.3)+
    scale_fill_manual(values=c(blu,grb,red))+
    theme(legend.title = element_blank())+
    geom_errorbar(aes(ymin=ci_l, ymax=ci_h), width=.2, size=.2, position=position_dodge(width = 0.7))
ggsave("Output/COV/Interpop_cov_comparison_MD7000.pdf",width = 6, height = 4)



## MD2000 plot
covs<-read.csv("Output/Cov/Cov_interpop.csv")
vars<-covs[grep("var", covs$label),]
covs<-covs[grep("cov", covs$label),]
covs$label<-gsub("cov\\(", "", covs$label)
covs$label<-gsub("\\)", "", covs$label)
write.csv(covs,"Output/COV/Covarince_interpop.csv")

#replace the target samples
covs<-read.csv("Output/COV/interpop_cov_extracted.csv")


bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'
ggplot(covs, aes(x=period, y=cov, fill=pops))+
    geom_bar(stat="identity",position=position_dodge(width = 0.7), width=0.8)+
    ylab("Covariance")+xlab('')+theme_classic()+
    geom_hline(yintercept = 0,color="gray70", size=0.3)+
    scale_fill_manual(values=c(blu,grb,red))+
    theme(legend.title = element_blank())+
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2, size=.2, position=position_dodge(width = 0.7))
ggsave("Output/COV/Interpop_cov_comparison.pdf",width = 6, height = 4)



# Plot G(t) for PWS and TB

gt<-data.frame(year=c(0,1,2,0,1,2), pop=c("PWS","PWS","PWS","TB","TB","TB"), 
               g=c(0, 0.00137183,0.00183995 ,0, 0.0012538,0.00147975))


gt2<-data.frame(year=c(0,1,2,0,1,2), pop=c("PWS","PWS","PWS","TB","TB","TB"), 
               g=c(0, 0.00137183,0.00183995 ,0, 0.0012538,0.00147975))


library(scales)
show_col(hue_pal()(3))

ggplot(data=gt, aes(x=year, y=g*100, color=pop))+
    geom_point(size=3, position=position_dodge(width = 0.1,preserve ="total"))+
    geom_line(position=position_dodge(width = 0.1,preserve ="total"))+
    ylab("G(t) %")+xlab('')+theme_classic()+
    theme(axis.text.x = element_blank(),legend.title = element_blank())+
    scale_x_continuous(breaks = c(0,1,2))+
    scale_color_manual(values=c("#F8766D", "#619CFF"))
ggsave("Output/COV/MD7000_Gt_overtime.pdf",width = 4.7, height = 2.1)

# Plot G(t) for all pop for no91 data 

gt<-data.frame(year=c(0,1,2,0,1,2,0,1,2), pop=c("PWS","PWS","PWS","SS","SS","SS","TB","TB","TB"), 
               g=c(0, 0.00368883,0.00419642,0,0.0032369,0.00447537, 0, 0.00292324,0.0034146))

library(scales)
show_col(hue_pal()(3))

ggplot(data=gt, aes(x=year, y=g*100, color=pop))+
    geom_point(size=3, position=position_dodge(width = 0.1,preserve ="total"))+
    geom_line(position=position_dodge(width = 0.1,preserve ="total"))+
    ylab("G(t) %")+xlab('')+theme_classic()+
    theme(axis.text.x = element_blank(),legend.title = element_blank())+
    scale_x_continuous(breaks = c(0,1,2))
    scale_color_manual(values=c("#F8766D", "#619CFF"))
ggsave("Output/COV/MD7000_Gt_overtime.pdf",width = 4.7, height = 2.1)





### Read tiled temporal covariance data PWS (100k)
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_100kwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_100kwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_100kwindow.csv", header = F)
iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_intervals_100Kwindow.csv", row.names = 1)

#1M window
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_1mwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_1mwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_1mwindow.csv", header = F)
iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_intervals_1mwindow.csv", row.names = 1)

#250k window
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_250kwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_250kwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_250kwindow.csv", header = F)
iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_intervals_250kwindow.csv", row.names = 1)



covs<-cbind(iv, cov12, cov23,cov13)
colnames(covs)[4:6]<-c("cov12","cov23","cov13")
covs$index=1:nrow(covs)

evens<-paste0("chr",seq(2,26, by=2))
covs$color<-"col1"
covs$color[covs$chrom %in% evens]<-"col2"

covs$cov12[is.nan(covs$cov12)]<-NA
covs$cov12[is.infinite(covs$cov12)]<-NA
write.csv(covs,"Output/COV/PWS_tempCovs_md7000_100k.csv")
write.csv(covs,"Output/COV/PWS_tempCovs_md7000_1m.csv")


ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=1, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())+
    ggtitle("PWS cov12")
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_1mWindow.pdf", width = 6, height = 3)    
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_100kWindow.pdf", width = 6, height = 3)    
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_250kWindow.pdf", width = 6, height = 3)    

ggplot(covs, aes(x=index, y=cov23, color=color))+
    geom_point(size=1, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())+
    ggtitle("PWS cov23")
ggsave("Output/COV/PWS_tempCovs23_acrossGenome_1mWindow.pdf", width = 6, height = 3)    
ggsave("Output/COV/PWS_tempCovs23_acrossGenome_250kWindow.pdf", width = 6, height = 3)    



# Where are the largest cov regions?

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



## Run snpEff
# Create a new vcf file containing the loci in p07_loc.
## at terminal
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_tempcov12_1M_outliers.bed --out Output/COV/annotation/PWS_cov12_outlier_1m --recode --keep-INFO-all
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_tempcov23_1M_outliers.bed --out Output/COV/annotation/PWS_cov23_outlier_1m --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12_outlier_1m.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12 > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov12_outlier_1m.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23_outlier_1m.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23 > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov23_outlier_1m.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov12_outlier_1m.vcf > Output/COV/annotation/PWS_cov12_1m_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov23_outlier_1m.vcf > Output/COV/annotation/PWS_cov23_1m_annotation


#snpEff results
compa<-c("cov12","cov23")

for (i in 1:2){
    df<-read.table(paste0("Output/COV/annotation/PWS_",compa[i],"_1m_annotation"), header = F)
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
    write.csv(Ano2, paste0("Output/COV/annotation/PWS_",compa[i],"_1m_outlier_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids1<-geneids[nchar(geneids)<=18]
    geneids1<-unique(geneids1)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_1m_outlier_geneid_list1.txt"))
    cat(paste0(geneids1,"; "))
    sink(NULL)
    
    #split the intergenic ids into two
    geneids2<-geneids[nchar(geneids)>18]
    ids2<-unlist(str_split(geneids2, "-", 2))
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-c(geneids, ids2)
    geneids<-unique(geneids)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_1m_outlier_geneid_list2.txt"))
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
    
    write.table(genenames2, paste0("Output/COV/annotation/PWS_",compa[i],"_1m_outlier_genenames_list_test.txt"), quote=F, row.names = F, col.names = F)
}


##### MD7000 10k windows
### Read tiled temporal covariance data PWS (100k)
cov12<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov12_1996-1991_2006-1996_md7000_10kwindow.csv", header = F)
cov23<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov23_2017-2006_2006-1996_md7000_10kwindow.csv", header = F)
cov13<-read.csv("~/Projects/Pacherring_Vincent/MD7000/PWS_cov13_2017-2006_1996-1991_md7000_10kwindow.csv", header = F)

iv<-read.csv("~/Projects/Pacherring_Vincent/MD7000/MD7000_intervals_10Kwindow.csv", row.names = 1)

covs<-cbind(iv, cov12, cov23,cov13)
colnames(covs)[4:6]<-c("cov12","cov23","cov13")

covs$ch<-gsub("chr","", covs$chrom)
covs$ch<-as.integer(covs$ch)
covs<-covs[order(covs$ch),]

covs$index=1:nrow(covs)

evens<-paste0("chr",seq(2,26, by=2))
covs$color<-"col1"
covs$color[covs$chrom %in% evens]<-"col2"

covs$cov12[is.nan(covs$cov12)]<-NA
covs$cov12[is.infinite(covs$cov12)]<-NA
covs$cov23[is.nan(covs$cov23)]<-NA
covs$cov23[is.infinite(covs$cov23)]<-NA
covs$cov13[is.nan(covs$cov13)]<-NA
covs$cov13[is.infinite(covs$cov13)]<-NA

write.csv(covs,"Output/COV/PWS_tempCovs_md7000_10k.csv")

ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_10k.pdf", width = 6, height = 3)    

ggplot(covs, aes(x=index, y=cov12, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+ylim(-0.2,0.2)+xlab('Chromosome')+
    theme(axis.text.x = element_blank())
ggsave("Output/COV/PWS_tempCovs12_acrossGenome_10k_zoomed.pdf", width = 6, height = 3)    


# Where are the largest cov regions?
covs2<-covs[order(covs$cov12, decreasing=T),]


covs2_1<-covs2[covs2$cov12>=0.1,]
covs2_1<-covs2_1[!is.na(covs2_1$cov12),] #226 windows above 
#top 1%
covs2_top<-covs2[1:705,c("chrom","start","end","cov12","ch")]
covs2_top<-covs2_top[order(covs2_top$ch, covs2_top$start),]
#create a bed file to find genes in these regions
write.table(format(covs2_top[,1:3], scientific=FALSE), "Output/COV/pws_cov12_outliers_10k.bed", quote = F, row.names = F, col.names = F,sep = "\t")


# 2nd time period
ggplot(covs, aes(x=index, y=cov23, color=color))+
    geom_point(size=.5, alpha=0.5)+
    theme_classic()+
    scale_color_manual(values=c("gray70","steelblue"), guide="none")+
    ylab("Covariance")+xlab('Chromosome')+
    theme(axis.text.x = element_blank())
ggsave("Output/COV/PWS_tempCovs23_acrossGenome_10k.pdf", width = 6, height = 3)    


covs3<-covs[order(covs$cov23, decreasing=T),]
#top 1%
covs3_top<-covs3[1:705,c("chrom","start","end","cov23","ch")]

covs3_top<-covs3_top[order(covs3_top$ch, covs3_top$start),]
write.table(format(covs3_top[,1:3], scientific=FALSE), "Output/COV/pws_cov23_outliers_10k.bed", quote = F, row.names = F, col.names = F,sep = "\t")

## Run snpEff
# Create a new vcf file containing the loci in p07_loc.
## at terminal
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_cov12_outliers_10k.bed --out Output/COV/annotation/PWS_cov12_outlier_10k --recode --keep-INFO-all
vcftools --gzvcf Data/new_vcf/population/PWS07_maf05.vcf.gz --bed Output/COV/pws_cov23_outliers_10k.bed --out Output/COV/annotation/PWS_cov23_outlier_10k --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12_outlier_10k.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov12_10k > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov12_outlier_10k.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23_outlier_10k.recode.vcf -stats ~/Projects/PacHerring/Output/COV/annotation/PWS_cov23_10k > ~/Projects/PacHerring/Output/COV/annotation/Anno.PWS_cov23_outlier_10k.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov12_outlier_10k.vcf > Output/COV/annotation/PWS_cov12_10k_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/COV/annotation/Anno.PWS_cov23_outlier_10k.vcf > Output/COV/annotation/PWS_cov23_10k_annotation


#snpEff results
compa<-c("cov12","cov23")

for (i in 1:2){
    df<-read.table(paste0("Output/COV/annotation/PWS_",compa[i],"_10k_annotation"), header = F)
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
    write.csv(Ano2, paste0("Output/COV/annotation/PWS_",compa[i],"_10k_outlier_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids1<-geneids[nchar(geneids)<=18]
    geneids1<-unique(geneids1)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_10k_outlier_geneid_list1.txt"))
    cat(paste0(geneids1,"; "))
    sink(NULL)
    
    #split the intergenic ids into two
    geneids2<-geneids[nchar(geneids)>18]
    ids2<-unlist(str_split(geneids2, "-", 2))
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-c(geneids, ids2)
    geneids<-unique(geneids)
    sink(paste0("Output/COV/annotation/PWS_",compa[i],"_10k_outlier_geneid_list2.txt"))
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
    
    write.table(genenames2, paste0("Output/COV/annotation/PWS_",compa[i],"_10k_outlier_genenames_list.txt"), quote=F, row.names = F, col.names = F)
}


#Background genes 
df<-read.table(paste0("Data/new_vcf/md7000_maf01_annotation"), header = F)

df3<-df %>%
    separate(V4, paste0("c",1:5), "\\|")
colnames(df3)[4:8]<-c("Allele","Annotation","Putative_impact","Gene_name", "Gene_ID")
remove<-!duplicated(df3[4:8])
anno<-df3[remove,] #274403
write.csv(anno, paste0("~/Projects/PacHerring/Data/new_vcf/md7000_maf01_annotation.csv"))

#anno<-read.csv("Data/new_vcf/md7000_maf01_annotation.csv", row.names = 1)

geneids<-c(anno$Gene_ID)
geneids<-unique(geneids) #33338
#split the intergenic ids into two
geneids2<-geneids[nchar(geneids)>18]
ids2<-unlist(str_split(geneids2, "-", 2))
geneids<-geneids[nchar(geneids)<=18]
geneids<-c(geneids, ids2)
geneids<-unique(geneids) #26982
write.table(geneids, "Data/new_vcf//md7000_maf01_background_geneids.xt", quote=F, row.names = F, col.names = F)



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

write.table(genenames2, paste0("Output/COV/annotation/PWS_",compa[i],"_10k_outlier_genenames_list.txt"), quote=F, row.names = F, col.names = F)
