#import temporal covariance values generated from cvtkpy
library(ggplot2)

pops<-c("PWS","TB","SS")
covs<-data.frame()
Variance<-data.frame()

#Use 1mb_window covariances

for (p in 1: length(pops)){
    #covariance output file
    cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/notebooks/temp_cov_matrix_",pops[p],"_1mwindow.csv"))
    cov<-cov[,-1]
    
    ci<-read.csv(paste0("~/Projects/Pacherring_Vincent/notebooks/",pops[p],"_Cov_CIs_bootstrap5000_1M.csv"))
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
    scale_shape_manual(values=c(16,17),labels=c("1991-","1996-"))
ggsave("Output/COV/Cov_overtime_CI_1Mwindow.pdf",width = 4.7, height = 3)



ggplot(covs, aes(x=year, y=value, color=series))+
    geom_point()+
    geom_path(aes(color=series))+
    xlim(1,3)+ylab("Covariance")+xlab('')+theme_classic()+ylim(-0.0035,0.001)+
    geom_hline(yintercept = 0,color="gray70", size=0.3)+
    theme(axis.text.x = element_blank(),legend.title = element_blank())+
    ggtitle("TB")
ggsave("Output/COV/TB_cov_overtime.pdf",width = 4.7, height = 3)



library(ggplot2)

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
