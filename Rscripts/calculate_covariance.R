library(ggplot2)
library(stringr)
## Explore the difference between N100 vs. N1000 with subsampled 100 from SLiM
files<-list()
flist<-list.files("~/Projects/pacherring_Vincent/slim/", pattern="flip.csv")
for (i in 1: length(flist)){
    files[[i]]<-read.csv(paste0("~/Projects/pacherring_Vincent/slim/",flist[i]), header = F)
    names(files)[i]<-str_match(flist[i], "WF\\_(.*?)\\_flip")[2]
}


#calculate delta AF
for (f in 1: length(files)){
    df<-files[[f]]
    #calculate freq difference
    dAf<-data.frame(matrix(nrow=nrow(df), ncol=0))
    for (i in 1:9){
        dAf[,paste0("d",i)]<-df[,i]-df[,(i+1)]
    }
    write.csv(dAf, paste0("Output/COV/Sim/dAF_WF_",names(files)[f]))
    
    cov<-data.frame(matrix(nrow=0, ncol=4))
    for (i in 1:(ncol(dAf)-1)){
        for(j in (i+1):9){
            cv<-c(i, j, paste0("time",i,"-","time",j), cov(dAf[,i], dAf[,j], use="pairwise.complete.obs"))
            cov<-rbind(cov, cv)
        }
    }
    colnames(cov)<-c("start","end", "time", "cov")
    cov$start<-as.integer(cov$start)
    cov$end<-as.integer(cov$end)
    cov$cov<-as.numeric(cov$cov)
    ggplot(data=cov, aes(x=end, y=cov,group=start, color=factor(start)))+
        geom_point(size=3)+
        geom_line()+
        ylab("Covariance")+xlab('')+theme_classic()+
        theme(axis.text.x = element_blank(),legend.title = element_blank())+
        geom_hline(yintercept = 0,color="gray70", size=0.3)+ggtitle(names(files)[f])+
        scale_y_continuous(labels = scales::comma)
    ggsave(paste0("Output/COV/Sim/Cov_", names(files)[f],".png"), width = 5, height = 3, dpi=300)
    write.csv(cov, paste0("Output/COV/Sim/Covariance_", names(files)[f],".csv"))
    ggplot(data=cov, aes(x=end, y=cov,group=start, color=factor(start)))+
        geom_point(size=3)+
        geom_line()+
        ylab("Covariance")+xlab('')+theme_classic()+
        theme(axis.text.x = element_blank(),legend.title = element_blank())+
        geom_hline(yintercept = 0,color="gray70", size=0.3)+ggtitle(names(files)[f])+
        scale_y_continuous(labels = scales::comma, limits=c(-0.001, 0.001))
    ggsave(paste0("Output/COV/Sim/Cov_", names(files)[f],"_fixedYlim.png"), width = 5, height = 3, dpi=300)
    
}


## Compare with no bias correction & standardization cov output from CVTK

fname="Slim_WF_N1000_noCnoS_temp_cov_matrix_250kwin.csv"
fname="Slim_WF_N1000_sub100_temp_cov_matrix_250kwin.csv"
fname="Slim_WF_N1000_sub100_noCnoS_temp_cov_matrix_250kwin.csv"
fname="Slim_WF_N1000_sub100_noC_temp_cov_matrix_250kwin.csv"
fname="Slim_WF_N1000_sub100_wCnoS_temp_cov_matrix_250kwin.csv"
fname="Slim_WF_N100_noCwS_temp_cov_matrix_250kwin.csv"
model="WF N100 NoC wS"
model="WF N1000 sub100 wC noS"

t=10

    cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/slim/",fname))
    cov<-cov[,-1]
    #reshape the matrix
    mat1<-cov[1:(t-1),]
    mat2<-cov[t:(2*t-2),]
    
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
    ids<-str_match(covar$label, "sim:\\s(.*?)\\,\\ssim:\\s(.*?)\\)")
    ids<-ids[,2:3]
    remove<-duplicated(lapply(1:nrow(ids), function(x) {
        A<-ids[x,]
        A[order(A)]
    }  ))
    covar<-covar[!remove,] 
    
    #assign the starting time period and covering period values
    vecn<-1:(t-2)
    syr<-c()
    for (i in 1:length(vecn)){
        syr<-c(syr, rep(i, times=(t-i-1)))
    }
    covar$start_year<-syr
    
    eyr<-c()
    for (i in 1:(t-2)){
        v<-(i+1):(t-1)
        eyr<-c(eyr,v)
    }
    covar$end_year<-eyr
    
    newfile<-gsub(".csv","", fname)
    newfile<-gsub("_temp_cov_matrix","", newfile)
    ggplot(data=covar, aes(x=end_year, y=value,group=start_year, color=factor(start_year)))+
        geom_point(size=3)+
        geom_line()+
        ylab("Covariance")+xlab('')+theme_classic()+
        theme(axis.text.x = element_blank(),legend.title = element_blank())+
        geom_hline(yintercept = 0,color="gray70", size=0.3)+ggtitle(paste0(model))+
        scale_y_continuous(labels = scales::comma)
        #scale_y_continuous(labels = scales::comma, limits = c(-0.001,0.001))
    ggsave(paste0("Output/COV/Sim/", newfile, "_neutral.png"),width = 6, height = 3, dpi=300)

    
#Estimate Ne from Vars
files=list.files("~/Projects/pacherring_Vincent/slim/",pattern="Slim_WF_N10")

ne.est<-data.frame()
rows<-c()
for (f in 1: length(files)){
    cov<-read.csv(paste0("~/Projects/Pacherring_Vincent/slim/",files[f]))
    cov<-cov[,-1]
    #reshape the matrix
    mat1<-cov[1:(t-1),]
    mat2<-cov[t:(2*t-2),]
    
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
    
    
    ne=1/(2*vars$value)*10*1/6
    newfile<-gsub(".csv","", files[f])
    newfile<-gsub("_temp_cov_matrix","", newfile)
    rows<-c(rows, newfile)
    ne.est<-rbind(ne.est, ne)
}
   
rownames(ne.est)<-rows
    
#Select the ones without any corrections
ne.est2<-ne.est[c(1,4,6,10,13),]
colnames(ne.est2)<-1:9

rownames(ne.est2)<-gsub("Slim_",'',rownames(ne.est2) )
rownames(ne.est2)<-gsub("_250kwin",'',rownames(ne.est2) )

ne.est2$Mean<-rowMeans(ne.est2)

write.csv(ne.est2,"Output/COV/Sim/Ne_estimates_from_SlimOutputs_WFmodel.csv")

ne.est2[,c("Mean", 1,2,3)]

#calculate Ne from Var
    #var(p1-p0)/p0(1-p0) = 1/2N.
    # ~t/2N.(1/2N)
    # Thus, N =~ 1/(2Var)
    vars$value
    for (i in 1: nrow(vars)){
        ne=1/(2*vars$value)*10*1/6
    



ex<-read.table("Output/COV/Sim/ex.txt", header=FALSE)

x<-apply(ex, 1, function(x) str_split(x, ":"))
x1<-x[[3]]
l<-unlist(lapply(x1,as.integer))
l
