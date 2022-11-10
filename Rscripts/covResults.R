covResults<-function(fname, model, t){
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
        scale_y_continuous(labels = scales::comma, limits = c(-0.001,0.001))
    ggsave(paste0("../Output/COV/Sim/", newfile, "_neutral.png"),width = 6, height = 3, dpi=300)
    print(paste0("../Output/COV/Sim/", newfile, "_neutral.png"))
}


#No ylim
covResults2<-function(fname, model, t){
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
    ggsave(paste0("../Output/COV/Sim/", newfile, "_neutral.png"),width = 6, height = 3, dpi=300)
    print(paste0("../Output/COV/Sim/", newfile, "_neutral.png"))
}