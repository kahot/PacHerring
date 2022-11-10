library(edgeR)

pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws96<-pops_info$Sample[pops_info$Population.Year=="PWS96"]
pws91<-pops_info$Sample[pops_info$Population.Year=="PWS91"]
pws07<-pops_info$Sample[pops_info$Population.Year=="PWS07"]
ca<-pops_info$Sample[pops_info$Population.Year=="CA17"]

pops<-c("pws91","pws96","pws07","ca")


#x <- read.delim("TableOfCounts.txt",row.names="Symbol")
# group <- factor(c(1,1,2,2))
# y <- DGEList(counts=x,group=group)
# keep <- filterByExpr(y)
# y <- y[keep,,keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# design <- model.matrix(~group)
# y <- estimateDisp(y,design)

#df<-read.csv(paste0("data/bam_depth/",plist[i],"_chr6-cyp.txt"), header=F, sep="\t")

#for (p in 1:length(pops)){
#    plist<-get(paste0(pops[p]))
#    plots<-list()
    depth1<-data.frame(pos=2000001:3500000)
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",pws91[i],"_chr6-cyp.txt"), header=F, sep="\t")
        depth1<-merge(depth1, df[,2:3], by.x="pos",by.y="V2")
    }        
    
    depth2<-data.frame(pos=2000001:3500000)
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",pws07[i],"_chr6-cyp.txt"), header=F, sep="\t")
        depth2<-merge(depth2, df[,2:3], by.x="pos",by.y="V2")
    }      
    
group<-factor(c(rep(1,times=10), rep(2, times=10)))

depth<-merge(depth1, depth2, by="pos")
depth<-depth[,-1]
y<-DGEList(counts=depth, group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

colnames(depth)<-c(paste0("PWS91_",1:10),paste0("PWS07_",1:10) )

depth2<-mapply("*",depth,y$samples["norm.factors"])


plots<-list()
depth2$index<-1:nrow(depth2)
for (i in 1: 10){
    plots[[i]]<-ggplot(depth[,c(1)], aes(x=V2, y=V3))+
        geom_point(size=0.3, alpha=0.4, color="blue")+
        ylab("Read depth")+xlab("chr6 pos")+theme_bw()+
        ggtitle(plist[i])
}



