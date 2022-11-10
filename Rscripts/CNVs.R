#Copy Number Variant Check
source("Rscripts/BaseScripts.R")


library(gridExtra)
library(zoo)
library(vcfR)

pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

pops_info$Sample[pops_info$Population.Year=="PWS96"]

## Create 2 bed files for each chromosome
chr<-read.table("Data/vcfs/chr_sizes.bed", header = F)

for (i in 1:26){
    df<-chr[chr$V1==paste0("chr",i),]
    mid<-ceiling(df$V3/2)
    new<-data.frame(ch=rep(paste0("chr",i),times=2), st=c(1,mid+1),en=c(mid, df$V3[1]))
    write.table(new[1,], paste0("Data/bam_depth/bed/chr",i,"_1.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(new[2,], paste0("Data/bam_depth/bed/chr",i,"_2.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
}

#Run samptools to extract read depth from individual bam files
samtools depth -b /home/ktist/ph/data/bam_depth/bed/chr6_2.bed /home/eoziolor/phpopg/data/align/1151_PWS07.bam > /home/ktist/ph/data/bam_depth/1151_PWS07_chr6-2.txt
gzip /home/ktist/ph/data/bam_depth/1151_PWS07_chr6-2.txt

# transfer to local and uncompress the file

#Read the file to plot for PWS pops near CYP1A

pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws96<-pops_info$Sample[pops_info$Population.Year=="PWS96"]
pws91<-pops_info$Sample[pops_info$Population.Year=="PWS91"]
pws07<-pops_info$Sample[pops_info$Population.Year=="PWS07"]
ca<-pops_info$Sample[pops_info$Population.Year=="CA17"]

pops<-c("pws91","pws96","pws07","ca")

plots<-list()

for (p in 1:length(pops)){
    plist<-get(paste0(pops[p]))
    plots<-list()
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",plist[i],"_chr6-cyp.txt"), header=F, sep="\t")
        plots[[i]]<-ggplot(df[df$V2>=2850000&df$V2<=3020000,], aes(x=V2, y=V3))+
            geom_point(size=0.3, alpha=0.4, color="blue")+
            ylab("Read depth")+xlab("chr6 pos")+theme_bw()+
            ggtitle(plist[i])
    }
    
    pdf(paste0("Output/CNV/chr6/", pops[p],".pdf"), width =8, height = 10)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}

#plot the entire length 2-3.5mb
for (p in 1:3){
    plist<-get(paste0(pws[p]))
    plots<-list()
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",plist[i],"_chr6-cyp.txt"), header=F, sep="\t")
        plots[[i]]<-ggplot(df, aes(x=V2, y=V3))+
            geom_point(size=0.1, alpha=0.4, color="gray60")+
            ylab("Read depth")+xlab("chr6 pos")+theme_bw()+
            ggtitle(plist[i])
    }
    pdf(paste0("Output/CNV/chr6/", pws[p],"_2-3.5mb.pdf"), width =12, height = 15)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}

### Look for ahr pathway genes

geneset<-read.csv("Output/ahR/ahRpathway_geneset.csv", row.names = 1)

#create bed files for the genes in the genest
for (i in 1:nrow(geneset)){
    st<-geneset$start[i]-100000
    en<-geneset$end[i]+100000
    df<-data.frame(chr=geneset$chr[i], st=st,end=en)
    write.table(df, paste0("Output/CNV/ahr/bed/",geneset$gene_name[i],".bed"), row.names = F, col.names = F, quote=F, sep="\t")
}




pops<-c("pws91","pws96","pws07","ca")


#ahr1b
for (p in 1:length(pops)){
    plist<-get(paste0(pops[p]))
    plots<-list()
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",plist[i],"_ahr1b.txt"), header=F, sep="\t")
        plots[[i]]<-ggplot(df, aes(x=V2, y=V3))+
            geom_point(size=0.3, alpha=0.4, color="blue")+
            ylab("Read depth")+xlab("chr2 pos")+theme_bw()+
            ggtitle(plist[i])+
            geom_vline(xintercept = 15770067, color="gray70", size=0.5)+
            geom_vline(xintercept = 15790033, color="gray70", size=0.5)
            
    }
    
    #pdf(paste0("Output/CNV/ahr/ahr1b_", pops[p],".pdf"), width =8, height = 10)
    #do.call(grid.arrange, c(plots, ncol=2))
    #dev.off()
    png(paste0("Output/CNV/ahr/ahr1b_", pops[p],".png"), width =8, height = 10, unit="in", res=100)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
}

#ahrrb
for (p in 1:length(pops)){
    plist<-get(paste0(pops[p]))
    plots<-list()
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",plist[i],"_ahrrb.txt"), header=F, sep="\t")
        plots[[i]]<-ggplot(df, aes(x=V2, y=V3))+
            geom_point(size=0.3, alpha=0.4, color="steelblue")+
            ylab("Read depth")+xlab("chr22 pos")+theme_bw()+
            ggtitle(plist[i])+
            geom_vline(xintercept = 18844730, color="gray70", size=0.5)+
            geom_vline(xintercept = 18872212, color="gray70", size=0.5)
    }
    
    png(paste0("Output/CNV/ahr/ahrrb_", pops[p],".png"), width =8, height = 10,unit="in", res=100)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}


#aip
for (p in 1:length(pops)){
    plist<-get(paste0(pops[p]))
    plots<-list()
    for (i in 1: 10){
        df<-read.csv(paste0("data/bam_depth/",plist[i],"_aip.txt"), header=F, sep="\t")
        plots[[i]]<-ggplot(df, aes(x=V2, y=V3))+
            geom_point(size=0.3, alpha=0.4, color="steelblue")+
            ylab("Read depth")+xlab("chr22 pos")+theme_bw()+
            ggtitle(plist[i])+
            geom_vline(xintercept = 21751751, color="gray70", size=0.5)+
            geom_vline(xintercept = 21758253, color="gray70", size=0.5)
    }
    
    png(paste0("Output/CNV/ahr/ahrrb_", pops[p],".png"), width =8, height = 10,unit="in", res=100)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}


geneset

for (g in 1: nrow(geneset)){
    gene<-geneset$gene_name[g]
    for (p in 1:length(pops)){
        plist<-get(paste0(pops[p]))
        plots<-list()
        for (i in 1: 10){
            df<-read.csv(paste0("data/bam_depth/",plist[i],"_",gene,".txt"), header=F, sep="\t")
            plots[[i]]<-ggplot(df, aes(x=V2, y=V3))+
                geom_point(size=0.3, alpha=0.4, color="steelblue")+
                ylab("Read depth")+xlab(paste(geneset$chr[g],"pos"))+theme_bw()+
                ggtitle(plist[i])+
                geom_vline(xintercept = geneset$start[g], color="gray70", size=0.5)+
                geom_vline(xintercept = geneset$end[g], color="gray70", size=0.5)
        }
        
        png(paste0("Output/CNV/ahr/",gene,"_", pops[p],".png"), width =8, height = 10,unit="in", res=100)
        do.call(grid.arrange, c(plots, ncol=2))
        dev.off()
        
    }
    
}

### Use vcf files to detect CNV?
library(vcfR)
#Look at if copy number variation exist near CYP1A like the killifish
vcf <- read.vcfR('Data/new_vcf/PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz')

#extract the read depth information 
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
write.csv(dp, "Output/CNV/MD7000_330k_maf05_read_depth.csv")

dp<-read.csv("Output/CNV/MD7000_maf05_read_depth.csv", row.names = 1)
#row.names(dp)<-dp[,1]
#dp<-dp[,-1]
colnames(dp)<-gsub("X",'',colnames(dp))

#divide by chromosome?
for (i in 1:26){
    df<-dp[grep(paste0("chr",i,"_"),rownames(dp)),]
    df<-as.data.frame(df)
    df$pos<-as.integer(gsub(paste0("chr",i,"_"), "", rownames(df)))
    #dname<-paste0("dp",i)
    #assign(dname,df)
    write.csv(df,paste0("Output/CNV/chr",i,"_330k_depth.csv"))
}

## read in the files

for (i in 1:26){
    df<-read.csv(paste0("Output/CNV/chr",i,"_depth.csv"), row.names = 1)
    dname<-paste0("dp",i)
    assign(dname,df)
}

## 1. look at the PWS populations
#pws1<-dp1[,grep("PWS",colnames(dp1))]
#pws1$pos<-dp1$pos
##select only PWS91
#p91<-pws1[,grep("PWS91", colnames(pws1))]
#p91$pos<-dp1$pos

#roll_df<-data.frame()
#for (i in 1:40){
#    rollm<-rollmean(p91[,i], k=1000, na.rm=T, align="center")
#    roll<-data.frame(roll_ave=rollm, sample=rep(colnames(p91[i]), times=length(rollm)))
#    roll_df<-rbind(roll_df, roll)
#}


#roll_df$position<-p91$pos[500:(nrow(p91)-500)]
#ggplot(roll_df, aes(x=position, y=roll_ave))+
#    facet_wrap(~sample, ncol=3)+
#    geom_point(size=0.5, color="gray60")+
#    theme_bw()
#ggsave("Output/CNV/pws91_chr1.png", height = 20, width=12, dpi=100, units="in")


pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)
#for (j in c(1,7,8,9,14)){

ch=2
for (j in 2: length(pops)){
    #dt<-dp1[,grep(pops[j],colnames(dp1))]
    dt<-dp2[,grep(pops[j],colnames(dp2))]
    
    #calculate rolling average
    mvAve<-apply(dt[,1:40],2, function(x) rollmean(x, k=200, na.rm=T, align="center"))
    mvAve<-data.frame(mvAve)
    #dt$pos<-dp1$pos
    dt$pos<-dp2$pos
    
    #plot the rolling average individually
    mvAve$position<-dt$pos[100:(nrow(dt)-100)]
    mvAvem<-melt(mvAve,id.vars="position")
    mvAvem$variable<-gsub("X",'',mvAvem$variable)
    ggplot(mvAvem, aes(x=position, y=value))+
        facet_wrap(~variable, ncol=3)+
        geom_point(size=0.3, color="gray60")+
        theme_bw()
    ggsave(paste0("Output/CNV/chr2_mv_", pops[j],".png"), height = 20, width=12, dpi=100, units="in")
    
    #plot the depth individually
    #png(paste0("Output/CNV/chr1_dp_",pops[j], ".png"), width = 16, height=16,units = "in", res=100)
    #par(mfrow=c(10,4))
    #par(mar = c(2, 4, 1, 1))
    #for (i in 1:40){
    #    plot(dt$pos, dt[,i], pch=".", xlab="", ylab=colnames(dt)[i])
    #}
    #dev.off()
}

for (j in c(3:6)){
    dt<-dp1[,grep(pops[j],colnames(dp1))]
    #calculate rolling average
    mvAve<-apply(dt[,1:40],2, function(x) rollmean(x, k=200, na.rm=T, align="center"))
    mvAve<-data.frame(mvAve)
    dt$pos<-dp1$pos
    #plot the rolling average individually
    mvAve$position<-dt$pos[100:(nrow(dt)-100)]
    mvAvem<-melt(mvAve,id.vars="position")
    mvAvem$variable<-gsub("X",'',mvAvem$variable)
    ggplot(mvAvem, aes(x=position, y=value))+
        facet_wrap(~variable, ncol=3)+
        geom_point(size=0.3, color="gray60")+
        theme_bw()
    ggsave(paste0("Output/CNV/chr1_mv_", pops[j],".png"), height = 20, width=12, dpi=100, units="in")
}



for (j in 1: length(pops)){
    dt<-dp1[,grep(pops[j],colnames(dp1))]
    dt$pos<-dp1$pos
    png(paste0("Output/CNV/chr1_dp_",pops[j], ".png"), width = 16, height=16,units = "in", res=100)
    par(mfrow=c(10,4))
    par(mar = c(2, 4, 1, 1))
    for (i in 1:40){
        plot(dt$pos, dt[,i], pch=".", xlab="", ylab=colnames(dt)[i])
    }
    dev.off()
}




dt$pos[1000]-dt$pos[1]#1,908,909
dt$pos[200]-dt$pos[1]#186,638

png(paste0("Output/CNV/chr1_dp_",pops[j], ".png"), width = 16, height=16,units = "in", res=100)
par(mfrow=c(10,4))
par(mar = c(2, 4, 1, 1))
for (i in 1:40){
    plot(dt$pos, dt[,i], pch=".", xlab="", ylab=colnames(dt)[i])
}
dev.off()
plot(dt$pos,dt$`0756_PWS96`, pch=".", xlab="", yalb=)


