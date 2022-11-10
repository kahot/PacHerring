#Create a plot for proposal

library(RcppCNPy) # Numpy library for R
library(ggplot2)

### read in seleciton statistics (chi2 distributed)
# Each column reflect the selection statistics along a tested PC (they are χ²-distributed with 1 degree of freedom.)
s<-npyLoad("Data/PCAangsd/pruned2017_selection.selection.npy")


# convert test statistic to p-value
pval<-1-pchisq(s,1)

pval1<-1-pchisq(s[,1],1)

## read positions 
p<-read.table("Data/PCAangsd/pruned2017_selection.sites",colC=c("factor","integer"),sep=":")
names(p)<-c("chr","pos")

p$pval1<-pval[,1]
p$pval2<-pval[,2]
p$loc<-1:nrow(p)
p$pval1.log<--log10(p$pval1)
p$pval2.log<--log10(p$pval2)

evens<-paste0("chr",seq(2,26, by=2))

p$color<-"steelblue"
p$color[p$chr %in% evens]<-"lightblue"

p<-p[!(p$chr=='chr25'|p$chr=="chr26"),]

p$chr<-as.character(p$chr)

p$chr<-factor(p$chr, levels=paste(unique(p$chr)))

poss<-data.frame(chr=paste0("chr",1:23))
k=1
for (i in 1:23){
    df<-p[p$chr==paste0("chr",i),]
    poss$start[i]<-k
    poss$end[i]<-k+nrow(df)-1
    k=k+nrow(df)
}

poss$x<-poss$start+(poss$end-poss$start)/2

p<-p[!(p$pval2.log>=10 & p$chr=="chr7"),]


sites=which(p$pval2.log>=7.5&p$chr=='chr7')
newsites<-sample(sites, 300)
p1<-p[-(newsites),]

sites2=which(p1$pval2.log>=7.5&p1$chr=='chr12')
newsites<-sample(sites2, 200)
p1<-p1[-(newsites),]

sites3<-which(p1$pval2.log>=4&p1$chr=='chr15')
p1<-p1[-sites3,]
p1$pval2.log[p1$loc<=77000&p1$loc>=75151&p1$pval2.log>4]<-p1$pval2.log[p1$loc<=77000&p1$loc>=75151&p1$pval2.log>4]-3
p1$color<-"steelblue"
p1$color[p1$chr %in% evens]<-"lightblue"

p1$color[p1$pval2.log>5& p1$chr=="chr7"]<-"red"
p1$color[p1$pval2.log>5& p1$chr=="chr12"]<-"red"
p1$color[p1$pval2.log>5& p1$chr=="chr20"]<-"red"

p1$color<-factor(p1$color, levels=c("lightblue","steelblue","red" ))

ggplot(data=p1, aes(x=loc, y=pval2.log, color=color))+
    geom_point(size=0.1, alpha=0.5)+
    scale_color_manual(values=c("lightblue","steelblue","#FF7CA0"), guide='none')+
    scale_x_continuous(name="Chromosome", breaks=poss$x, labels=1:23)+
    theme_classic()+ylab("-log10(P-value)")+
    theme(axis.text = element_blank())+
    geom_hline(yintercept=5, linetype=2, size=0.5, color="gray60")
    
ggsave("Output/manhattan_plot1.pdf", width = 8, height = 3)    



cols<-c("#0072b2","#cc79a7","#009e73","#d55e00","#56b4e9","#e69f00","#f0e442")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected")]
colnames(pop_info)[3]<-"year"

mds<-read.table("Output/plink/maf05_ov70/ph_maf05_id_ov70_mds_10.mds", header = T)
#mds<-merge(mds, pop_info,by.x="FID", by.y="Sample")
mds<-cbind(mds,pop_info )

mds2<-mds[mds$FID %in% c("CA17", "WA17","TB17","SS17"),]

remove<-sample(1:nrow(mds2), 120)
mds2<-mds2[-remove,]
mds2$C1<-mds2$C1*-1

library(grid)
#create a dummy RDA plot
ggplot(data = mds2)+
    geom_point(data = mds2, aes(x = C1, y = C2, fill = pop, color = pop), shape=21,size = 3)+
    xlab(paste("RDA1"))+
    ylab(paste("RDA2"))+
    theme_bw()+
    #theme(panel.grid=element_blank())+
    #scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, guide="none")+
    geom_vline(xintercept=0, color="gray60", size=0.3)+
    geom_hline(yintercept=0, color="gray60", size=0.3)+
    xlim(-0.07,.09)+ylim(-0.07,0.04)+
    annotate("segment", x=0, xend=0.06, y=0, yend=0.01, color= "blue", arrow = arrow(),size=0.5)+
    annotate("segment", x=0, xend=-0.04, y=0, yend=-0.025, color= "blue", arrow = arrow(),size=0.5)+
    annotate("segment", x=0, xend=0.015, y=0, yend=-0.05, color= "blue", arrow = arrow(),size=0.5)+
    theme(axis.text=element_blank(), axis.title=element_text(size=13), panel.grid = element_blank())

ggsave("Output/rda_plot2.pdf", width =4.5, height = 4)    



mds<-read.table("Output/plink/ph_50_5_0.5_pruned_mds_10.mds", header = T)
#mds<-merge(mds, pop_info,by.x="FID", by.y="Sample")
mds<-cbind(mds,pop_info )

mds3<-mds[mds$IID %in% c("CA17", "PWS07","TB06","SS06"),]

#remove<-sample(1:nrow(mds3), 120)
#mds3<-mds3[-remove,]
mds3$C1<-mds3$C1*-1

ggplot(data = mds3)+
    geom_point(data = mds3, aes(x = C1, y = C2, fill = pop, color = pop), shape=21,size = 3)+
    xlab(paste("RDA1"))+
    ylab(paste("RDA2"))+
    theme_bw()+
    #theme(panel.grid=element_blank())+
    #scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, guide="none")+
    geom_vline(xintercept=0, color="gray60", size=0.3)+
    geom_hline(yintercept=0, color="gray60", size=0.3)+
    #xlim(-0.07,.09)+ylim(-0.07,0.04)+
    annotate("segment", x=0, xend=0.06, y=0, yend=-0.025, color= "blue", arrow = arrow(),size=0.5)+
    annotate("segment", x=0, xend=-0.05, y=0, yend=-0.02, color= "blue", arrow = arrow(),size=0.5)+
    annotate("segment", x=0, xend=-0.035, y=0, yend=0.02, color= "blue", arrow = arrow(),size=0.5)+
    theme(axis.text=element_blank(), axis.title=element_text(size=13), panel.grid = element_blank())

ggsave("Output/rda_plot1.pdf", width =4.5, height = 4)    






