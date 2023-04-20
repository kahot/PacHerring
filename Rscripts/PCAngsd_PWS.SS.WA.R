## Plot PCAnsgd results for PWS17, SS17, and WA17

#Colors I use for the herring populations
cols<-c("#56b4e9", "#cc79a7","#009e73","#0072b2","#d55e00","#e69f00","#f0e442")

library(ggplot2)
library(cowlplot)

#### PWS SS WA 2017 populations only ####
#using maf05 pruned vcf as a source
pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop3<-pop_info[pop_info$Population.Year %in% c("PWS17","SS17","WA17"),]
pop3<-pop3[,c("Sample","pop","Year.Collected")]
colnames(pop3)[3]<-"year"


C <- as.matrix(read.table(paste0("Data/PCAangsd/PWS.SS.WA2017.cov")))
e <- eigen(C)
pca <-data.frame(Sample=pop3$Sample, 
                 pop= pop3$pop,
                 year=pop3$year,
                 PC1=e$vectors[,1],PC2=e$vectors[,2],
                 PC3=e$vectors[,3],PC4=e$vectors[,4],
                 PC5=e$vectors[,5],PC6=e$vectors[,6],
                 PC7=e$vectors[,7],PC8=e$vectors[,8],
                 stringsAsFactors=FALSE)

prop_explained <- c()
for (s in e$values[1:10]) {
    #print(s / sum(e$values))
    prop_explained <- c(prop_explained,round(((s / sum(e$values))*100),2))
}
#write.csv(pca,"Output/PCA/pca_PWS.WA.SS2017.csv", row.names = F)

pca$pop<-factor(pca$pop, levels=c("PWS","SS","WA"))

ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop), size = 2.5, shape=21)+
    scale_fill_manual(values=paste0(cols[c(2,3,5)],"99"), guide="none")+
    scale_color_manual(values=cols[c(2,3,5)], name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()+
    guides(color = guide_legend(override.aes = list(size=3.5,shape=16, fills=cols[c(2,3,5)])))
ggsave("Output/PCA/pca_PWS.WA.SS2017.png", height = 6, width = 7, dpi=300)

#Plot PC1 to PC3
p1<-ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop), size = 2.5, shape=21)+
    scale_fill_manual(values=paste0(cols[c(2,3,5)],"99"), guide="none")+
    scale_color_manual(values=cols[c(2,3,5)], name="Population", guide="none")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()

p2<- ggplot()+
    geom_point(data = pca, aes(x = PC2, y = PC3, fill = pop, color = pop), size = 2.5, shape=21)+
    scale_fill_manual(values=paste0(cols[c(2,3,5)],"99"), guide="none")+
    scale_color_manual(values=cols[c(2,3,5)], name="Population",guide="none")+
    xlab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()


p3<- ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC3, fill = pop, color = pop), size = 2.5, shape=21)+
    scale_fill_manual(values=paste0(cols[c(2,3,5)],"99"), guide="none")+
    scale_color_manual(values=cols[c(2,3,5)], name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()+
    guides(color = guide_legend(override.aes = list(size=3.5,shape=16, fills=cols[c(2,3,5)])))


png("Output/PCA/PWS.SS.WA.2017_PC1-3.png",width=10, height=3, units="in",res=300)
ggdraw()+
    draw_plot(p1,x=0,y=0, width=0.295,height=1)+
    draw_plot(p2,x=0.3,y=0, width=0.295,height=1)+
    draw_plot(p3,x=0.6,y=0, width=0.41,height=1)
dev.off()    
