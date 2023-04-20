#PCA results from angsd
source("Rscripts/BaseScripts.R")
library(stringr)
library(gridExtra)


pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected")]
colnames(pop_info)[3]<-"year"

cols
"#56b4e9" "#cc79a7" "#009e73" "#0072b2" "#d55e00" "#e69f00" "#f0e442"

##### Look at PCA for each chromosome
Plots<-list()
#no results for chr24
chs<-c(1:23,25:26)
for (i in 1:length(chs)){
    C <- as.matrix(read.table(paste0("Data/PCAangsd/MD7000_maf05_chr",chs[i],".cov")))
    e <- eigen(C)
    pca <-data.frame(Sample=pop_info$Sample, 
                     pop=pop_info$pop,
                     year=pop_info$year,
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
    #write.csv(pca, paste0("Output/PCA/pca_chr",chs[i],".csv"), row.names = F)
    pca$pop<-factor(pca$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))
    Plots[[i]]<-ggplot()+
        geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 3)+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+
        xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
        ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
        theme_bw()+
        ggtitle(paste0("Chr",chs[i]))
    print(i)
}

pdf("Output/PCA/PCA_byChromosome_all.pdf",width = 35, height = 30)
do.call(grid.arrange, c(Plots, ncol=5))
dev.off()


for (i in 1:length(Plots)){
    p<-Plots[[i]]
    if (i>=24) ggsave(paste0("Output/PCA/PCA_Chr",i+1,".png"), plot=p, width=4,height=3.5, dpi=300)
    else ggsave(paste0("Output/PCA/PCA_Chr",i,".png"), plot=p, width=6,height=5, dpi=300)
}



#chr 12
Plots[[12]]
ggsave("Output/PCA/Chr12_allpops.png", width = 5, height = 4, dpi=150)

i=12
ggplot()+
    geom_point(data = pca, aes(x = PC2, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 3)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()+
    ggtitle(paste0("Chr",chs[i]))
ggsave("Output/PCA/Chr12_allpops_PC2.PC3.png", width = 5, height = 4, dpi=150)
ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 3)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()+
    ggtitle(paste0("Chr",chs[i]))
ggsave("Output/PCA/Chr12_allpops_PC1.PC3.png", width = 5, height = 4, dpi=150)


# chr15
Plots[[15]]
ggsave("Output/PCA/Chr15_allpops.png", width = 5, height = 4, dpi=300)

# chr8
Plots[[8]]+
ggsave("Output/PCA/Chr8_allpops.png", width = 5, height = 4, dpi=300)


########
#PCA without TB

pop1<-pop_info[pop_info$pop!="TB",]
plots<-list()
chs<-c(1:23,25:26)
for (i in 1:length(chs)){
    C <- as.matrix(read.table(paste0("Data/PCAangsd/noTB_pruned_chr",chs[i],".cov")))
    e <- eigen(C)
    pca <-data.frame(Sample=pop1$Sample, 
                     pop= pop1$pop,
                     year=pop1$year,
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
    #write.csv(pca, paste0("Output/PCA/noTB_pca_chr",chs[i],".csv"), row.names = F)
    
    pca$pop<-factor(pca$pop, levels=c("PWS","SS", "BC","WA","CA"))
    plots[[i]]<-ggplot()+
        geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 3)+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols[2:6],"99"), guide="none")+
        scale_color_manual(values=cols[2:6], name="Population")+
        xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
        ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
        theme_bw()+
        ggtitle(paste0("Chr",chs[i]))
}
#no sites left after pruning in chr24

pdf("Output/PCA/noTB_PCA_byChromosome.pdf",width = 35, height = 30)
do.call(grid.arrange, c(plots, ncol=5))
dev.off()

for (i in 1:length(plots)){
    p<-plots[[i]]
    if (i>=24) ggsave(paste0("Output/PCA/PCA_noTB_Chr",i+1,".png"), plot=p, width=4,height=3.5, dpi=300)
    else ggsave(paste0("Output/PCA/PCA_noTB_Chr",i,".png"), plot=p, width=6,height=5, dpi=300)
}


### All pops combined

C <- as.matrix(read.table(paste0("Data/PCAangsd/AllPH_maf05_pruned.cov")))
e <- eigen(C)
pca <-data.frame(Sample=pop_info$Sample, 
                 pop= pop_info$pop,
                 year=pop_info$year,
                 PC1=e$vectors[,1],PC2=e$vectors[,2],
                 PC3=e$vectors[,3],PC4=e$vectors[,4],
                 PC5=e$vectors[,5],PC6=e$vectors[,6],
                 PC7=e$vectors[,7],PC8=e$vectors[,8],
                 PC9=e$vectors[,9],PC10=e$vectors[,10],
                 stringsAsFactors=FALSE)

prop_explained <- c()
for (s in e$values[1:10]) {
    #print(s / sum(e$values))
    prop_explained <- c(prop_explained,round(((s / sum(e$values))*100),2))
}
#write.csv(pca,"Output/PCA/pca_allPops.csv", row.names = F)
#pca<-read.csv("Output/PCA/pca_allPops.csv")

pca$pop<-factor(pca$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))

ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"66"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()
ggsave("Output/PCA/pca_allPops.pdf", height = 6, width = 7)
ggsave("Output/PCA/pca_allPops.png", height = 6, width = 7, dpi=300)


#different PC axes:
pcas<-list()
pcas[[1]]<-ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 2)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"4D"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()

pcas[[2]]<-ggplot()+
    geom_point(data = pca, aes(x = PC2, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"4D"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()
pcas[[3]]<-ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"4D"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()
pcas[[4]]<-ggplot()+
    geom_point(data = pca, aes(x = PC3, y = PC4, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"4D"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    ylab(paste("PC 4: ", prop_explained[4],"%\n",sep = ""))+
    theme_bw()

pdf("Output/PCA/pca_axes_allPh.pdf", width = 20, height=4)
do.call(grid.arrange, c(pcas, ncol=4))
dev.off()

### 2017 populations
C <- as.matrix(read.table(paste0("Data/PCAangsd/Y2017_pruned.cov")))
e <- eigen(C)
pca <-data.frame(Sample=pop_info$Sample[pop_info$year==2017], 
                 pop= pop_info$pop[pop_info$year==2017],
                 year=pop_info$year[pop_info$year==2017],
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
write.csv(pca,"Output/PCA/pca_Y2017.csv", row.names = F)
pca$pop<-factor(pca$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))
ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"66"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()
ggsave("Output/PCA/pca_Y2017.pdf", height = 6, width = 7)


### PWS & SS
C <- as.matrix(read.table(paste0("Data/PCAangsd/PWS_SS_maf05_pruned.cov")))
e <- eigen(C)
pwsss<-pop_info[pop_info$pop=="PWS"|pop_info$pop=="SS",]
pca <-data.frame(Sample=pwsss$Sample, 
                 pop= pwsss$pop,
                 year=pwsss$year,
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
write.csv(pca,"Output/PCA/pca_PWS_SS_pruned.csv", row.names = F)
pca$pop<-factor(pca$pop, levels=c("PWS","SS"))

ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols[2:3],"66"), guide="none")+
    scale_color_manual(values=cols[2:3], name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()
ggsave("Output/PCA/pca_PWS_SS_pruned.pdf", height = 6, width = 7)

#different axes
pcas<-list()
pcas[[1]]<-ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 2)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols[2:3],"4D"), guide="none")+
    scale_color_manual(values=cols[2:3], name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    theme_bw()
pcas[[2]]<-ggplot()+
    geom_point(data = pca, aes(x = PC2, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols[2:3],"4D"), guide="none")+
    scale_color_manual(values=cols[2:3], name="Population")+
    xlab(paste("PC 2: ", prop_explained[2],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()
pcas[[3]]<-ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols[2:3],"4D"), guide="none")+
    scale_color_manual(values=cols[2:3], name="Population")+
    xlab(paste("PC 1: ", prop_explained[1],"%\n",sep = ""))+
    ylab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    theme_bw()
pcas[[4]]<-ggplot()+
    geom_point(data = pca, aes(x = PC3, y = PC4, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols[2:3],"4D"), guide="none")+
    scale_color_manual(values=cols[2:3], name="Population")+
    xlab(paste("PC 3: ", prop_explained[3],"%\n",sep = ""))+
    ylab(paste("PC 4: ", prop_explained[4],"%\n",sep = ""))+
    theme_bw()
pdf("Output/PCA/pca_axes_PWS_SS.pdf", width = 20, height=4)
do.call(grid.arrange, c(pcas, ncol=4))
dev.off()



#Find which individuals of PWS were outliers
outliers1<-pca[pca$PC1<=(-0.1),]
#write.csv(outliers1, "Output/PCA/PWS_SS_pcangsd_outliers.csv")



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
write.csv(pca,"Output/PCA/pca_PWS.WA.SS2017.csv", row.names = F)
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
ggsave("Output/PCA/pca_PWS.WA.SS2017.2.png", height = 5, width = 6.1, dpi=300)

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

library(cowplot)
png("Output/PCA/PWS.SS.WA.2017_PC1-3.png",width=10, height=3, units="in",res=300)
ggdraw()+
    draw_plot(p1,x=0,y=0, width=0.295,height=1)+
    draw_plot(p2,x=0.3,y=0, width=0.295,height=1)+
    draw_plot(p3,x=0.6,y=0, width=0.41,height=1)
dev.off()    
