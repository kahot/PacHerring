#Pacific Herring SFS and resulting Fst from ANGSD/realSFS 
# Fst was calculated using folded 2dsfs 
source("Rscripts/BaseScripts.R")

#1. PWS
### Fst from 2Dsfs
pops<-c("PWS91","PWS96","PWS07","PWS17")
comb<-t(combn(pops,2))

fst<-data.frame()
for (i in 1: nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    df<-read.delim(paste0("Data/new_vcf/angsd/fromVCF/2D/fst_",pop1,"_",pop2,"_50kWindow_maf00"))
    conames<-colnames(df)[2:4]
    colnames(df)[4]<-"Fst"
    colnames(df)[1:3]<-conames
    df$pop<-paste0(pop1,".vs.",pop2)
    df$ch=as.integer(gsub("chr","", df$chr))
    df<-df[order(df$ch),]
    df$loc<-1:nrow(df)
    fst<-rbind(fst, df)
}

# Plost Fst for all PWS years (distant as PWS91)
evens<-paste0("chr",seq(2,26, by=2))
fst$color<-"col1"
fst$color[fst$chr %in% evens]<-"col2"
fst$pop<-factor(fst$pop, levels=unique(fst$pop))

#add chromosome number
df<-fst[fst$pop=="PWS91.vs.PWS96",]
rows<-data.frame(chr=1:26)
for (i in 1:26){
    if (i ==1){
        rows$n[i]<-nrow(df[df$ch==i,])
        rows$middle[i]<--nrow(df[df$ch==i,])/2
    }
    if (i >1){
        rows$n[i]<-nrow(df[df$ch==i,])
        rows$middle[i]<-sum(rows$n[1:(i-1)])+rows$n[i]/2
    }
}

ggplot(fst, aes(x=loc, y=Fst, color=color))+
    facet_wrap(~pop, ncol = 1, strip.position="right")+
    geom_point(size=0.2)+
    scale_color_manual(values=c("gray50","steelblue"))+
    theme_bw()+
    ylab("Fst")+xlab('Genome position')+
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
    scale_x_continuous(breaks=rows$middle, labels=1:26)
#ggsave("Output/fst_pbs/PWS_Fst_pairwise_comparison.pdf", width = 20, height = 10)
ggsave("Output/fst_pbs/PWS_Fst_pairwise_comparison.png", width = 20, height = 10, dpi=150)





## Plot along each chromosome
fst$chr<-factor(fst$chr, levels=paste0("chr",1:26))
fstpw<-fst
plots<-list()
compare<-paste0(unique(fstpw$pop))
max(fstpw$Fst)
for (i in 1:6){ 
    fs<-gsub("vs.","",compare[i])
    pops <- unlist(strsplit(fs, "\\."))
    cat ("Comparison between PWS",pops[1]," vs.PWS", pops[2],"\n")
    maxy<-max(fstpw$Fst[fstpw$pop==compare[i]])
    # Fst with actual line to highlight the differences
    plots[[i]] <- ggplot(fstpw[fstpw$pop==compare[i],], aes(x =midPos, y =Fst )) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+ylim(0,maxy+0.02)+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste0("PWS",pops[1]," vs.", pops[2]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
}

png(paste0("Output/fst_pbs/MD7000/PWS_Fst_maf00.png"), height = 8, width = 18, res=150, units = "in")
grid.arrange(plots[[3]], plots[[2]], plots[[4]], plots[[1]],plots[[5]],plots[[6]], ncol=3)
dev.off()




#plot heatmap
fstmats<-list.files("Output/fst_pbs/", pattern="Fst_matrix.csv")
for (i in 1:3){
    df<-read.csv(paste0("Output/fst_pbs/",fstmats[i]), row.names = 1)
    diag(df)<-0
    df$pop<-rownames(df)
    dfm<-melt(df,na.rm=T, id.vars='pop')
    
    #NA to diagonal
    dfm$value[dfm$value==0]<-NA
    dfm$pop<-factor(dfm$pop, levels=pops)
    dfm$value<-round(dfm$value, 4)
    ggplot(data = dfm, aes(pop, variable, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradientn(colors=c("white", "#0C54FF"), limits=c(0, (max(dfm$value, na.rm=T)+0.005)),na.value="gray80", 
                             name="Fst")+
        theme_minimal()+ xlab("")+ylab("")+
        theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                         size = 12, hjust = 0.5))+
        theme(axis.text.y = element_text(size = 12))+
        coord_fixed()+
        geom_text(aes(pop, variable, label = value), color = "black", size = 5)
    fname<-gsub(".csv",'',fstmats[i])
    ggsave(paste0("Output/fst_pbs/",fname,".pdf"), width = 5, height = 5)
    
}

#chromosome by chromosome comparison
plots<-list()
Fst<-data.frame()
for (j in 1:26){
    fst.ch<-fst[fst$ch==j,]
    
    pairch<-data.frame(matrix(ncol=4, nrow=4), row.names=pops)
    colnames(pairch)<-pops
    for (i in 1:6){
        pop1<-comb[i,1]
        pop2<-comb[i,2]
        df<-fst.ch[fst.ch$pop==pwss[i],]
        pairch[pop1,pop2]<-mean(df$Fst, na.rm=T)
    }
    diag(pairch)<-0
    pairch$pop<-rownames(pairch)
    dfm<-melt(pairch,na.rm=T, id.vars='pop')
    
    #NA to diagonal
    dfm$value[dfm$value==0]<-NA
    dfm$pop<-factor(dfm$pop, levels=pops)
    dfm$value<-round(dfm$value, 4)
    plots[[j]]<-ggplot(data = dfm, aes(pop, variable, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradientn(colors=c("white", "#0C54FF"), limits=c(0, (max(dfm$value, na.rm=T)+0.005)),na.value="gray80", 
                             name="Fst")+
        theme_minimal()+ xlab("")+ylab("")+
        theme(axis.text.x = element_text(angle = 45, vjust = 0, 
                                         size = 9, hjust = 0.5))+
        theme(axis.text.y = element_text(size = 9))+
        coord_fixed()+ggtitle(paste0("Chr",j))+
        geom_text(aes(pop, variable, label = value), color = "black", size = 2.5)
    dfm$chr<-j
    Fst<-rbind(Fst, dfm)
    
}

pdf("Output/fst_pbs/PWS_pairwiseFst_byChrom.pdf", width = 18, height = 12)
do.call(grid.arrange, c(plots, ncol=6))
dev.off()

Fst$id<-paste0(Fst$pop," vs.",Fst$variable)
Fst<-Fst[!is.na(Fst$value),]
ggplot(Fst, aes(x=chr, y=value,color=id))+
    geom_point()+
    geom_path(stat="identity")+
    theme_minimal()+ylab("Fst")+
    scale_x_continuous(breaks=1:26, labels = 1:26)+
    theme(legend.title = element_blank(), panel.grid.minor.x = element_blank())
ggsave("Output/fst_pbs/PWS_Fst_byChromosome_dotplot.pdf", width = 13, height=6.5)








# MD7000 TB    

tb1 <- read.delim("Data/new_vcf/angsd/fromVCF/2D/fst_TB06_TB_91_TB96_50kWindow")
tb2 <- read.delim("Data/new_vcf/angsd/fromVCF/2D/fst_TB06_TB17_TB91_50kWindow")
tb3 <- read.delim("Data/new_vcf/angsd/fromVCF/2D/fst_TB06_TB17_TB96_50kWindow")

colnames(tb1)[5:7]<-c("Fst06.91","Fst06.96","Fst91.96")
colnames(tb2)[5:7]<-c("Fst06.17","Fst06.91","Fst17.91")
colnames(tb3)[5:7]<-c("Fst06.17","Fst06.96","Fst17.96")

colnames(tb1)[8:10]<-c("PBS06","PBS91","PBS96")
colnames(tb2)[8:10]<-c("PBS06","PBS17","PBS91")
colnames(tb3)[8:10]<-c("PBS06","PBS17","PBS96")


#Merge Fst values
fsttb<-merge(tb1[,1:7], tb2[,c(1,5,7)], by="region")
fsttb<-merge(fsttb, tb3[,c(1,7)], by="region")
fsttb$ch<-as.integer(gsub("chr","",fsttb$chr))

fsttb<-fsttb[order(fsttb$ch, fsttb$midPos),]

fsttb$chr<-factor(fsttb$chr, levels=paste0("chr",1:26))
plots<-list()
#TB91 vs 96 has extremely high FSt in chr3 -> ymax =0.3
#everything else ymax=0.2
for (i in 1:6){ 
    fs<-gsub("Fst","",colnames(fsttb)[i+4])
    pops <- unlist(strsplit(fs, "\\."))
    cat ("Comparison between TB",pops[1]," vs.TB", pops[2],"\n")
    
    # Fst with actual line to highlight the differences
    ma=ifelse (i == 3, 0.3, 0.2)
    
    plots[[i]] <- ggplot(fsttb, aes_string(x = 'midPos', y = paste0(colnames(fsttb)[i+4]))) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+ylim(0,ma)+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste0("TB",pops[1]," vs.", pops[2]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
}

png(paste0("Output/fst_pbs/MD7000/TB_Fst2.png"), height = 8, width = 18, res=150, units = "in")
grid.arrange(plots[[3]], plots[[2]], plots[[4]], plots[[1]],plots[[5]],plots[[6]], ncol=3)
dev.off()

#PBS
tbpbs<-merge(tb1, tb2[,c("region", "PBS17")], by="region")
tbpbs$ch<-as.integer(gsub("chr",'',tbpbs$chr))

tbpbs<-tbpbs[order(tbpbs$ch, tbpbs$midPos),]
tbpbs$loc<-1:nrow(tbpbs)

evens<-paste0("chr",seq(2,26, by=2))
tbpbs$color<-"col1"
tbpbs$color[tbpbs$chr %in% evens]<-"col2"

tbpbsm<-melt(tbpbs[c("chr","loc","color","PBS96","PBS06","PBS17")], id.vars=c("chr","loc","color"))
ggplot(tbpbsm, aes(x=loc, y=value, color=color))+
    facet_wrap(~variable, ncol = 1, strip.position="right")+
    geom_point(size=0.2)+
    scale_color_manual(values=c("gray30",grb))+
    theme_bw()+
    ylab("PBS")+xlab('Genome position')+
    theme(legend.position = "none")+
    theme(axis.text.x = element_blank())
ggsave("Output/fst_pbs/MD7000/TB_pbs_3yrs_stacked.pdf", width = 10, height = 7)




#MD7000 SS
ss1 <- read.delim("Data/new_vcf/angsd/Fst/fst_SS06_SS17_SS96_50kWindow")
colnames(ss1)[5:7]<-c("Fst06.17","Fst06.96","Fst17.96")
ss1$ch<-as.integer(gsub("chr","",ss1$chr))

fst<-ss1[order(ss1$ch, ss1$midPos),]
fst$loc<-1:nrow(fst)

fst$chr<-factor(fst$chr, levels=paste0("chr",1:26))
fstss<-fst
plots<-list()
for (i in 1:3){ 
    fs<-gsub("Fst","",colnames(fstss)[i+4])
    pops <- unlist(strsplit(fs, "\\."))
    cat ("Comparison between SS",pops[1]," vs.SS", pops[2],"\n")
    
    # Fst with actual line to highlight the differences
    plots[[i]] <- ggplot(fstss, aes_string(x = 'midPos', y = paste0(colnames(fstss)[i+4]))) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+ylim(0,0.155)+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste0("SS",pops[1]," vs.", pops[2]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
}

png(paste0("Output/fst_pbs/MD7000/SS_Fst.png"), height = 4, width = 16, res=150, units = "in")
grid.arrange(plots[[2]], plots[[1]], plots[[3]], ncol=3)
dev.off()



#Year2017
fst1 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_CA17_PWS17_50kWindow")
colnames(fst1)[5:7]<-c("Fst_BC.CA","Fst_BC.PWS","Fst_CA.PWS")

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_CA17_SS17_50kWindow")
colnames(fst2)[6:7]<-c("Fst_BC.SS","Fst_CA.SS")
fst1<-cbind(fst1[,1:7],fst2[,6:7])

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_CA17_TB17_50kWindow")
colnames(fst2)[6:7]<-c("Fst_BC.TB","Fst_CA.TB")
fst1<-cbind(fst1,fst2[,6:7])

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_CA17_WA17_50kWindow")
colnames(fst2)[6:7]<-c("Fst_BC.WA","Fst_CA.WA")
fst1<-cbind(fst1,fst2[,6:7])

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_PWS17_SS17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_PWS.SS"

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_PWS17_TB17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_PWS.TB"

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_PWS17_WA17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_PWS.WA"

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_CA17_TB17_WA17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_TB.WA"

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_BC17_SS17_TB17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_SS.TB"

fst2 <- read.delim("Data/new_vcf/angsd/Fst/fst_CA17_SS17_WA17_50kWindow")
fst1<-cbind(fst1,fst2[,7])
colnames(fst1)[ncol(fst1)]<-"Fst_SS.WA"

write.csv(fst1,"Output/fst_pbs/MD7000/Fst_window_year2017_combined.csv")

fst1$ch<-as.integer(gsub("chr","",fst1$chr))

fst<-fst1[order(fst1$ch, fst1$midPos),]

fst$chr<-factor(fst$chr, levels=paste0("chr",1:26))
fst17<-fst

plots<-list()
for (i in 1:nrow(comb4)){ 
    fs<-gsub("Fst_","",colnames(fst17)[i+4])
    pops <- unlist(strsplit(fs, "\\."))
    cat (paste0("Comparison between ",pops[1],"17 vs.", pops[2],"17\n"))
    
    # Fst with actual line to highlight the differences
    plots[[i]] <- ggplot(fst17, aes_string(x = 'midPos', y = paste0(colnames(fst17)[i+4]))) + 
        geom_point(size = 1, color = "gray",alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste0(pops[1],"17 vs.", pops[2],"17"))+
        geom_line(color="steelblue", size=0.2)+
        facet_wrap(~chr, ncol = 9)
}

png(paste0("Output/fst_pbs/MD7000/Year2017_Fst1.png"), height = 8, width = 20, res=150, units = "in")
do.call(grid.arrange, c(plots[1:6], ncol=3))
dev.off()


png(paste0("Output/fst_pbs/MD7000/Year2017_Fst2.png"), height = 12, width = 20, res=150, units = "in")
do.call(grid.arrange, c(plots[7:15], ncol=3))
dev.off()

#TB comb
colnames(fst17)[5:19]

png(paste0("Output/fst_pbs/MD7000/Year2017_Fst_TB.png"), height = 8, width = 20, res=150, units = "in")
grid.arrange(plots[[6]], plots[[7]], plots[[11]],plots[[13]],plots[[14]], ncol=3)
dev.off()
#CA comb
png(paste0("Output/fst_pbs/MD7000/Year2017_Fst_CA.png"), height = 8, width = 20, res=150, units = "in")
grid.arrange(plots[[1]], plots[[3]], plots[[5]],plots[[7]],plots[[9]], ncol=3)
dev.off()


#Plot MD7000 mean Fst
## Fst matrix temporal samples (within population comparisons)
# temporal samples per populations
pops<-c("PWS","SS","TB")
popsname<-c("pw","ss","tb")
for (i in 1:length(pops)){
    p<-pops[i]
    fst<-get(paste0("fst",popsname[i]))
    if (p=="PWS") y3="07"
    else y3="06"
    if (p=="SS"){
        p96xp06 <- round(mean(fst$Fst06.96),4)
        p96xp17 <- round(mean(fst$Fst17.96),4)
        p07xp17 <- round(mean(fst$Fst06.17),4)
        
        fst_vec <- c(0,p96xp06,p96xp17,
                     p96xp06,0,p07xp17,
                     p96xp17,p07xp17,0)
        fst_mat = matrix(fst_vec, nrow = 3, ncol = 3)
        colnames(fst_mat) <- c("SS96","SS06","SS17")
        rownames(fst_mat) <- c("SS96","SS06","SS17")
        
    }
    else{
        p91xp96 <- round(mean(fst$Fst91.96),4)
        p91xp17 <- round(mean(fst$Fst17.91),4)
        p96xp17 <- round(mean(fst$Fst17.96),4)
        if (i==3){
            p91xp07 <- round(mean(fst$Fst06.91),4)
            p96xp07 <- round(mean(fst$Fst06.96),4)
            p07xp17 <- round(mean(fst$Fst06.17),4)
        }
        if (i==1){
            p91xp07 <- round(mean(fst$Fst07.91),4)
            p96xp07 <- round(mean(fst$Fst07.96),4)
            p07xp17 <- round(mean(fst$Fst07.17),4)
        }
        
        fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
                     p91xp96,0,p96xp07,p96xp17,
                     p91xp07,p96xp07,0,p07xp17,
                     p91xp17,p96xp17,p07xp17,0)
        
        fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
        colnames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
        rownames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
    }
    
    fst_mat[lower.tri(fst_mat)]<- NA
    write.csv(fst_mat,paste0("Output/fst_pbs/MD7000/",p,"_Fst_matirx.csv"))
}


fstmats<-list.files("Output/fst_pbs/MD7000/", pattern="Fst_matirx.csv")
for (i in 1:3){
    df<-read.csv(paste0("Output/fst_pbs/MD7000/",fstmats[i]))
    df<-read.csv(paste0("Output/fst_pbs/MD7000/PWS_Fst_matirx.csv"))
    dfm<-melt(df,na.rm=T)
    #NA to diagonal
    dfm$value[dfm$value==0]<-NA
    
    dfm$X<-factor(dfm$X, levels=c(levels(dfm$variable)))
    ggplot(data = dfm, aes(X, variable, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradientn(colors=c("white", "blue"), limits=c(0, (max(dfm$value, na.rm=T)+0.005)),na.value="gray80", 
                             name="Fst")+
        theme_minimal()+ xlab("")+ylab("")+
        theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                         size = 12, hjust = 0.5))+
        theme(axis.text.y = element_text(size = 12))+
        coord_fixed()+
        geom_text(aes(X, variable, label = value), color = "black", size = 5)
    fname<-gsub(".csv",'',fstmats[i])
    ggsave(paste0("Output/fst_pbs/MD7000/",fname,".pdf"), width = 5, height = 5)
    
}




## 2017 between populations Fst matrix
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_CA17_PWS17_50kWindow")
bcxca <- round(mean(fst$Fst01),4)
bcxpw <- round(mean(fst$Fst02),4)
caxpw <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_CA17_WA17_50kWindow")
bcxwa <- round(mean(fst$Fst02),4)
caxwa <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_BC17_SS17_TB17_50kWindow")
bcxss <- round(mean(fst$Fst01),4)
bcxtb <- round(mean(fst$Fst02),4)
ssxtb <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_CA17_PWS17_SS17_50kWindow")
caxpw <- round(mean(fst$Fst01),4)
caxss <- round(mean(fst$Fst02),4)
pwxss <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_CA17_PWS17_TB17_50kWindow")
caxtb <- round(mean(fst$Fst02),4)
pwxtb <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_PWS17_SS17_WA17_50kWindow")
pwxwa <- round(mean(fst$Fst02),4)
ssxwa <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_SS17_TB17_WA17_50kWindow")
tbxwa <- round(mean(fst$Fst12),4)



fst_vec <- c(0,pwxtb,ssxtb,bcxtb,tbxwa,caxtb,
             pwxtb,0,pwxss,bcxpw,pwxwa,caxpw,
             ssxtb,pwxss,0,bcxss,ssxwa,caxss,
             bcxtb,bcxpw,bcxss,0,bcxwa,bcxca,
             tbxwa,pwxwa,ssxwa,bcxwa,0,caxwa,
             caxtb,caxpw,caxss,bcxca,caxwa,0)

fst_mat = matrix(fst_vec, nrow = 6, ncol = 6)
colnames(fst_mat) <- c("TB17","PWS17","SS17","BC17","WA17","CA17")
rownames(fst_mat) <- c("TB17","PWS17","SS17","BC17","WA17","CA17")

fst_mat[lower.tri(fst_mat, diag = F)]<-NA
write.csv(fst_mat, "Output/fst_pbs/MD7000/Fst_matrix_2017_all.csv")

# Melt the correlation matrix
melted_cormat <- melt(fst_mat, na.rm = TRUE)

melted_cormat[melted_cormat==0]<-NA
# Heatmap
melted_cormat$color<-"a"
melted_cormat$color[melted_cormat$value>=0.1]<-"b"

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colors=c("white", "blue"), limits=c(0, (max(melted_cormat$value, na.rm=T)+0.005)),na.value="gray80", 
                         name="Fst")+
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 12, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value, color=color),  size = 5)+
    scale_color_manual(values=c("black", "white"), guide='none')
    
ggsave("Output/fst_pbs/MD7000/Fst_matrix_2017_all.pdf", height = 6, width = 6)

## The pattern is somewhat consistent with isolation by distance, with CA and TB unique, PWS, SS, BC, WA are very genetically similar.


## 2D SFS
# The output from ANGSD is a flatten matrix: each value is the count of sites with the corresponding joint frequency ordered as
# [0,0] [0,1] [0,2] ..

#PWS17 = 56, PWS07=46
sfs2d<-read.table("Data/new_vcf/angsd/unfolded_PWS07_PWS17.sfs")
n1<-46 #pop1 # of individuals
n2<-56
sfs<-t(matrix(sfs2d, nrow=n2*2+1, ncol=n1*2+1))
sfsdf<-data.frame(lapply(data.frame(t(sfs)), unlist), stringsAsFactors=FALSE)


#plot 2DSFS as heatmap
#create a matrix from ANGSD output
vec2mat<-function(vec, n1,n2, pop1, pop2){
    n1<-n1
    n2<-n2
    pop1<-pop1
    pop2<-pop2
    ANGSD.2D.SFS <- scan(paste(vec, sep=""), quiet=T)
    ANGSD.2D.SFS <- t(matrix(ANGSD.2D.SFS, nrow=n2*2+1, ncol=n1*2+1))
    # mask non-variant sites
    ANGSD.2D.SFS[1,1] <- 0
    ANGSD.2D.SFS[nrow(ANGSD.2D.SFS),ncol(ANGSD.2D.SFS)] <- 0
    df<-data.frame(ANGSD.2D.SFS)
    colnames(df)<-0:(ncol(df)-1)
    df$count<-0:(nrow(df)-1)
        
    return(df)
}

#Plot heatmap     

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops.info$yr<-''
pops.info$yr[pops.info$year==96|pops.info$year==91]<-paste0(19,pops.info$year[pops.info$year==96|pops.info$year==91])
pops.info$yr[pops.info$year==07|pops.info$year==06|pops.info$year==17]<-paste0(20,pops.info$year[pops.info$year==07|pops.info$year==06|pops.info$year==17])
pops.info$yr<-apply(pops.info["yr"], 1, function(x) {if(x==206) x=2006
                                        if (x==207) x=2007
                                        else x=x})
pops.info$yr<-as.integer(pops.info$yr)

pops<-unique(pops.info$Population.Year)

pwss<-pops[grep("PWS",pops)]
tbs<-pops[grep("TB",pops)]
sss<-pops[grep("SS",pops)]
y17<-pops[grep("17",pops)]

comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-t(comb4)

#PWS
Plots<-list()
sfs.pws<-data.frame()
for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/unfolded_",pop1,"_",pop2,".sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
    if (pops.info$yr[pops_info$Population.Year==pop1][1]>pops.info$yr[pops_info$Population.Year==pop2][1]){
        sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
        colnames(sfs1)<-0:(ncol(sfs1)-1)
        sfs1$count<-0:(nrow(sfs1)-1)
        p2<-pop2
        pop2<-pop1
        pop1<-p2
    }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.pws<-rbind(sfs.pws, sfsm2)
    }

sfs.pws$pop1<-factor(sfs.pws$pop1, levels=c("PWS91","PWS96","PWS07","PWS17"))
sfs.pws$pop2<-factor(sfs.pws$pop2, levels=c("PWS91","PWS96","PWS07","PWS17"))

ggplot(sfs.pws, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_PWS.png", width = 10, height = 8, dpi=300)
   
#TB
sfs.tb<-data.frame()
for (i in 1: nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/unfolded_",pop1,"_",pop2,".sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
    if (pops.info$yr[pops_info$Population.Year==pop1][1]>pops.info$yr[pops_info$Population.Year==pop2][1]){
        sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
        colnames(sfs1)<-0:(ncol(sfs1)-1)
        sfs1$count<-0:(nrow(sfs1)-1)
        p2<-pop2
        pop2<-pop1
        pop1<-p2
    }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.tb<-rbind(sfs.tb, sfsm2)
}

sfs.tb$pop1<-factor(sfs.tb$pop1, levels=c("TB91","TB96","TB06","TB17"))
sfs.tb$pop2<-factor(sfs.tb$pop2, levels=c("TB91","TB96","TB06","TB17"))

ggplot(sfs.tb, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_TB.png", width = 10, height = 8, dpi=300)

#SS
sfs.ss<-data.frame()
for (i in 1: nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/unfolded_",pop1,"_",pop2,".sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
    if (pops.info$yr[pops_info$Population.Year==pop1][1]>pops.info$yr[pops_info$Population.Year==pop2][1]){
        sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
        colnames(sfs1)<-0:(ncol(sfs1)-1)
        sfs1$count<-0:(nrow(sfs1)-1)
        p2<-pop2
        pop2<-pop1
        pop1<-p2
    }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.ss<-rbind(sfs.ss, sfsm2)
}
sfs.ss$pop1<-factor(sfs.ss$pop1, levels=c("SS96","SS06","SS17"))
sfs.ss$pop2<-factor(sfs.ss$pop2, levels=c("SS96","SS06","SS17"))

ggplot(sfs.ss, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_SS.png", width = 10, height = 8, dpi=300)


##2017

Plots<-list()
sfs17<-data.frame()
for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/unfolded_",pop1,"_",pop2,".sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
   # if (pops.info$yr[pops_info$Population.Year==pop1][1]>pops.info$yr[pops_info$Population.Year==pop2][1]){
   #     sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
   #     colnames(sfs1)<-0:(ncol(sfs1)-1)
   #     sfs1$count<-0:(nrow(sfs1)-1)
   #     p2<-pop2
   #     pop2<-pop1
   #     pop1<-p2
   # }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs17<-rbind(sfs17, sfsm2)
}

#sfs.pws$pop1<-factor(sfs.pws$pop1, levels=c("PWS91","PWS96","PWS07","PWS17"))
sfs.pws$pop2<-factor(sfs.pws$pop2, levels=c("PWS91","PWS96","PWS07","PWS17"))

ggplot(sfs17, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_2017.png", width = 20, height = 16, dpi=300)







#####
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)

#Use same Y axis
sfs1D<-data.frame()
for (i in 1:length(pops)){
    sfs <- scan(paste("Data/sfs/downsample_41/unfolded/",pops[i],"_minQ20_minMQ30_unfolded.sfs",sep = ""))
    sfs1 <- data.frame(ac=sfs)
    sfs1$count<-0:(nrow(sfs1)-1)
    #remove the invariable sites
    sfs1<-sfs1[-c(1,nrow(sfs1)),]
    sfs1$pop<-pops[i]
    sfs1D<-rbind(sfs1D, sfs1)
}
    paste0(pops,",")
sfs1D$pop<-factor(sfs1D$pop, levels=c("PWS91","PWS96","PWS07","PWS17","TB91","TB96", "TB06","TB17","SS96","SS06", "SS17","BC17","WA17","CA17"))
ggplot(data=sfs1D, aes(x=count, y=ac))+
    facet_wrap(~pop, ncol=4)+
        geom_bar(stat="identity", color="gray")+xlab("Frequency bin")+ ylab("Number of alleles")+
        theme_classic()+
        scale_y_continuous(labels=scales::comma)+
        theme(strip.background = element_rect(
            color="black", fill="gray80", size=0.5, linetype="solid"))
ggsave("Output/SFS/1DSFS_all.pdf", width = 12, height = 8)


#######
# Compare Joe's downsampled SFS (capped to 100 read depth and 100M sites) vs. SFS based on all individuals and sites (capped 500 depths, >590M sites)

sfsj<-scan("Data/new_vcf/angsd/fromBam/PWS07.Joe.chr1.sfs")
sfsk<-scan("Data/new_vcf/angsd/fromBam/PWS07_unfolded_chr1.sfs")
#plot variable sites
barplot(sfsj[-c(1,length(sfsj))])
barplot(sfsk[-c(1,length(sfsk))])

ch1sfs<-data.frame(joe=sfsj, kaho=sfsk[1:length(sfsj)])
#standardized to Joe's'
p<-sfsj[1]/sfsk[1] #ratio of invariable sites
#0.957
p<-sfsj[2]/sfsk[2] #ratio of a singleton
#0.4529
ch1sfs$k_std<-ch1sfs$kaho*p
#remove the invariable sites
ch1<-ch1sfs[-c(1,nrow(ch1sfs)),]
ch1$count<-1:nrow(ch1)

ch1m<-melt(ch1, id.vars = "count")
ggplot(ch1m[ch1m$count>=1&ch1m$count<=40,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.5))

ch1m2<-ch1m[ch1m$variable!="kaho",]
ggplot(ch1m2[ch1m2$count>=1&ch1m2$count<=40,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))

p2<-sfsj[3]/sfsk[3] #ratio of a doubleton
#0.618
ch1$k_std2<-ch1$kaho*p2
#remove the invariable sites
ch12<-melt(ch1[,c("count","joe","k_std2")], id.vars="count")
ggplot(ch12[ch12$count>=1&ch12$count<=40,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))

#Shape is different.


#Compare with Joe's estimate for the whole genome
pw07_joe<-scan("Data/sfs/downsample_41/folded/PWS07_minQ20_minMQ30_folded.sfs")
#Compare with Joe's estimate for the whole genome
pw7<-pws07.sfs[1:length(pw07_joe),]
pw7$count<-1:nrow(pw7)
pw7<-pw7[,c("count","sum")]
pw7$joe<-pw07_joe

#standardized to Joe's'
p<-pw7$joe[1]/pw7$sum[1] #ratio of invariable sites
#0.637
p<-pw7$joe[2]/pw7$sum[2] #ratio of a singleton
#0.46489
pw7$k<-pw7$sum*p

#remove the invariable sites
pw7<-pw7[-c(1,nrow(pw7)),]
pw72<-pw7[,c("count","joe","k")]
pw7m<-melt(pw72, id.vars = "count")
ggplot(pw7m, aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.5))

ggplot(pw7m[pw7m$count>=1&pw7m$count<=20,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))



s<-scan('Data/sfs/downsample_41/folded/PWS07_minQ20_minMQ30_folded.sfs')
s<-s[-c(1,length(s))]
s<-s/sum(s)
barplot(s,names=1:length(s),main='SFS')




## Calculate thetas








#### annotation of potentially selected sites #####
# top 1% outliers of PBSs (Fst/PBS from ANGSD were based on windows)
pb.out<-pbsm[order(abs(pbsm$value), decreasing = T),] #210222 windows
pb.out<-pb.out[1:1913,]

pb.out<-merge(pb.out, pbs[,c("loc","region","midPos")], by="loc", all.x=T)
table(pb.out$variable)
#PBS96 PBS07 PBS17 
#268  1084   561  

length(unique(pb.out$loc))
#1913/2102

p17<-pb.out[pb.out$variable=="PBS17",]
p07<-pb.out[pb.out$variable=="PBS07",]
p96<-pb.out[pb.out$variable=="PBS96",]

sum<-data.frame(table(pb.out$variable,pb.out$chr))
#       1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
#PBS96 15  0 15 14  4 20 20 12  7 14  8 28  5 15 15  9  4 18 15  7  8 11 26 11  1  0
#PBS07 66 38 58 54 27 59 45 39 19 87 45 65 35 57 48 59 46 47 61 38 58 28 40 35 11 12
#PBS17 23 28 31 25 26 32 21 19  9 34 28 38 24 33  5 42 32 20 37 27 15  8 40  8  8 10

chr<-aggregate(sum$Freq, by=list(sum$Var2), sum)
chr[order(chr$x),]

#extract the regions and find the annotations
#create the bed file
pops<-c("PWS96","PWS07","PWS17")
files<-c("p96","p07","p17")
for (i in 1: length(pops)){
    loc<-data.frame()
    n=1
    pbs<-get(files[i])
    while (n <=nrow(pbs)){
        x<-pbs$loc[n]
        
        if (!is.na(pbs$loc[n+1])& pbs$loc[n+1]!=(x+1)) {
            newrow=c(pbs$loc[n], pbs$chr[n],pbs$midPos[n]-25000, pbs$midPos[n]+25000 )
            loc<-rbind(loc, newrow)
            n=n+1
        }
        else if (is.na(pbs$loc[n+1])) n=n+1
        else if (pbs$loc[n+1]==(x+1)){
            k=0
            while (!is.na(pbs$loc[n+k]) & pbs$loc[n+k]==(x+k)) k=k+1
            newrow=c(pbs$loc[n], pbs$ch[n],pbs$midPos[n]-25000, pbs$midPos[n+k-1]+25000)
            loc<-rbind(loc, newrow)
            n=n+k
        }
    }        
    
    loc<-loc[,-1]
    #convert the numbers to non-scientific
    loc[,2]<-as.integer(loc[,2])
    loc[,3]<-as.integer(loc[,3])
    write.table(loc, paste0("Output/fst_pbs/PWS/", pops[i],"_OutlierPBS_loci.bed"), quote=F, sep="\t",row.names=F, col.names = F)
}





# Create a new vcf file containing the loci in p07_loc.
 ## at terminal
vcftools --gzvcf Data/vcfs/population_PWS07_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS07_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS07_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS17_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS17_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS96_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS96_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS96_pbsoutlier --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS07_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS17_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS96_PBS_outlier.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS07_PBS_outlier.vcf > PWS07_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS17_PBS_outlier.vcf > PWS17_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS96_PBS_outlier.vcf > PWS96_annotation


#snpEff gene results
pops<-c("PWS96","PWS07","PWS17")
for (f in 1:3){
    #read the annotation info
    df<-read.table(paste0("Output/fst_pbs/PWS/",pops[f],"_annotation"), header = F)
    annotations<-data.frame()
    for (i in 1: nrow(df)){
        anns<-unlist(strsplit(df$V4[i], "\\|"))
        anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
        annotations<-rbind(annotations, anns)
    }     

    colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
    Ano<-cbind(df[,1:3], annotations)
    colnames(Ano)[1:3]<-c("chr","pos","AF")
    #remove the duplicated annotations for deeper digging
    remove<-!duplicated(annotations)
    Ano2<-Ano[remove,]
    write.csv(Ano2, paste0("Output/fst_pbs/PWS/", pops[f],"_highPBS_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-unique(geneids)
    
    
    sink(paste0("Output/fst_pbs/PWS/", pops[f],"_geneid_list.txt"))
    cat(paste0(geneids,"; "))
    sink(NULL)
}



annotations<-data.frame()
for (i in 1: nrow(df)){
    anns<-unlist(strsplit(df$V4[i], "\\|"))
    anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
    annotations<-rbind(annotations, anns)
}      

colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
Ano<-cbind(df[,1:3], annotations)
colnames(Ano)[1:3]<-c("chr","pos","AF")


#remove the duplicated annotations for deeper digging
remove<-!duplicated(annotations)
Ano2<-Ano[remove,]
write.csv(Ano2, "Output/fst_pbs/PWS/PWS07_highPBS_genelist.csv")            

geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
geneids<-unique(geneids)
geneids<-geneids[nchar(geneids)<=18]
geneids<-unique(geneids)


sink("Output/fst_pbs/PWS/PWS07_geneid_list.txt")
cat(paste0(geneids,"; "))
sink(NULL)



#####  from snpEff results ###
# Read the gene file

p07<-read.delim("Output/fst_pbs/PWS/PWS07.genes.txt",skip=1 )

#unique genes

#remove the no-annotated genes
p07a<-p07[grep("ENSCHA", p07$X.GeneName,invert = T),]
length(unique(p07a$X.GeneName)) #570 genes

write.table(unique(p07a$X.GeneName), "Output/fst_pbs/PWS/pws07_unique_genes.txt", quote = F, row.names = F)


## 
Ano2[Ano2$Gene_name=="FHIT",]
Ano2[Ano2$Gene_name=="styx",]

findGene<-function(gene){
    df<-Ano2[Ano2$Gene_name==gene,]
    
    return(df)
}

findGene("DDAH1")
findGene("vaspa")
findGene("hoxa2b")
Ano2[Ano2$Gene_name2=="hoxa2b",]
