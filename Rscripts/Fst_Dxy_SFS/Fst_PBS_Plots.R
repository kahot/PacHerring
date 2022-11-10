#Pacific Herring Fst all comparisons chromosome view 
#Visualize Fst and population branch statistic across all comparisons for each chromosome.

library(ggplot2)
library(tidyverse)
#library(reticulate)
library(reshape2)
library(plyranges)
library(seqinr)
library(gridExtra)
#BiocManager::install("visFuns.R")


# color-blind friendly 
# Wong, B. Points of view: Color blindness. Nat Methods (2011).
bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'


# A population’s population branch statistic value at a given locus corresponds to the magnitude of allele frequency change relative to its divergence from the other two populations.

pop_names = c("PWS91","PWS96","PWS07","PWS17","TB91","TB96","TB06","TB17","SS96","SS06","SS17","BC17","WA17","CA17")


#fst/pbs calculated for 50k windows by Joe
fsts <- list.files(path="Data/fst_pbs/maf05/")
fsts
i=8
fsts<-list.files(path="Output/fst_pbs/", pattern="Fst_pbs_50kb")
for (i in 1: length(fsts)){  

    pops <- strsplit(fsts[i], "_")
    pops <- c(pops[[1]][8], pops[[1]][9],gsub("\\.txt","",pops[[1]][10]))
    cat ("Comparison between ",pops,"\n")
    fst <-read.delim(paste0("Data/fst_pbs/maf05/",fsts[i]))
    
    fst <- fst %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
    #fst$chr_num <- factor(fst$chr, levels = c(1:26))
    
    # Fst with smooth line
    p1 <- ggplot(fst, aes(x = midPos, y = Fst01)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[1],"vs. ", pops[2]))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
        facet_wrap(~chr, ncol = 9)
    # Fst with actual line to highlight the differences
    p11 <- ggplot(fst, aes(x = midPos, y = Fst01)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[1],"vs. ", pops[2]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    p2 <- ggplot(fst, aes(x = midPos, y = Fst02)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[1],"vs. ", pops[3]))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
        facet_wrap(~chr, ncol = 9)
    p22 <- ggplot(fst, aes(x = midPos, y = Fst02)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[1],"vs. ", pops[3]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    
    p3 <- ggplot(fst, aes(x = midPos, y = Fst12)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[2],"vs. ", pops[3]))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
        facet_wrap(~chr, ncol = 9)
    p33 <- ggplot(fst, aes(x = midPos, y = Fst12)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Fst\n")+ xlab("")+ 
        ggtitle(paste(pops[2],"vs. ", pops[3]))+
        geom_line(color=blu, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    pdf(paste0("Output/fst_pbs/Fst_compare_", pops[1],"_",pops[2],"_",pops[3],".pdf"), height = 12, width = 8)
    grid.arrange(p1,p2,p3, ncol=1)
    dev.off()
    pdf(paste0("Output/fst_pbs/Fst_compare2_", pops[1],"_",pops[2],"_",pops[3],".pdf"), height = 12, width = 8)
    grid.arrange(p11,p22,p33, ncol=1)
    dev.off()
    
    
    #PBS
    p4 <- ggplot(fst, aes(x = midPos, y = PBS0)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+ 
        ggtitle(paste0(pops[1]," (",pops[2],' & ',pops[3],')'))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = red)+
        facet_wrap(~chr, ncol = 9)
    
    p44 <-ggplot(fst, aes(x = midPos, y = PBS0)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+ 
        ggtitle(paste0(pops[1]," (",pops[2],' & ',pops[3],')'))+
        geom_line(color=red, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    p5 <- ggplot(fst, aes(x = midPos, y = PBS1)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+
        ggtitle(paste0(pops[2]," (",pops[1],' & ',pops[3],')'))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = red)+
        facet_wrap(~chr, ncol = 9)
    
    p55 <- ggplot(fst, aes(x = midPos, y = PBS1)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+ 
        ggtitle(paste0(pops[2]," (",pops[1],' & ',pops[3],')'))+
        geom_line(color=red, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    p6 <- ggplot(fst, aes(x = midPos, y = PBS2)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+ 
        ggtitle(paste0(pops[3]," (",pops[1],' & ',pops[2],')'))+
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = red)+
        facet_wrap(~chr, ncol = 9)

    p66<- ggplot(fst, aes(x = midPos, y = PBS2)) + 
        geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("PBS\n")+ xlab("")+ 
        ggtitle(paste0(pops[3]," (",pops[1],' & ',pops[2],')'))+
        geom_line(color=red, size=0.2)+
        facet_wrap(~chr, ncol = 9)
    
    pdf(paste0("Output/fst_pbs/PBS_compare_", pops[1],"_",pops[2],"_",pops[3],".pdf"), height = 12, width = 8)
    grid.arrange(p4,p5,p6, ncol=1)
    dev.off()
    
    pdf(paste0("Output/fst_pbs/PBS_compare2_", pops[1],"_",pops[2],"_",pops[3],".pdf"), height = 12, width = 8)
    grid.arrange(p44,p55,p66, ncol=1)
    dev.off()
    
    
}



#MD7000 mean Fst
##########
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







##########
## Fst matrix temporal samples (within population comparisons)
# temporal samples per populations
pops<-c("PWS","SS","TB")

for (i in 1:length(pops)){
    p<-pops[i]
    
    if (p=="PWS") y3="07"
    else y3="06"
    if (p=="SS"){
        fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_SS96_SS06_SS17.txt")
        p96xp07 <- round(mean(fst$Fst01),4)
        p96xp17 <- round(mean(fst$Fst02),4)
        p07xp17 <- round(mean(fst$Fst12),4)
        
        fst_vec <- c(0,p96xp07,p96xp17,
                     p96xp07,0,p07xp17,
                     p96xp17,p07xp17,0)
        fst_mat = matrix(fst_vec, nrow = 3, ncol = 3)
        colnames(fst_mat) <- c("SS96","SS06","SS17")
        rownames(fst_mat) <- c("SS96","SS06","SS17")
        
    }
    else{
        fst<-read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p, y3,".txt"))
        p91xp96 <- round(mean(fst$Fst01),4)
        p91xp07 <- round(mean(fst$Fst02),4)
        p96xp07 <- round(mean(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p,"17.txt"))
        p91xp17 <- round(mean(fst$Fst02),4)
        p96xp17 <- round(mean(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,y3,"_",p,"17.txt"))
        p07xp17 <- round(mean(fst$Fst12),4)
        fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
                 p91xp96,0,p96xp07,p96xp17,
                 p91xp07,p96xp07,0,p07xp17,
                 p91xp17,p96xp17,p07xp17,0)
    
        fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
        colnames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
        rownames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
    }
    
    fst_mat[lower.tri(fst_mat)]<- NA
    write.csv(fst_mat,paste0("Output/fst_pbs/",p,"_Fst_matirx.csv"))
}
    
##calculate SE
library(plotrix)
for (i in 1:length(pops)){
    p<-pops[i]
    
    if (p=="PWS") y3="07"
    else y3="06"
    if (p=="SS"){
        fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_SS96_SS06_SS17.txt")
        p96xp07 <- round(mean(fst$Fst01),4)
        p96xp17 <- round(mean(fst$Fst02),4)
        p07xp17 <- round(mean(fst$Fst12),4)
        
        fst_vec <- c(0,p96xp07,p96xp17,
                     p96xp07,0,p07xp17,
                     p96xp17,p07xp17,0)
        fst_mat = matrix(fst_vec, nrow = 3, ncol = 3)
        colnames(fst_mat) <- c("SS96","SS06","SS17")
        rownames(fst_mat) <- c("SS96","SS06","SS17")
        
    }
    else{
        fst<-read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p, y3,".txt"))
        p91xp96 <- std.error(fst$Fst01)
        p91xp07 <- round(mean(fst$Fst02),4)
        p96xp07 <- round(mean(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p,"17.txt"))
        p91xp17 <- round(mean(fst$Fst02),4)
        p96xp17 <- round(mean(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,y3,"_",p,"17.txt"))
        p07xp17 <- round(mean(fst$Fst12),4)
        fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
                     p91xp96,0,p96xp07,p96xp17,
                     p91xp07,p96xp07,0,p07xp17,
                     p91xp17,p96xp17,p07xp17,0)
        
        fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
        colnames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
        rownames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
    }
    
    fst_mat[lower.tri(fst_mat)]<- NA
    write.csv(fst_mat,paste0("Output/fst_pbs/",p,"_Fst_matirx.csv"))
}



##  Use median 
pops<-c("PWS","SS","TB")
for (i in 1:length(pops)){
    p<-pops[i]
    
    if (p=="PWS") y3="07"
    else y3="06"
    if (p=="SS"){
        fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_SS96_SS06_SS17.txt")
        p96xp07 <- round(median(fst$Fst01),4)
        p96xp17 <- round(median(fst$Fst02),4)
        p07xp17 <- round(median(fst$Fst12),4)
        
        fst_vec <- c(0,p96xp07,p96xp17,
                     p96xp07,0,p07xp17,
                     p96xp17,p07xp17,0)
        fst_mat = matrix(fst_vec, nrow = 3, ncol = 3)
        colnames(fst_mat) <- c("SS96","SS06","SS17")
        rownames(fst_mat) <- c("SS96","SS06","SS17")
        
    }
    else{
        fst<-read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p, y3,".txt"))
        p91xp96 <- round(median(fst$Fst01),4)
        p91xp07 <- round(median(fst$Fst02),4)
        p96xp07 <- round(median(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p,"17.txt"))
        p91xp17 <- round(median(fst$Fst02),4)
        p96xp17 <- round(median(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,y3,"_",p,"17.txt"))
        p07xp17 <- round(median(fst$Fst12),4)
        fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
                     p91xp96,0,p96xp07,p96xp17,
                     p91xp07,p96xp07,0,p07xp17,
                     p91xp17,p96xp17,p07xp17,0)
        
        fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
        colnames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
        rownames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
    }
    
    fst_mat[lower.tri(fst_mat)]<- NA
    write.csv(fst_mat,paste0("Output/fst_pbs/",p,"_Fst_median_matirx.csv"))
}








#fstmats<-list.files("Output/fst_pbs/", pattern="_Fst_matirx.csv")
#fstmats<-list.files("Output/fst_pbs/", pattern="_Fst_median_matirx.csv")

fstmats<-list.files("Output/fst_pbs/MD7000/", pattern="Fst_matirx.csv")


for (i in 1:3){
    df<-read.csv(paste0("Output/fst_pbs/MD7000/",fstmats[i]))

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
    ggsave(paste0("Output/fst_pbs/MD7000/",fname,"2.pdf"), width = 5, height = 5)
    
}

#Create a data table containing all Fst points 
pops<-c("PWS","SS","TB")
for (i in 1:length(pops)){
    p<-pops[i]
    
    if (p=="PWS") y3="07"
    else y3="06"
    if (p=="SS"){
        fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_SS96_SS06_SS17.txt")
        p96xp07 <- round(median(fst$Fst01),4)
        p96xp17 <- round(median(fst$Fst02),4)
        p07xp17 <- round(median(fst$Fst12),4)
        
        fst_vec <- c(0,p96xp07,p96xp17,
                     p96xp07,0,p07xp17,
                     p96xp17,p07xp17,0)
        fst_mat = matrix(fst_vec, nrow = 3, ncol = 3)
        colnames(fst_mat) <- c("SS96","SS06","SS17")
        rownames(fst_mat) <- c("SS96","SS06","SS17")
        
    }
    else{
        fst<-read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p, y3,".txt"))
        pws <- data.frame(fst$Fst01)
        pws$yr1<-1991
        pws$Year <-1991
        
        p91xp07 <- round(median(fst$Fst02),4)
        p96xp07 <- round(median(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,"96_",p,"17.txt"))
        p91xp17 <- round(median(fst$Fst02),4)
        p96xp17 <- round(median(fst$Fst12),4)
        fst <- read.delim(paste0("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_",p,"91_",p,y3,"_",p,"17.txt"))
        p07xp17 <- round(median(fst$Fst12),4)
        fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
                     p91xp96,0,p96xp07,p96xp17,
                     p91xp07,p96xp07,0,p07xp17,
                     p91xp17,p96xp17,p07xp17,0)
        
        fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
        colnames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
        rownames(fst_mat) <- c(paste0(p,"91"),paste0(p,"96"),paste0(p,y3),paste0(p,"17"))
    }
    
    fst_mat[lower.tri(fst_mat)]<- NA
    write.csv(fst_mat,paste0("Output/fst_pbs/",p,"_Fst_median_matirx.csv"))
}





fs<-read.csv("Output/fst_pbs/Fst_values_3pop.csv")
fs$X[fs$X==2007]<-2006
fs$Year[fs$Year==2007]<-2006


ggplot()+
    geom_point(data=fs[fs$X==1991,], aes(x=Year, y=Fst, color=pop), size=2)+
    geom_line(data=fs[fs$X==1991,], aes(x=Year, y=Fst, color=pop))+
    scale_color_manual(values=c(red,blu))+
    xlab('')+ylab("Fst")+
    theme_classic()+theme(legend.title = element_blank())
ggsvae()

ggplot()+
    geom_point(data=fs[fs$X==1991,], aes(x=Year, y=Fst, color=pop), size=2)+
    geom_line(data=fs[fs$X==1991,], aes(x=Year, y=Fst, color=pop))+
    geom_point(data=fs[fs$pop=="SS",], aes(x=Year, y=Fst), color=yel, size=2)+
    scale_color_manual(values=c(red,blu))
xlab('')+ylab("Fst")+
    theme_classic()


    



## 2017 between populations Fst matrix
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB17_PWS17_SS17.txt")
tbxpw <- round(mean(fst$Fst01),4)
tbxss <- round(mean(fst$Fst02),4)
pwxss <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB17_PWS17_BC17.txt")
tbxbc <- round(mean(fst$Fst02),4)
pwxbc <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB17_PWS17_WA17.txt")
tbxwa <- round(mean(fst$Fst02),4)
pwxwa <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB17_PWS17_CA17.txt")
tbxca <- round(mean(fst$Fst02),4)
pwxca <- round(mean(fst$Fst12),4)

fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS17_SS17_BC17.txt")
ssxbc <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS17_SS17_WA17.txt")
ssxwa <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS17_SS17_CA17.txt")
ssxca <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_BC17_WA17_CA17.txt")
bcxwa <- round(mean(fst$Fst01),4)
bcxca <- round(mean(fst$Fst02),4)
waxca <- round(mean(fst$Fst12),4)

fst_vec <- c(0,tbxpw,tbxss,tbxbc,tbxwa,tbxca,
             tbxpw,0,pwxss,pwxbc,pwxwa,pwxca,
             tbxss,pwxss,0,ssxbc,ssxwa,ssxca,
             tbxbc,pwxbc,ssxbc,0,bcxwa,bcxca,
             tbxwa,pwxwa,ssxwa,bcxwa,0,waxca,
             tbxca,pwxca,ssxca,bcxca,waxca,0)

fst_mat = matrix(fst_vec, nrow = 6, ncol = 6)
colnames(fst_mat) <- c("TB17","PWS17","SS17","BC17","WA17","CA17")
rownames(fst_mat) <- c("TB17","PWS17","SS17","BC17","WA17","CA17")

fst_mat[lower.tri(fst_mat, diag = F)]<-NA
write.csv(fst_mat, "Output/fst_pbs/Fst_matrix_2017_all.csv")

# Melt the correlation matrix
melted_cormat <- melt(fst_mat, na.rm = TRUE)

melted_cormat[melted_cormat==0]<-NA
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colors=c("white", "blue"), limits=c(0, (max(melted_cormat$value, na.rm=T)+0.005)),na.value="gray80", 
                        name="Fst")+
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 12, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 5)
ggsave("Output/fst_pbs/Fst_matrix_2017_all.pdf", height = 6, width = 6)
