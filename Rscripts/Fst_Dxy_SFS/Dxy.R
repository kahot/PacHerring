## Calculate Dxy
# ngsTools website (https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)
#"dxy been shown to be over-estimated and should be used only for inspecting the distribution and not to make inferences based on its absolute values"


library(gridExtra)
library(windowscanr)
source("Rscripts/BaseScripts.R")

#Use calcDxy.R (https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/calcDxy.R) using mafs from angsd

#create a script to run calcDxy.R for population pairs 
pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","TB17")
sss<-c("SS96","SS06","SS17")
y17<-c("TB17","PWS17","SS17","BC17","WA17","CA17")
comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-data.frame(t(comb4))
comb<-data.frame(rbind(comb1, comb2, comb3))
comb$Y2017<-"N"
comb4$Y2017<-"Y"
comb<-rbind(comb, comb4)


sink("Scripts/calculateDxy/run_calcDxy.sh")
cat("#!/bin/bash \n\n")
for (i in 1: nrow(comb)){
    cat(paste0("Rscript calcDxy.R  -p ",comb[i,1]," --popA='",comb[i,1],".mafs' -q ",comb[i,2], " --popB='",comb[i,2],".mafs'\n"))
    cat(paste0("mv Dxy_persite.txt Dxy_persite_",comb[i,1],"_",comb[i,2],".txt \n"))
}
sink(NULL)
#


#estimate summary statistics for two populations

### run ngstools in angsd
zcat Results/TSI.saf.gz > Results/TSI.saf
zcat Results/PEL.saf.gz > Results/PEL.saf
NSITES=`wc -l Data/intersect.txt | cut -f 1 -d " "` # if not already done
$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles Results/TSI.saf Results/PEL.saf -nsites $NSITES -nind 10 10 -outfile Results/TSI.PEL.stats.txt
####

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
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
comb4<-t(comb4)#PWS

for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))

        
    ### Manipulating the table and print dxy table
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # -> Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = c("dxy"), 
                       win_size = 50000,
                       win_step = 10000,
                       funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}
Plots1<-list()
for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    
    Plots1[[i]]<-ggplot(rol_win,  aes(x = win_mid, y = dxy)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.1)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_wrap(~ch, ncol = 9)
}

png("Output/Dxy//Dxy_PWS_md7000_line.png", height = 10, width = 22, res=150, units = "in")
grid.arrange(Plots1[[6]], Plots1[[3]], Plots1[[1]], Plots1[[2]],Plots1[[4]],Plots1[[5]], ncol=3)
dev.off()

#png("Output/Dxy//Dxy_PWS_md7000.png", height = 8, width = 18, res=150, units = "in")
#grid.arrange(Plots1[[6]], Plots1[[3]], Plots1[[1]], Plots1[[2]],Plots1[[4]],Plots1[[5]], ncol=3)
#dev.off()

for (i in 1: nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))
    
    
    ### Manipulating the table and print dxy table
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # -> Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = c("dxy"), 
                       win_size = 50000,
                       win_step = 10000,
                       funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}

for (i in 1: nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, groups = "chromo",position = "position", 
                       values = c("dxy"), win_size = 50000,win_step = 10000,funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}

for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, groups = "chromo",position = "position", 
                       values = c("dxy"), win_size = 50000,win_step = 10000,funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}




#Dxy comparison
pop1<-comb1[1,1]
pop2<-comb1[1,2]
rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
dxy<-rol_win[,c("chromo","ch","win_mid","dxy")]
colnames(dxy)[4]<-paste0(pop1,"_",pop2)

for (i in 2: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    dxy<-cbind(dxy, rol_win$dxy)
    colnames(dxy)[i+3]<-paste0(pop1,"_",pop2)
}

ggplot()+
    geom_point(data=dxy[dxy$ch==1,], aes(x=win_mid, y=PWS91_PWS96), color=blu,size=0.5, alpha=0.5)+
    geom_point(data=dxy[dxy$ch==1,], aes(x=win_mid, y=PWS07_PWS17), color=red,size=0.5, alpha=0.5)
    
ggplot()+
    geom_point(data=dxy[dxy$ch==1,], aes(x=win_mid, y=PWS91_PWS96), color=blu,size=0.5, alpha=0.5)+
    geom_point(data=dxy[dxy$ch==1,], aes(x=win_mid, y=PWS07_PWS96), color="red",size=0.5, alpha=0.5)

# TB
for (i in 1: nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))
    
    
    ### calculate dxy 
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # -> Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = c("dxy"), 
                       win_size = 50000,
                       win_step = 10000,
                       funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}
Plots2<-list()
for (i in 1: nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    
    rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    
    Plots2[[i]]<-ggplot(rol_win,  aes(x = win_mid, y = dxy)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.8)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_wrap(~ch, ncol = 9)
}
    
Plots2<-Plots
png("Output/Dxy//Dxy_TB_md7000_line.png", height = 10, width = 22, res=150, units = "in")
grid.arrange(Plots2[[6]], Plots2[[3]], Plots2[[1]], Plots2[[2]],Plots2[[4]],Plots2[[5]], ncol=3)
dev.off()

#SS
for (i in 1: nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))

    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = c("dxy"), 
                       win_size = 50000,
                       win_step = 10000,
                       funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}

Plots3<-list()
for (i in 1: nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    
    rol_win<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    
    Plots3[[i]]<-ggplot(rol_win,  aes(x = win_mid, y = dxy)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.1)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_wrap(~ch, ncol = 9)
}

png("Output/Dxy//Dxy_SS_md7000_line.png", height = 4, width = 18, res=150, units = "in")
grid.arrange(Plots3[[2]], Plots3[[3]], ncol=3)
dev.off()


#2017 pops (used files created by Joe ) -> slightly different 
Plots4<-list()
for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    
    rol_win<-read.delim(paste0("Data/dxy/windows/",pop1,"_",pop2,"_Dxy_50kb_win_10kb_step.txt"))
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    
    Plots4[[i]]<-ggplot(rol_win,  aes(x = win_mid, y = dxy)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.1)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_wrap(~ch, ncol = 9)
}
n=1:6
png("Output/Dxy//Dxy_Y2017_md7000_line1.png", height = 10, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()

n=7:12
png("Output/Dxy//Dxy_Y2017_md7000_line2.png", height = 10, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()

n=13:15
png("Output/Dxy//Dxy_Y2017_md7000_line3.png", height = 5, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()



for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/',pop1,'.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/',pop2,'.mafs'))
    
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = c("dxy"), 
                       win_size = 50000,
                       win_step = 10000,
                       funs = c("sum"))
    rol_win$dxy <- rol_win$dxy_sum/50000
    
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"))
}

Plots4<-list()
for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    
    rol_win<-read.delim(paste0("Data/dxy/windows/",pop1,"_",pop2,"_Dxy_50kb_win_10kb_step.txt"))
    rol_win$ch<-as.integer(gsub("chr", "", rol_win$chromo))
    
    Plots4[[i]]<-ggplot(rol_win,  aes(x = win_mid, y = dxy)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.1)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_wrap(~ch, ncol = 9)
}
n=1:6
png("Output/Dxy/Dxy_Y2017_md7000_line1.png", height = 10, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()

n=7:12
png("Output/Dxy/Dxy_Y2017_md7000_line2.png", height = 10, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()

n=13:15
png("Output/Dxy//Dxy_Y2017_md7000_line3.png", height = 5, width = 22, res=150, units = "in")
do.call(grid.arrange, c(Plots4[n], ncol=3))
dev.off()


