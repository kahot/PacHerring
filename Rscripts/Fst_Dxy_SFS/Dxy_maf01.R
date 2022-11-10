## Calculate Dxy

library(gridExtra)
library(windowscanr)
source("Rscripts/BaseScripts.R")


#estimate summary statistics for two populations

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)
pwss<-pops[grep("PWS",pops)]
tbs<-pops[grep("TB",pops)]
sss<-pops[grep("SS",pops)]
y17<-pops[grep("17",pops)]

comb1<-combn(pwss, 2)
comb1<-t(comb1)
#comb2<-combn(tbs, 2)
#comb2<-t(comb2)
#comb3<-combn(sss, 2)
#comb3<-t(comb3)
#comb4<-combn(y17, 2)
#comb4<-t(comb4)


#PWS
Plots<-list()
for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    allfreq1 <- read.delim(paste0('Data/new_vcf/AF/maf01/',pop1,'_maf01.mafs'))
    allfreq2 <- read.delim(paste0('Data/new_vcf/AF/maf01/',pop2,'_maf01.mafs'))

        
    ### Manipulating the table and print dxy table
    allfreq <- merge(allfreq1, allfreq2, by=c("chromo","position"))
    allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
    # -> Actual dxy calculation
    allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
    
    rol_win <- winScan(x = allfreq, 
                       groups = "chromo", 
                       position = "position", 
                       values = "dxy", 
                       win_size = 50000,
                       win_step = 10000,
                       funs = "sum")
    rol_win$dxy <- rol_win$dxy_sum/50000
    
    write.csv(rol_win, paste0("Output/Dxy/",pop1,"_",pop2,"_maf01_Dxy_50kwin_10kstep.csv"))
}


#compare Dxy estimated from maf05 and maf01 results
Plots<-list()
for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    dxy1<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_Dxy_50kwin_10kstep.csv"), row.names = 1)
    dxy2<-read.csv(paste0("Output/Dxy/",pop1,"_",pop2,"_maf01_Dxy_50kwin_10kstep.csv"), row.names = 1)
    dxy1$data<-"maf05"
    dxy2$data<-"maf01"
    dxy<-rbind(dxy1, dxy2)
    dxy$ch<-as.integer(gsub("chr", "", dxy$chromo))
    
    library(gridExtra)
    pdf("Output/Dxy/Dxy_maf01.vs.maf05_comparison_chr11-20.pdf", width = 15, height = 6)
    grid.arrange(
        ggplot(dxy[dxy$ch %in% c(11:15),],  aes(x = win_mid, y = dxy, color=data)) + 
        geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab("Dxy\n")+ xlab("")+ 
        ggtitle(paste(pop1,"vs. ", pop2))+
        geom_line( color = "steelblue", size=0.1)+
        #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
        facet_grid(data~ch),
        ggplot(dxy[dxy$ch %in% c(16:20),],  aes(x = win_mid, y = dxy, color=data)) + 
            geom_point(size = 1, color = "gray",alpha = 0.6, shape = 1)+
            theme_minimal()+
            theme(axis.text.x=element_blank())+
            ylab("Dxy\n")+ xlab("")+ 
            #ggtitle(paste(pop1,"vs. ", pop2))+
            geom_line( color = "steelblue", size=0.1)+
            #geom_smooth(method = "loess", se = FALSE, span = 1/100, color = "steelblue", size=0.5)+
            facet_grid(data~ch),nrow=2)
    dev.off()
    
    #ggsave("Output/Dxy/Dxy_maf01.vs.maf05_comparison_chr1-10.pdf", width = 15, height = 3)
    #ggsave("Output/Dxy/Dxy_maf01.vs.maf05_comparison_chr6-10.pdf", width = 15, height = 3)

}

#Two are highly similar but maf01 values are generally higher