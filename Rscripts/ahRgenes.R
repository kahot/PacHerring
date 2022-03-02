
library(ggplot2)
library(tidyverse)
library(plyranges)
library(ggpubr)
library(reshape2)

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


# Load data
#made smaller windows size than Joe's 1mb
win <- "5kb"
step <- "500b"
freqs_win_p <- read.csv(paste0("Output/freqs/PW_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
freqs_win_t <- read.csv(paste0("Output/freqs/TB_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
freqs_win_s <- read.csv(paste0("Output/freqs/SS_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)

win <- "50kb"
step <- "5kb"
freqs_win_p <- read.csv(paste0("Output/freqs/PW_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
freqs_win_t <- read.csv(paste0("Output/freqs/TB_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
freqs_win_s <- read.csv(paste0("Output/freqs/SS_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)



all_pops <- data.frame(chr = freqs_win_p$chr,
                       win_mid = freqs_win_p$win_mid,
                       zt01_p = freqs_win_p$zt01_mean,
                       zt02_p = freqs_win_p$zt02_mean,
                       zt03_p = freqs_win_p$zt03_mean,
                       zt12_p = freqs_win_p$zt12_mean,
                       zt13_p = freqs_win_p$zt13_mean,
                       zt23_p = freqs_win_p$zt23_mean,
                       zt01_t = freqs_win_t$zt01_mean,
                       zt02_t = freqs_win_t$zt02_mean,
                       zt03_t = freqs_win_t$zt03_mean,
                       zt12_t = freqs_win_t$zt12_mean,
                       zt13_t = freqs_win_t$zt13_mean,
                       zt23_t = freqs_win_t$zt23_mean,
                       zt12_s = freqs_win_s$zt12_mean,
                       zt13_s = freqs_win_s$zt13_mean,
                       zt23_s = freqs_win_s$zt23_mean)


all_pops$chr_num <- factor(all_pops$chr, levels = c(1:26))


# Look at shifts near ahr genes
gff <- read.delim("Data/annotations/Clupea_harengus.Ch_v2.0.2.100.gff3",header = FALSE, sep = "\t",comment.char = '#')
#head(gff)
colnames(gff) <- c("seqnames","source","type","start","end","score","strand","phase","attributes")
#gff$seqnames <- paste("chr" , as.character(gff$seqnames),sep = "")

genes <- gff[gff$type=="gene",]
genes <- genes %>% separate(attributes, c("at1","at2"), remove = FALSE,sep = "Name=") %>% separate(at2, c("gene_name","at3"), extra = "drop",sep = ";")
genes <- select(genes, -at1,-at3)
#g <- genes[genes$strand == "\\.",]
#head(genes)
genes <- droplevels(genes)


#genes in human Arylhydrocarbon Receptor (AhR) Signaling Pathway (https://maayanlab.cloud/Harmonizome/gene_set/Arylhydrocarbon+receptor+%28AhR%29+signaling+pathway%28Homo+sapiens%29/Wikipathways+Pathways)
ahr<-read.csv("Output/ahR/ahRsigPathway_geneset.csv")
ahrgenes<-tolower(ahr$Symbol)
#remove the 1 from cyp
ahrgenes[ahrgenes=="cyp1a1"]<-"cyp1a"
ahrgenes[ahrgenes=="cyp1b1"]<-"cyp1b"

ahrgenes<-c(ahrgenes, toupper(ahrgenes)) #add the upper capital names just in case
ahrgenes<-unique(ahrgenes)

geneset<-data.frame()
for (i in 1:length(ahrgenes)){
    df<-genes[grep(paste0('^',ahrgenes[i]), genes$gene_name),]
    geneset<-rbind(geneset, df)
}

#which chromosomes?
unique(geneset$seqnames)
# [1]  2 22 20  9 11  3  6  1 12 24 25 15 10 13 16

#remove the duplicates
geneset<-geneset[!duplicated(geneset),]




#Start with ahR 
i <- 2
#for (i in chrs){
all_pops1 <- all_pops[all_pops$chr_num %in% 2,]
ahr1b <- geneset[geneset$seqnames %in% 2,]

ggplot(all_pops1, aes_string(x = "win_mid", y = "zt01_p")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    facet_wrap(~chr_num, ncol = 9)+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt01_p", color = "blu"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt01_t", color = "red"),
                method = "loess", se = FALSE, span = 1/10)+
    scale_color_manual(values=c(red,blu),name = "population",
                       labels = c("PWS","TB"))+
    geom_vline(aes(xintercept = start), ahr1b)+
    geom_vline(aes(xintercept = end), ahr1b)+
    #ylim(-0.05, 0.05)+
    ggtitle("1991-1996")

#more focused to the region
ggplot(all_pops1, aes_string(x = "win_mid", y = "zt01_p")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    facet_wrap(~chr_num, ncol = 9)+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(aes(x = win_mid, y = zt01_p), color = blu)+
    geom_line(aes(x = win_mid, y = zt01_t), color = red)+
    scale_color_manual(values=c(red,blu),name = "population",
                       labels = c("PWS","TB"))+
    geom_vline(xintercept = ahr1b$start)+
    geom_vline(xintercept = ahr1b$end, color="yellow")+
    xlim(ahr1b$start-200000, ahr1b$start+200000)+
    ggtitle("1991-1996")


#PWS populations
colors <- c("1991-1996" =  blu, "1996-2007" = red, "2007-2017" = org)

ggplot(all_pops1, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    #theme(axis.text.x=element_blank())+
    #facet_wrap(~chr_num, ncol = 9)+
    geom_line(aes(x = win_mid, y = zt01_p, color="1991-1996"))+
    geom_line(aes(x = win_mid, y = zt12_p, color="1996-2007"))+
    geom_line(aes(x = win_mid, y = zt23_p, color="2007-2017"))+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahr1b$start, color="gray", size=0.5)+
    geom_vline(xintercept = ahr1b$end, color="gray", size=0.5)+
    xlim(ahr1b$start-50000, ahr1b$start+50000)+
    ggtitle("PWS")+
    annotate(geom="text", x=ahr1b$start+100, y=0.45, label="ahr1b",color ='black', size=4, hjust =0)
ggsave("Output/ahR/AHR1b_PWS_5k_AFshift.pdf", width=9.5, height =6)


ggplot(all_pops1, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_smooth(aes(x = win_mid, y = zt01_p, color = "1991-1996"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt12_p, color = "1996-2007"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt23_p, color = "2007-2017"),
                method = "loess", se = FALSE, span = 1/10)+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period") +
    scale_color_manual(values = colors[1:3])+
    geom_vline(xintercept = ahr1b$start, color="gray", size=0.5)+
    geom_vline(xintercept = ahr1b$end, color="gray", size=0.5)+
    xlim(ahr1b$start-100000, ahr1b$start+100000)+
    ggtitle("PWS")+
    annotate(geom="text", x=ahr1b$start+100, y=0.1, label="ahr1b",color ='black', size=4, hjust =0)
ggsave("Output/ahR/AHR1b_PWS_5k_AFshift_smooth.pdf", width=9.5, height =6)

### compare with other populations
all_pops2m<-melt(all_pops1[,c(2,3,6,8,9,12,14,15,17)], id.vars="win_mid")

colors <- c("1991-1996" =  blu, "1996-2007" = red, "2007-2017" = org)


#1. PWS vs. TB
all_pops2m$loc<-substr(all_pops2m$variable, 6,6)
all_pops2m$year<-substr(all_pops2m$variable, 1,4)

unique(all_pops2m$year)
unique(all_pops2m$loc)

year.labs<-c("1991-1996","1996-2007","2007-2017")
names(year.labs)<-c("zt01", "zt12","zt23")

#plots[[1]]<-
ggplot(all_pops2m[all_pops2m$loc %in% c("p","t"), ], aes(x = win_mid, y = value, color=year, linetype=loc)) + 
    theme_minimal()+
    facet_wrap(~year, ncol=1, labeller = labeller(year = year.labs), strip.position = "right")+
    geom_line()+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period", linetype="") +
    scale_color_manual(values = c(blu,red,org), labels= c("1991-1996","1996-2007","2007-2017") )+
    scale_linetype_manual(values = c("solid","dashed"), labels= c("PWS","TB", "SS"))+
    geom_vline(xintercept = ahr1b$start, color="gray", size=0.5)+
    geom_vline(xintercept = ahr1b$end, color="gray", size=0.5)+
    xlim(ahr1b$start-50000, ahr1b$start+50000)+
    ggtitle("AHR1b")
ggsave("Output/ahR/AHR1b_PWS_TB_shifts.pdf", width = 9, height = 8)

all_pops2m$loc<-factor(all_pops2m$loc, levels=c("p","t","s"))
ggplot(all_pops2m, aes(x = win_mid, y = value, color=year, linetype=loc)) + 
    theme_minimal()+
    facet_wrap(~year, ncol=1, labeller = labeller(year = year.labs), strip.position = "right")+
    geom_line()+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period", linetype="") +
    scale_color_manual(values = c(blu,red,org), labels= c("1991-1996","1996-2007","2007-2017") )+
    scale_linetype_manual(values = c("solid","dashed","dotted"), labels= c("PWS","TB", "SS"))+
    geom_vline(xintercept = ahr1b$start, color="gray", size=0.5)+
    geom_vline(xintercept = ahr1b$end, color="gray", size=0.5)+
    xlim(ahr1b$start-50000, ahr1b$start+50000)+
    ggtitle("AHR1b")
ggsave("Output/ahR/AHR1b_PWS_TB_SS_shifts.pdf", width = 9, height = 8)


ggplot(all_pops2m, aes(x = win_mid, y = value, color=year, linetype=loc)) + 
    theme_minimal()+
    facet_wrap(~year, ncol=1, labeller = labeller(year = year.labs), strip.position = "right")+
    geom_smooth(method = "loess", span = 1/10)+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period", linetype="") +
    scale_color_manual(values = c(blu,red,org), labels= c("1991-1996","1996-2007","2007-2017") )+
    scale_linetype_manual(values = c("solid","dashed","dotted"), labels= c("PWS","TB", "SS"))+
    geom_vline(xintercept = ahr1b$start, color="gray", size=0.5)+
    geom_vline(xintercept = ahr1b$end, color="gray", size=0.5)+
    xlim(ahr1b$start-50000, ahr1b$start+50000)+
    ggtitle("AHR1b")
ggsave("Output/ahR/AHR1b_PWS_TB_SS_shifts_smooth.pdf", width = 9, height = 8)





#ahrrb

all_pops22 <- all_pops[all_pops$chr_num %in% 22,]
ahrrb <- geneset[geneset$seqnames %in% 22,]

ggplot(all_pops22, aes_string(x = "win_mid", y = "zt01_p")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    facet_wrap(~chr_num, ncol = 9)+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt12_p", color = "blu"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_t", color = "red"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_s", color = "yel"),
                method = "loess", se = FALSE, span = 1/10)+
    scale_color_manual(values=c(red,blu,yel),name = "population",
                       labels = c("PWS","TB","SS"))
#ggsave("Output/ahR/CHR22_1996-2007_5k_AFshift_smooth.pdf", width=9.5, height =6)


ggplot(all_pops22, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_smooth(aes(x = win_mid, y = zt01_p, color = "1991-1996"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt12_p, color = "1996-2007"),
                method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt23_p, color = "2007-2017"),
                method = "loess", se = FALSE, span = 1/10)+
    labs(x = "",
         y = expression(paste(Delta, "z")),
         color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahrrb$start[1], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$end[1], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$start[3], color="gray70", size=0.5)+
    geom_vline(xintercept = ahrrb$end[3], color="gray70", size=0.5)+
    xlim(ahrrb$start[1]-100000, ahrrb$end[3]+100000)+
    ggtitle("PWS chr22")+
    annotate(geom="text", x=ahrrb$start[1]+100, y=0.065, label="ahrrb",color ='black', size=4, hjust =0)+
    annotate(geom="text", x=ahrrb$start[3]+100, y=0.065, label="aip",color ='black', size=4, hjust =0)
ggsave("Output/ahR/chr22_PWS_5k_AFshift_smooth.pdf", width=9.5, height =6)

#Focus on ahrrb only 5k no smooth
ggplot(all_pops22, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_line(aes(x = win_mid, y = zt01_p, color="1991-1996"))+
    geom_line(aes(x = win_mid, y = zt12_p, color="1996-2007"))+
    geom_line(aes(x = win_mid, y = zt23_p, color="2007-2017"))+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahrrb$start[1], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$end[1], color="gray", size=0.5)+
    xlim(ahrrb$start[1]-100000, ahrrb$end[1]+100000)+
    ggtitle("PWS chr22")+
    annotate(geom="text", x=ahrrb$start[1]+100, y=0.3, label="ahrrb",color ='black', size=4, hjust =0)
ggsave("Output/ahR/AHRRb_PWS_5k_AFshift.pdf", width=9.5, height =6)

ggplot(all_pops22, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_smooth(aes(x = win_mid, y = zt01_p, color = "1991-1996"), method = "loess", se = T, span = 1/15)+
    geom_smooth(aes(x = win_mid, y = zt12_p, color = "1996-2007"), method = "loess", se = T, span = 1/15)+
    geom_smooth(aes(x = win_mid, y = zt23_p, color = "2007-2017"), method = "loess", se = T, span = 1/15)+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahrrb$start[1], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$end[1], color="gray", size=0.5)+
    xlim(ahrrb$start[1]-100000, ahrrb$end[1]+100000)+
    ggtitle("PWS chr22")+
    annotate(geom="text", x=ahrrb$start[1]+100, y=0.3, label="ahrrb",color ='black', size=4, hjust =0)
ggsave("Output/ahR/AHRRb_PWS_5k_AFshift_smooth.pdf", width=9.5, height =6)

    

#Focus on AIP (Aryl hydrocarbon receptor interacting protein) only 5k no smooth
ggplot(all_pops22, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_line(aes(x = win_mid, y = zt01_p, color="1991-1996"))+
    geom_line(aes(x = win_mid, y = zt12_p, color="1996-2007"))+
    geom_line(aes(x = win_mid, y = zt23_p, color="2007-2017"))+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahrrb$start[3], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$end[3], color="gray", size=0.5)+
    xlim(ahrrb$start[3]-50000, ahrrb$end[3]+50000)+
    ggtitle("PWS chr22")+
    annotate(geom="text", x=ahrrb$start[3]+100, y=0.3, label="AIP", size=4, hjust =0)
ggsave("Output/ahR/AIP_PWS_5k_AFshift.pdf", width=9.5, height =6)

ggplot(all_pops22, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_smooth(aes(x = win_mid, y = zt01_p, color = "1991-1996"), method = "loess", se = T, span = 1/15)+
    geom_smooth(aes(x = win_mid, y = zt12_p, color = "1996-2007"), method = "loess", se = T, span = 1/15)+
    geom_smooth(aes(x = win_mid, y = zt23_p, color = "2007-2017"), method = "loess", se = T, span = 1/15)+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = ahrrb$start[3], color="gray", size=0.5)+
    geom_vline(xintercept = ahrrb$end[3], color="gray", size=0.5)+
    xlim(ahrrb$start[3]-50000, ahrrb$end[3]+50000)+
    ggtitle("PWS chr22")+
    annotate(geom="text", x=ahrrb$start[3]+100, y=0.3, label="AIP", size=4, hjust =0)
ggsave("Output/ahR/AIP_PWS_5k_AFshift_smooth.pdf", width=9.5, height =6)

#CYP1A
all_pops6 <- all_pops[all_pops$chr_num %in% 6,]
cyp <- geneset[geneset$seqnames %in% 6,]

ggplot(all_pops6, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_line(aes(x = win_mid, y = zt01_p, color="1991-1996"))+
    geom_line(aes(x = win_mid, y = zt12_p, color="1996-2007"))+
    geom_line(aes(x = win_mid, y = zt23_p, color="2007-2017"))+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = cyp$start[2], color="gray", size=0.5)+
    geom_vline(xintercept = cyp$end[2], color="gray", size=0.5)+
    xlim(cyp$start[2]-100000, cyp$end[2]+100000)+
    ggtitle("PWS chr6")+
    annotate(geom="text", x=ahrrb$start[3]+100, y=0.3, label="CYP1A", size=4, hjust =0)
ggsave("Output/ahR/CYP1A_PWS_5k_AFshift.pdf", width=9.5, height =6)
#no snps are in CYP1A 

ggplot(all_pops6, aes(x = win_mid, y = zt01_p,)) + 
    theme_minimal()+
    geom_smooth(aes(x = win_mid, y = zt01_p, color = "1991-1996"), method = "loess", se = T, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt12_p, color = "1996-2007"), method = "loess", se = T, span = 1/10)+
    geom_smooth(aes(x = win_mid, y = zt23_p, color = "2007-2017"), method = "loess", se = T, span = 1/10)+
    labs(x = "", y = expression(paste(Delta, "z")), color = "Period") +
    scale_color_manual(values = colors)+
    geom_vline(xintercept = cyp$start[2], color="gray", size=0.5)+
    geom_vline(xintercept = cyp$end[2], color="gray", size=0.5)+
    xlim(cyp$start[2]-100000, cyp$end[2]+100000)+
    ggtitle("PWS chr6")+
    annotate(geom="text", x=ahrrb$start[3]+100, y=0.3, label="CYP1A", size=4, hjust =0)
ggsave("Output/ahR/AIP_PWS_5k_AFshift_smooth.pdf", width=9.5, height =6)

#why the big dip?
cp<-all_pops6[all_pops6$win_mid>(cyp$start[2]-50000) & all_pops6$win_mid<(cyp$end[2]+50000),] #no snps are in CYP1A 





### create the figure for all genes in geneset 

colors <- c("1991-1996" =  blu, "1996-2007" = red, "2007-2017" = org)
year.labs<-c("1991-1996","1996-2007","2007-2017")
names(year.labs)<-c("zt01", "zt12","zt23")

for (i in 1: nrow(geneset)){
    print(i)
    chr<-geneset$seqnames[i]
    
    dt <- all_pops[all_pops$chr_num %in% chr,]
    ginfo <- geneset[i,]
    
    #melt the dt
    dtm<-melt(dt[,c(2,3,6,8,9,12,14,15,17)], id.vars="win_mid")
    
    #extract location and year info
    dtm$loc<- substr(dtm$variable, 6,6)
    dtm$year<-substr(dtm$variable, 1,4)
    
    dtm$loc<-factor(dtm$loc, levels=c("p","t","s"))
    
    ggplot(dtm[dtm$win_mid>=(ginfo$start-100000)& dtm$win_mid<=(ginfo$end+100000),], aes(x = win_mid, y = value, color=year, linetype=loc)) + 
        theme_minimal()+
        facet_wrap(~year, ncol=1, labeller = labeller(year = year.labs), strip.position = "right")+
        geom_line()+
        labs(x = "",y = expression(paste(Delta, "z")), color = "Period", linetype="") +
        scale_color_manual(values = c(blu,red,org), labels= c("1991-1996","1996-2007","2007-2017") )+
        scale_linetype_manual(values = c("solid","dashed","dotted"), labels= c("PWS","TB", "SS"))+
        geom_vline(xintercept = ginfo$start, color="gray", size=0.3)+
        geom_vline(xintercept = ginfo$end, color="gray", size=0.3)+
        #xlim(ginfo$start-100000, ginfo$end+100000)+
        ggtitle(ginfo$gene_name)
    ggsave(paste0("Output/ahR/FIg_5k/",ginfo$gene_name, "_shifts.pdf"), width = 9, height = 8)
    
    
    ggplot(dtm[dtm$win_mid>=(ginfo$start-100000)& dtm$win_mid<=(ginfo$end+100000),], aes(x = win_mid, y = value, color=year, linetype=loc)) + 
        theme_minimal()+
        facet_wrap(~year, ncol=1, labeller = labeller(year = year.labs), strip.position = "right")+
        geom_smooth(method = "loess", span = 0.4)+
        labs(x = "",
             y = expression(paste(Delta, "z")),
             color = "Period", linetype="") +
        scale_color_manual(values = c(blu,red,org), labels= c("1991-1996","1996-2007","2007-2017") )+
        scale_linetype_manual(values = c("solid","dashed","dotted"), labels= c("PWS","TB", "SS"))+
        geom_vline(xintercept = ginfo$start, color="gray", size=0.3)+
        geom_vline(xintercept = ginfo$end, color="gray", size=0.3)+
        ggtitle(paste0("chr",chr,": ",ginfo$gene_name))
    ggsave(paste0("Output/ahR/FIg_5k/",ginfo$gene_name, "_shifts_smooth.pdf"), width = 9, height = 8)

}
