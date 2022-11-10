
library(ggplot2)
library(tidyverse)
library(reticulate)
library(reshape2)
library(plyranges)
library(seqinr)
library(windowscanr)
library(ggpubr)
library(venn)
library(gridExtra)

bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'

pops<-c("PW","SS","TB")

files<-list.files("Data/freqs/persite/")

#calculate smaller window size
#window_size <- 100000
#step_size <- 10000

window_size <- 50000
step_size <- 5000



#change the window size
for (i in 1:length(files)){
    print(i)
    df<-read.table(paste0("Data/freqs/persite/", files[i]))
    
    if (i ==1|i==3) vals<-c("zt01","zt02","zt03","zt12","zt13","zt23")
    if (i==2) vals<-c("zt12","zt13","zt23")
    freqs_win <- winScan(x = df, 
                         groups = "chr", 
                         position = "pos", 
                         values = vals, 
                         win_size = window_size,
                         win_step = step_size,
                         funs = c("mean"),
                         cores=10)
    freqs_win <- na.omit(freqs_win)
    freqs_win <- freqs_win %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    
    #write.csv(freqs_win, paste0("Output/freqs/",substr(files[i],1,2), "_shifts_100kb_win_10kb_step.csv"))
    write.csv(freqs_win, paste0("Output/freqs/",substr(files[i],1,2), "_shifts_50kb_win_5kb_step.csv"))
}


#AF=Allele Frequency 
#AF is calculated as AC/AN, where AC is total ALT allele count and AN is total allele called in genotypes


# Compare magnitude of AF shifts between populations

# the distribution of ∆z for each population from
#- 1991 --> 1996
#- 1996 --> 2007
#- 2007 --> 2017


### different parameters
win <- "100kb"
step <- "10kb"

#for (i in 1:3){
#    pop<-pops[i]
#    df<-read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
#}


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

write.csv(all_pops,"Output/freqs/freq_shifts_all_populations_100kb_10kb.csv")

#max and minimum delta 
apply(all_pops[,c(3,6,8,9,12,14,15,17)],2,FUN =  max) #0.4315
apply(all_pops[,c(3,6,8,9,12,14,15,17)],2,FUN =  min) #-0.4397
#1991-1996 change


zp <- data.frame(z91_96 = all_pops$zt01_p,pop = "PWS")
zt <- data.frame(z91_96 = all_pops$zt01_t,pop = "TB")
z1 <- rbind(zp,zt)

p1<-ggplot(z1, aes(x=z91_96, fill=pop)) +
    geom_density(alpha=0.4)+
    scale_fill_manual(name = "Population",values=c(red,blu))+
    ylab("density\n") + xlab(expression(paste(Delta, "z 1991-1996")))+
    xlim(-0.45,0.45)+
    theme_minimal()
p1

zp <- data.frame(z96_07 = all_pops$zt12_p,pop = "PWS")
zt <- data.frame(z96_07 = all_pops$zt12_t,pop = "TB")
zs <- data.frame(z96_07 = all_pops$zt12_s,pop = "SS")

z2 <- rbind(zp,zt,zs)

p2<-ggplot(z2, aes(x=z96_07, fill=pop)) +
    geom_density(alpha=0.4)+
    scale_fill_manual(name = "Population",values=c(red,yel,blu))+
    ylab("density\n") + xlab(expression(paste(Delta, "z 1996-2007")))+
    xlim(-0.45,0.45)+
    theme_minimal()
p2

zp <- data.frame(z07_17 = all_pops$zt23_p,pop = "PWS")
zt <- data.frame(z07_17 = all_pops$zt23_t,pop = "TB")
zs <- data.frame(z07_17 = all_pops$zt23_s,pop = "SS")

z3 <- rbind(zp,zt,zs)

p3<-ggplot(z3, aes(x=z07_17, fill=pop)) +
    geom_density(alpha=0.4)+
    scale_fill_manual(name = "Population",values=c(red,yel,blu))+
    ylab("density\n") + xlab(expression(paste(Delta, "z 2007-2017")))+
    xlim(-0.45,0.45)+
    theme_minimal()
p3

pdf("Output/freqs/deltaZ_density_100kb_10kb.pdf", height = 8,width = )
grid.arrange(p1,p2,p3, ncol=1)
dev.off()



#################################
# Genome-wide allele frequency shifts

## Within population all chromosomes
win="100kb"
step="10kb"

plot_af_chromosomes <- function(pop,years,win,step){
    freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    
    yr<-substr(years,3,4)
    if (yr=="01") period="1991-1996"
    if (yr=="12") period="1996-2007"
    if (yr=="23") period="2007-2017"
    
    p <- ggplot(freqs_win, aes_string(x = "win_mid", y = years)) + 
        geom_point(size = 1, color = gry,alpha = 0.7, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = red)+
        facet_wrap(~chr_num, ncol = 9)+
        ggtitle(paste(pop,period, sep = " "))
    print(p)
    ggsave(plot=p,paste0("Output/freqs/",pop,"_AF-shifts_perChromosome_window_",win,"_year_",period,".pdf"),width = 10, height = 8 )
    
}


plot_af_chromosomes_within <- function(pop,win,step){
    freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
        geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
        geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
        facet_wrap(~chr_num, ncol = 9)+
        scale_color_manual(values=c(lir,grb,org),name = "contrast",
                           labels = c("1991-1996", "1996-2007", "2007-2017"))+
        ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
        ggtitle(pop)
    print(p)
    ggsave(plot=p,paste0("Output/freqs/",pop,"_AF-shifts_perChromosome_window_",win,".pdf"),width = 10, height = 8 )
    
}


plot_af_chromosomes_within2 <- function(pop,win,step){
    freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

    p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt12_mean")) + 
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
        geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
        facet_wrap(~chr_num, ncol = 9)+
        scale_color_manual(values=c(grb,org),name = "contrast",
                           labels = c("1996-2007", "2007-2017"))+
        ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
        ggtitle(pop)
    print(p)
    ggsave(plot=p,paste0("Output/freqs/",pop,"_AF-shifts_perChromosome_window_",win,".pdf"),width = 10, height = 8 )
}

# Within Population
for (i in 1:length(pops)){
    pop<-pops[i]
    if (i ==2) {
        plot_af_chromosomes_within2(pop,win,step)
        plot_af_chromosomes(pop,"zt12_mean",win,step)
        plot_af_chromosomes(pop,"zt23_mean",win,step)
    }
    else {
        plot_af_chromosomes_within(pop,win,step)
        plot_af_chromosomes(pop,"zt01_mean",win,step)
        plot_af_chromosomes(pop,"zt12_mean",win,step)
        plot_af_chromosomes(pop,"zt23_mean",win,step)
    }
}


pop<-"PW"
freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

freqs<-melt(freqs_win[,c("chr_num","win_mid","zt01_mean" ,"zt12_mean","zt12_mean")], id.vars=c("chr_num", "win_mid"))

ggplot(freqs, aes(x = win_mid, y =value, color=variable)) +
    facet_wrap(~Group, ncol=2)

ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
    geom_line(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
    geom_line(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
    facet_wrap(~chr_num, ncol = 9)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
    ggtitle(pop)



    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line()
    geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
    facet_wrap(~chr_num, ncol = 9)+
    scale_color_manual(values=c(lir,grb,org),name = "contrast",
                       labels = c("1991-1996", "1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+theme(legend.position = c(0.95, 0.2))+
    ggtitle(pop)
print(p)
ggsave(plot=p,paste0("Output/freqs/",pop,"_AF-shifts_perChromosome_window_",win,".pdf"),width = 10, height = 8 )






## Within population chromosomes with putative inversions (7,12,15,16,20)

#```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

plot_af_chromosome <- function(pop,years,win,step,chrs){
    
    freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    freqs_win <- freqs_win[freqs_win$chr_num %in% chrs,]
    
    p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
        geom_point(size = 1, color = gry,alpha = 0.7, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_smooth(method = "loess", se = FALSE, span = 1/10, color = red)+
        facet_wrap(~chr_num, ncol = 1)+
        ggtitle(paste(pop,years, sep = " "))
    print(p)
}


plot_af_chromosome_within <- function(pop,win,step,chrs){
    freqs_win <- read.csv(paste0("Output/freqs/",pop,"_shifts_",win,"_win_",step,"_step.csv"), row.names = 1)
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    freqs_win <- freqs_win[freqs_win$chr_num %in% chrs,]
    
    p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_smooth(aes_string(x = "win_mid", y = "zt01_mean", color = "grb"),method = "loess", se = FALSE, span = 1/10)+
        geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
        geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
        facet_wrap(~chr_num, ncol = 3)+
        scale_color_manual(values=c(lir,grb,org),name = "contrast",
                           labels = c("1991-1996", "1996-2007", "2007-2017"))+
        ylim(-0.05, 0.05)+
        ggtitle(pop)
    print(p)
}

# chr 7
pop <- "PW"
chrs <- c(7)

plot_af_chromosome(pop,"zt01_mean",win,step,chrs)
#plot_af_chromosome(pop,"zt12_mean",win,step,chrs)
#plot_af_chromosome(pop,"zt23_mean",win,step,chrs)
plot_af_chromosome_within(pop,win,step,chrs)


ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
    geom_point(size = 1, color = gry,alpha = 0.7, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_line(color = red)+
    facet_wrap(~chr_num, ncol = 1)
    ggtitle(paste(pop,years, sep = " "))

ggplot(freqs_win, aes_string(x = "win_mid", y = "zt01_mean")) + 
        geom_point(size = 1, color = gry,alpha = 0.7, shape = 1)+
        theme_minimal()+
        theme(axis.text.x=element_blank())+
        ylab(expression(paste(Delta, "z")))+ xlab("")+ 
        geom_line(color = red)+
        facet_wrap(~chr_num, ncol = 1)
    ggtitle(paste(pop,years, sep = " "))


pop <- "TB"
plot_af_chromosome_within(pop,win,step,chrs)


pop <- "SS"
freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
freqs_win <- freqs_win[freqs_win$chr_num %in% chrs,]

p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt12_mean")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
    facet_wrap(~chr_num, ncol = 3)+
    scale_color_manual(values=c(grb,org),name = "contrast",
                       labels = c("1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+
    ggtitle(pop)
print(p)


# all inversions
pop <- "PWS"
years <- "zt01_mean"
win <- "1mb"
step <- "100kb"
chrs <- c(7,12,15,16,20)

plot_af_chromosome(pop,"zt01_mean",win,step,chrs)
plot_af_chromosome_within(pop,win,step,chrs)

pop <- "TB"
plot_af_chromosome_within(pop,win,step,chrs)


pop <- "SS"
freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
freqs_win <- freqs_win[freqs_win$chr_num %in% chrs,]

p <- ggplot(freqs_win, aes_string(x = "win_mid", y = "zt12_mean")) + 
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+ 
    geom_smooth(aes_string(x = "win_mid", y = "zt12_mean", color = "lir"),method = "loess", se = FALSE, span = 1/10)+
    geom_smooth(aes_string(x = "win_mid", y = "zt23_mean", color = "org"),method = "loess", se = FALSE, span = 1/10)+
    facet_wrap(~chr_num, ncol = 3)+
    scale_color_manual(values=c(grb,org),name = "contrast",
                       labels = c("1996-2007", "2007-2017"))+
    ylim(-0.05, 0.05)+
    ggtitle(pop)
print(p)


##########################
## Identify the regions with big/opposite allele freq changes over time
pop <- "PWS"

freqs_win <- read.table(paste("Data/freqs/window/PWS_shifts_50kb_win_10kb_step.txt"), sep = "")
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")
# identify windows where allele frequencies increased between each sampling period
pws_increased_freq <- freqs_win[freqs_win$zt01_mean >= quantile(freqs_win$zt01_mean, .9)[[1]] &
                                    freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                                    freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]

# identify windows where allele frequencies decreased between each sampling period
pws_decreased_freq <- freqs_win[freqs_win$zt01_mean <= quantile(freqs_win$zt01_mean, .1)[[1]] &
                                    freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                                    freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]

#print(paste(pop,"increased"))
#print(nrow(top_increased))# / nrow(freqs_win))
#print(paste(pop,"decreased"))
#print(nrow(top_decreased))# / nrow(freqs_win))

#hist(c(top_increased$zt01_mean,
#       top_increased$zt12_mean,
#       top_increased$zt23_mean,
#       top_decreased$zt01_mean,
#       top_decreased$zt12_mean,
#       top_decreased$zt23_mean), xlab = expression(paste(Delta, "z")),
#       main = pop, xlim = c(-0.25, 0.3))


pop <- "TB"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")
# identify windows where allele frequencies increased between each sampling period
tb_increased_freq <- freqs_win[freqs_win$zt01_mean >= quantile(freqs_win$zt01_mean, .9)[[1]] &
                                   freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                                   freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]

# identify windows where allele frequencies decreased between each sampling period
tb_decreased_freq <- freqs_win[freqs_win$zt01_mean <= quantile(freqs_win$zt01_mean, .1)[[1]] &
                                   freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                                   freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]

#print(paste(pop,"increased"))
#print(nrow(top_increased))# / nrow(freqs_win))
#print(paste(pop,"decreased"))
#print(nrow(top_decreased))# / nrow(freqs_win))

#hist(c(top_increased$zt01_mean,
#       top_increased$zt12_mean,
#       top_increased$zt23_mean,
#       top_decreased$zt01_mean,
#       top_decreased$zt12_mean,
#       top_decreased$zt23_mean), xlab = expression(paste(Delta, "z")),
#       main = pop, xlim = c(-0.25, 0.3))


pop <- "SS"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")
# identify windows where allele frequencies increased between each sampling period
ss_increased_freq <- freqs_win[freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                                   freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]

# identify windows where allele frequencies decreased between each sampling period
ss_decreased_freq <- freqs_win[freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                                   freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]



pws_increased_freq$windex <- paste(pws_increased_freq$seqnames, pws_increased_freq$win_mid, sep = ":")
pws_decreased_freq$windex <- paste(pws_decreased_freq$seqnames, pws_decreased_freq$win_mid, sep = ":")

tb_increased_freq$windex <- paste(tb_increased_freq$seqnames, tb_increased_freq$win_mid, sep = ":")
tb_decreased_freq$windex <- paste(tb_decreased_freq$seqnames, tb_decreased_freq$win_mid, sep = ":")

ss_increased_freq$windex <- paste(ss_increased_freq$seqnames, ss_increased_freq$win_mid, sep = ":")
ss_decreased_freq$windex <- paste(ss_decreased_freq$seqnames, ss_decreased_freq$win_mid, sep = ":")


venn(list(PWS_up=pws_increased_freq$windex,
          TB_up=tb_increased_freq$windex,
          SS_up=ss_increased_freq$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)

venn(list(PWS_down=pws_decreased_freq$windex,
          TB_down=tb_decreased_freq$windex,
          SS_down=ss_decreased_freq$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)



```

## PWS chr 18 (windows with strongest positive shift across all sampling periods)

```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

pop <- "PWS"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

# identify windows where allele frequencies increased between each sampling period
top_increased <- freqs_win[freqs_win$zt01_mean >= quantile(freqs_win$zt01_mean, .9)[[1]] &
                               freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                               freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]

# identify windows where allele frequencies decreased between each sampling period
top_decreased <- freqs_win[freqs_win$zt01_mean <= quantile(freqs_win$zt01_mean, .1)[[1]] &
                               freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                               freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]

# PWS chr with strongest positive shift (18)

freqs_win1 <- freqs_win[freqs_win$chr == "18",]
top_increased1 <- top_increased[top_increased$chr == "18",]

p01 <- ggplot(freqs_win1, aes(x = win_mid, y = zt01_mean)) + 
    geom_point(size = 2, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+
    geom_vline(aes(xintercept = win_mid), top_increased1)
p12 <- ggplot(freqs_win1, aes(x = win_mid, y = zt12_mean)) + 
    geom_point(size = 2, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+
    geom_vline(aes(xintercept = win_mid), top_increased1)
p23 <- ggplot(freqs_win1, aes(x = win_mid, y = zt23_mean)) + 
    geom_point(size = 2, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab(expression(paste(Delta, "z")))+ xlab("")+
    geom_vline(aes(xintercept = win_mid), top_increased1)

#p + geom_vline(aes(xintercept = win_mid), tops)

ggarrange(p01,p12,p23, ncol = 1, nrow = 3)#,
#labels = c("A", "B", "C"))

# Find genes in regions

gff <- read.delim("Data/annotations/Clupea_harengus.Ch_v2.0.2.100.gff3",header = FALSE, sep = "\t",comment.char = '#')
#head(gff)
colnames(gff) <- c("seqnames","source","type","start","end","score","strand","phase","attributes")
genes <- gff[gff$type=="gene",]
genes <- genes %>% separate(attributes, c("at1","at2"), remove = FALSE,sep = "Name=") %>% separate(at2, c("gene_name","at3"), extra = "drop",sep = ";")
genes <- select(genes, -at1,-at3)
g <- genes[genes$strand == "\\.",]
#head(genes)
genes <- droplevels(genes)

genes <- genes %>% as_granges()

names(top_increased)[names(top_increased) == "chr"] <- "seqnames"
names(top_increased)[names(top_increased) == "win_start"] <- "start"
names(top_increased)[names(top_increased) == "win_end"] <- "end"

names(top_decreased)[names(top_decreased) == "chr"] <- "seqnames"
names(top_decreased)[names(top_decreased) == "win_start"] <- "start"
names(top_decreased)[names(top_decreased) == "win_end"] <- "end"


tops <- rbind(top_increased,top_decreased)

tops<- tops %>% as_granges()

genes_overlap <- join_overlap_intersect(tops, genes) %>% as.data.frame()

#write.table(as.data.frame(unique(na.omit(focal_chr_genes$gene_name))), "Data/GO/fst_outliers/PWS17_v_TB17_Fst_99percentile_genes.txt", quote= FALSE, row.names = FALSE, col.names = FALSE)


```


# Identify windows where allele frequencies showed shifts after collapse and reversed over time:

- increased between 91 and 96, decreased between 96 and 07, and continued to decrease between 07 and 17
- decreased between 91 and 96, increased between 96 and 07, and continued to increase between 07 and 17


```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

pop <- "PWS"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
pws_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= 0 &
                              freqs_win$zt12_mean <= 0 &
                              freqs_win$zt23_mean <= 0,]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
pws_dn_up_up <- freqs_win[freqs_win$zt01_mean <= 0 &
                              freqs_win$zt12_mean >= 0 &
                              freqs_win$zt23_mean >= 0,]

#print(paste(pop,"up_dn_dn"))
#up_dn_dn# / nrow(freqs_win)
#print(paste(pop,"dn_up_up"))
#dn_up_up# / nrow(freqs_win)
#
#print(hist(c(up_dn_dn$zt01_mean,
#             up_dn_dn$zt12_mean,
#             up_dn_dn$zt23_mean,
#             dn_up_up$zt01_mean,
#             dn_up_up$zt12_mean,
#             dn_up_up$zt23_mean), xlab = expression(paste(Delta, "z")),
#             main = pop))#, xlim = c(-0.25, 0.3)))

pop <- "TB"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
tb_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= 0 &
                             freqs_win$zt12_mean <= 0 &
                             freqs_win$zt23_mean <= 0,]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
tb_dn_up_up <- freqs_win[freqs_win$zt01_mean <= 0 &
                             freqs_win$zt12_mean >= 0 &
                             freqs_win$zt23_mean >= 0,]


pws_up_dn_dn$windex <- paste(pws_up_dn_dn$seqnames, pws_up_dn_dn$win_mid, sep = ":")
pws_dn_up_up$windex <- paste(pws_dn_up_up$seqnames, pws_dn_up_up$win_mid, sep = ":")

tb_up_dn_dn$windex <- paste(tb_up_dn_dn$seqnames, tb_up_dn_dn$win_mid, sep = ":")
tb_dn_up_up$windex <- paste(tb_dn_up_up$seqnames, tb_dn_up_up$win_mid, sep = ":")


venn(list(PWS_up_dn_dn=pws_up_dn_dn$windex,
          TB_up_dn_dn=tb_up_dn_dn$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)

venn(list(PWS_dn_up_up=pws_dn_up_up$windex,
          TB_dn_up_up=tb_dn_up_up$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)




```


# Identify windows where allele frequencies showed LARGE shifts after collapse and reversed over time:

- increased between 91 and 96, decreased between 96 and 07, and continued to decrease between 07 and 17
- decreased between 91 and 96, increased between 96 and 07, and continued to increase between 07 and 17


```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

pop <- "PWS"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
pws_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= quantile(freqs_win$zt01_mean, .9)[[1]] &
                              freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                              freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
pws_dn_up_up <- freqs_win[freqs_win$zt01_mean <= quantile(freqs_win$zt01_mean, .1)[[1]] &
                              freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                              freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]

#print(paste(pop,"up_dn_dn"))
#up_dn_dn# / nrow(freqs_win)
#print(paste(pop,"dn_up_up"))
#dn_up_up# / nrow(freqs_win)
#
#print(hist(c(up_dn_dn$zt01_mean,
#             up_dn_dn$zt12_mean,
#             up_dn_dn$zt23_mean,
#             dn_up_up$zt01_mean,
#             dn_up_up$zt12_mean,
#             dn_up_up$zt23_mean), xlab = expression(paste(Delta, "z")),
#             main = pop))#, xlim = c(-0.25, 0.3)))

pop <- "TB"

freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_",win,"_win_",step,"_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
tb_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= quantile(freqs_win$zt01_mean, .9)[[1]] &
                             freqs_win$zt12_mean <= quantile(freqs_win$zt12_mean, .1)[[1]] &
                             freqs_win$zt23_mean <= quantile(freqs_win$zt23_mean, .1)[[1]],]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
tb_dn_up_up <- freqs_win[freqs_win$zt01_mean <= quantile(freqs_win$zt01_mean, .1)[[1]] &
                             freqs_win$zt12_mean >= quantile(freqs_win$zt12_mean, .9)[[1]] &
                             freqs_win$zt23_mean >= quantile(freqs_win$zt23_mean, .9)[[1]],]


pws_up_dn_dn$windex <- paste(pws_up_dn_dn$seqnames, pws_up_dn_dn$win_mid, sep = ":")
pws_dn_up_up$windex <- paste(pws_dn_up_up$seqnames, pws_dn_up_up$win_mid, sep = ":")

tb_up_dn_dn$windex <- paste(tb_up_dn_dn$seqnames, tb_up_dn_dn$win_mid, sep = ":")
tb_dn_up_up$windex <- paste(tb_dn_up_up$seqnames, tb_dn_up_up$win_mid, sep = ":")


venn(list(PWS_up_dn_dn=pws_up_dn_dn$windex,
          TB_up_dn_dn=tb_up_dn_dn$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)

venn(list(PWS_dn_up_up=pws_dn_up_up$windex,
          TB_dn_up_up=tb_dn_up_up$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)




```

# Permutation tests to determine significant outlier windows.

Randomly sample individuals from t0 and t1 and measure change in allele frequencies in 50kb windows. Identify genome-wide mean of means across windows and quantiles across 1000 permutations. See [AF_permutations.R](https://github.com/joemcgirr/pac_herring/blob/master/Rmarkdown/AF_shifts/AF_permutations.R)


```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7,class.source = 'fold-show'}

pop <- "PWS"


perm_t01 <- read.table(paste("Data/freqs/perms/",pop,"91_",pop,"96_perms.txt", sep = ""), header = TRUE)
perm_t12 <- read.table(paste("Data/freqs/perms/",pop,"96_",pop,"07_perms.txt", sep = ""), header = TRUE)
perm_t23 <- read.table(paste("Data/freqs/perms/",pop,"07_",pop,"17_perms.txt", sep = ""), header = TRUE)


freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_50kb_win_10kb_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

c(min(freqs_win$zt01_mean),max(freqs_win$zt01_mean))
c(mean(perm_t01$min),mean(perm_t01$max))

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
pws_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= quantile(perm_t01$q_90, .9)[[1]] &
                              freqs_win$zt12_mean <= quantile(perm_t12$q_10, .1)[[1]] &
                              freqs_win$zt23_mean <= quantile(perm_t23$q_10, .1)[[1]],]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
pws_dn_up_up <- freqs_win[freqs_win$zt01_mean <= quantile(perm_t01$q_10, .1)[[1]] &
                              freqs_win$zt12_mean >= quantile(perm_t12$q_90, .9)[[1]] &
                              freqs_win$zt23_mean >= quantile(perm_t23$q_90, .9)[[1]],]

#print(paste(pop,"up_dn_dn"))
#up_dn_dn# / nrow(freqs_win)
#print(paste(pop,"dn_up_up"))
#dn_up_up# / nrow(freqs_win)
#
#print(hist(c(up_dn_dn$zt01_mean,
#             up_dn_dn$zt12_mean,
#             up_dn_dn$zt23_mean,
#             dn_up_up$zt01_mean,
#             dn_up_up$zt12_mean,
#             dn_up_up$zt23_mean), xlab = expression(paste(Delta, "z")),
#             main = pop))#, xlim = c(-0.25, 0.3)))

pop <- "TB"

perm_t01 <- read.table(paste("Data/freqs/perms/",pop,"91_",pop,"96_perms.txt", sep = ""), header = TRUE)
perm_t12 <- read.table(paste("Data/freqs/perms/",pop,"96_",pop,"06_perms.txt", sep = ""), header = TRUE)
perm_t23 <- read.table(paste("Data/freqs/perms/",pop,"06_",pop,"17_perms.txt", sep = ""), header = TRUE)


freqs_win <- read.table(paste("Data/freqs/window/",pop,"_shifts_50kb_win_10kb_step.txt", sep = ""))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
colnames(freqs_win)[c(1:3)] <- c("seqnames","start","end")

# identify windows where allele frequencies increased between 91 and 96
# and decreased between 96 and 07
# and continued to decrease between 07 and 17
tb_up_dn_dn <- freqs_win[freqs_win$zt01_mean >= quantile(perm_t01$q_90, .9)[[1]] &
                             freqs_win$zt12_mean <= quantile(perm_t12$q_10, .1)[[1]] &
                             freqs_win$zt23_mean <= quantile(perm_t23$q_10, .1)[[1]],]

# identify windows where allele frequencies decreased between 91 and 96
# and increased between 96 and 07
# and continued to increase between 07 and 17
tb_dn_up_up <- freqs_win[freqs_win$zt01_mean <= quantile(perm_t01$q_10, .1)[[1]] &
                             freqs_win$zt12_mean >= quantile(perm_t12$q_90, .9)[[1]] &
                             freqs_win$zt23_mean >= quantile(perm_t23$q_90, .9)[[1]],]


pws_up_dn_dn$windex <- paste(pws_up_dn_dn$seqnames, pws_up_dn_dn$win_mid, sep = ":")
pws_dn_up_up$windex <- paste(pws_dn_up_up$seqnames, pws_dn_up_up$win_mid, sep = ":")

tb_up_dn_dn$windex <- paste(tb_up_dn_dn$seqnames, tb_up_dn_dn$win_mid, sep = ":")
tb_dn_up_up$windex <- paste(tb_dn_up_up$seqnames, tb_dn_up_up$win_mid, sep = ":")


p<-ggplot(all_pops, aes(x=zt01_t)) +
    geom_density(alpha=0.4)+
    scale_fill_manual(name = "Population",values=c(red,yel,blu))+
    ylab("density\n") + xlab(expression(paste("TB",Delta, "z 1991-1996")))+
    theme_minimal()+
    geom_vline(xintercept = quantile(perm_t01$q_10, .1), color = red)+
    geom_vline(xintercept = quantile(perm_t01$q_90, .9), color = red)+
    geom_vline(xintercept = quantile(freqs_win$zt01_mean, .1), color = blu)+
    geom_vline(xintercept = quantile(freqs_win$zt01_mean, .9), color = blu)
p

venn(list(PWS_up_dn_dn=pws_up_dn_dn$windex,
          TB_up_dn_dn=tb_up_dn_dn$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)

venn(list(PWS_dn_up_up=pws_dn_up_up$windex,
          TB_dn_up_up=tb_dn_up_up$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)


```

# Fst windows
Which loci are consistently differentiated between sampling periods?
    
    Find windows that show high Fst in 1991 v 1996, 1996 v 2007, and 2007 - 2017 comparisons
```{r, show = FALSE,message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7}
fsts <- list.files(path="Data/fst_pbs/maf05/", full.names=T)
#fsts <- fsts[1]

fst <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt", header = TRUE, stringsAsFactors = FALSE)
fst2 <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt", header = TRUE, stringsAsFactors = FALSE)
names(fst2) <- c("region","chr","midPos","Nsites","Fst02","Fst03","Fst23","PBS0","PBS2","PBS3")

fst$Fst03 <- fst2$Fst03
fst$Fst23 <- fst2$Fst23

fst <- fst %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
fst$chr_num <- factor(fst$chr, levels = c(1:26))

high_fst <- fst[fst$Fst01 >= quantile(fst$Fst01, .90)[[1]] &
                    fst$Fst12 >= quantile(fst$Fst12, .90)[[1]] &
                    fst$Fst23 >= quantile(fst$Fst23, .90)[[1]],]
high_fst$windex <- paste(high_fst$chr,high_fst$midPos,sep = ":")
pws <- high_fst
head(high_fst)
nrow(high_fst)



fst <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB06.txt", header = TRUE, stringsAsFactors = FALSE)
fst2 <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB17.txt", header = TRUE, stringsAsFactors = FALSE)
names(fst2) <- c("region","chr","midPos","Nsites","Fst02","Fst03","Fst23","PBS0","PBS2","PBS3")

fst$Fst03 <- fst2$Fst03
fst$Fst23 <- fst2$Fst23

fst <- fst %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
fst$chr_num <- factor(fst$chr, levels = c(1:26))

high_fst <- fst[fst$Fst01 >= quantile(fst$Fst01, .90)[[1]] &
                    fst$Fst12 >= quantile(fst$Fst12, .90)[[1]] &
                    fst$Fst23 >= quantile(fst$Fst23, .90)[[1]],]
high_fst$windex <- paste(high_fst$chr,high_fst$midPos,sep = ":")
tb <- high_fst

head(high_fst)
nrow(high_fst)

venn(list(PWS=pws$windex,
          TB=tb$windex), 
     ilabels = TRUE,opacity = 0.5,box = FALSE,
     zcolor = c(red,blu,yel), ilcs = 1)


```


# Fst matrix temporal samples (within population comparisons)

## PWS
```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7}

fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt")
p91xp96 <- round(mean(fst$Fst01),4)
p91xp07 <- round(mean(fst$Fst02),4)
p96xp07 <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS17.txt")
p91xp17 <- round(mean(fst$Fst02),4)
p96xp17 <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt")
p07xp17 <- round(mean(fst$Fst12),4)

fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
             p91xp96,0,p96xp07,p96xp17,
             p91xp07,p96xp07,0,p07xp17,
             p91xp17,p96xp17,p07xp17,0)

fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
colnames(fst_mat) <- c("PWS91","PWS96","PWS07","PWS17")
rownames(fst_mat) <- c("PWS91","PWS96","PWS07","PWS17")
print(fst_mat)
upper_tri <- get_upper_tri(fst_mat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = blu, high = red, 
                        limit = c(min(c(fst_mat)[c(fst_mat) > 0]),max(fst_mat)), space = "Lab", 
                        name="Fst") +
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 12, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 5)

```

## TB

```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7}

fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB06.txt")
p91xp96 <- round(mean(fst$Fst01),4)
p91xp07 <- round(mean(fst$Fst02),4)
p96xp07 <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB17.txt")
p91xp17 <- round(mean(fst$Fst02),4)
p96xp17 <- round(mean(fst$Fst12),4)
fst <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB06_TB17.txt")
p07xp17 <- round(mean(fst$Fst12),4)

fst_vec <- c(0,p91xp96,p91xp07,p91xp17,
             p91xp96,0,p96xp07,p96xp17,
             p91xp07,p96xp07,0,p07xp17,
             p91xp17,p96xp17,p07xp17,0)

fst_mat = matrix(fst_vec, nrow = 4, ncol = 4)
colnames(fst_mat) <- c("TB91","TB96","TB06","TB17")
rownames(fst_mat) <- c("TB91","TB96","TB06","TB17")
print(fst_mat)
upper_tri <- get_upper_tri(fst_mat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = blu, high = red, 
                        limit = c(min(c(fst_mat)[c(fst_mat) > 0]),max(fst_mat)), space = "Lab", 
                        name="Fst") +
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 12, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 5)

```

## SS

```{r, message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7}

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
print(fst_mat)
upper_tri <- get_upper_tri(fst_mat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = blu, high = red, 
                        limit = c(min(c(fst_mat)[c(fst_mat) > 0]),max(fst_mat)), space = "Lab", 
                        name="Fst") +
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 12, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 5)

```



```{r, show = FALSE,message=FALSE, warning=FALSE, fig.width= 10, fig.height= 7}
fsts <- list.files(path="Data/fst_pbs/maf05/", full.names=T)
#fsts <- fsts[1]


fst <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt", header = TRUE, stringsAsFactors = FALSE)
fst2 <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt", header = TRUE, stringsAsFactors = FALSE)
names(fst2) <- c("region","chr","midPos","Nsites","Fst02","Fst03","Fst23","PBS0","PBS2","PBS3")
#fst1 <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB06.txt", header = TRUE, stringsAsFactors = FALSE)
#fst2 <- read.table("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_TB91_TB96_TB17.txt", header = TRUE, stringsAsFactors = FALSE)

fst$Fst03 <- fst2$Fst03
fst$Fst23 <- fst2$Fst23

fst <- fst %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
fst$chr_num <- factor(fst$chr, levels = c(1:26))


# identify windows where allele frequencies increased between each sampling period
high_fst <- fst[fst$Fst01 >= quantile(fst$Fst01, .99)[[1]] &
                    fst$Fst12 >= quantile(fst$Fst12, .99)[[1]] &
                    fst$Fst23 >= quantile(fst$Fst23, .99)[[1]],]

head(high_fst)
nrow(high_fst)

p <- ggplot(fst, aes(x = midPos, y = Fst01)) + 
    geom_point(size = 1, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab("Fst\n")+ xlab("")+ 
    geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
    facet_wrap(~chr_num, ncol = 9)+
    ggtitle("1991 - 1996")

p + geom_vline(aes(xintercept = midPos), high_fst)

p <- ggplot(fst, aes(x = midPos, y = Fst12)) + 
    geom_point(size = 1, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab("Fst\n")+ xlab("")+ 
    geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
    facet_wrap(~chr_num, ncol = 9)+
    ggtitle("1996 - 2007")

p + geom_vline(aes(xintercept = midPos), high_fst)

p <- ggplot(fst, aes(x = midPos, y = Fst23)) + 
    geom_point(size = 1, color = gry,alpha = 0.6, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab("Fst\n")+ xlab("")+ 
    geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
    facet_wrap(~chr_num, ncol = 9)+
    ggtitle("2007 - 2017")

p + geom_vline(aes(xintercept = midPos), high_fst)


```


