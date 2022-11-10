#Create sliding window allele freq files from angsd
source("Rscripts/BaseScripts.R")
library(reshape2)
library(windowscanr)


# Prince William Sound
t0 <- read.table("Data/new_vcf/AF/PWS91.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t1 <- read.table("Data/new_vcf/AF/PWS96.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t2 <- read.table("Data/new_vcf/AF/PWS07.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t3 <- read.table("Data/new_vcf/AF/PWS17.mafs",stringsAsFactors = FALSE, sep="\t", header = T)


freqs <- data.frame(chr = t0$chromo,
                    pos = t0$position,
                    t0_AF = t0$knownEM,
                    t1_AF = t1$knownEM,
                    t2_AF = t2$knownEM,
                    t3_AF = t3$knownEM)

freqs$t0_transformed_freq <- asin(sqrt(freqs$t0_AF))
freqs$t1_transformed_freq <- asin(sqrt(freqs$t1_AF))
freqs$t2_transformed_freq <- asin(sqrt(freqs$t2_AF))
freqs$t3_transformed_freq <- asin(sqrt(freqs$t3_AF))

#Allele frequency change over time
freqs$zt01 <- freqs$t1_transformed_freq - freqs$t0_transformed_freq
freqs$zt02 <- freqs$t2_transformed_freq - freqs$t0_transformed_freq
freqs$zt03 <- freqs$t3_transformed_freq - freqs$t0_transformed_freq
freqs$zt12 <- freqs$t2_transformed_freq - freqs$t1_transformed_freq
freqs$zt13 <- freqs$t3_transformed_freq - freqs$t1_transformed_freq
freqs$zt23 <- freqs$t3_transformed_freq - freqs$t2_transformed_freq
write.csv(freqs, "Data/freqs/PWS_md7000_shifts_persite.csv")

#TB
t0 <- read.table("Data/new_vcf/AF/TB91.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t1 <- read.table("Data/new_vcf/AF/TB96.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t2 <- read.table("Data/new_vcf/AF/TB06.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t3 <- read.table("Data/new_vcf/AF/TB17.mafs",stringsAsFactors = FALSE, sep="\t", header = T)

freqs <- data.frame(chr = t0$chromo,
                    pos = t0$position,
                    t0_AF = t0$knownEM,
                    t1_AF = t1$knownEM,
                    t2_AF = t2$knownEM,
                    t3_AF = t3$knownEM)

freqs$t0_transformed_freq <- asin(sqrt(freqs$t0_AF))
freqs$t1_transformed_freq <- asin(sqrt(freqs$t1_AF))
freqs$t2_transformed_freq <- asin(sqrt(freqs$t2_AF))
freqs$t3_transformed_freq <- asin(sqrt(freqs$t3_AF))

#Allele frequency change over time
freqs$zt01 <- freqs$t1_transformed_freq - freqs$t0_transformed_freq
freqs$zt02 <- freqs$t2_transformed_freq - freqs$t0_transformed_freq
freqs$zt03 <- freqs$t3_transformed_freq - freqs$t0_transformed_freq
freqs$zt12 <- freqs$t2_transformed_freq - freqs$t1_transformed_freq
freqs$zt13 <- freqs$t3_transformed_freq - freqs$t1_transformed_freq
freqs$zt23 <- freqs$t3_transformed_freq - freqs$t2_transformed_freq
write.csv(freqs, "Data/freqs/TB_md7000_shifts_persite.csv")


# Stika Sound
t1 <- read.table("Data/new_vcf/AF/SS96.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t2 <- read.table("Data/new_vcf/AF/SS06.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
t3 <- read.table("Data/new_vcf/AF/SS17.mafs",stringsAsFactors = FALSE, sep="\t", header = T)

freqs <- data.frame(chr = t0$chromo,
                    pos = t0$position,
                    t1_AF = t1$knownEM,
                    t2_AF = t2$knownEM,
                    t3_AF = t3$knownEM)

freqs$t1_transformed_freq <- asin(sqrt(freqs$t1_AF))
freqs$t2_transformed_freq <- asin(sqrt(freqs$t2_AF))
freqs$t3_transformed_freq <- asin(sqrt(freqs$t3_AF))

#Allele frequency change over time
freqs$zt12 <- freqs$t2_transformed_freq - freqs$t1_transformed_freq
freqs$zt13 <- freqs$t3_transformed_freq - freqs$t1_transformed_freq
freqs$zt23 <- freqs$t3_transformed_freq - freqs$t2_transformed_freq
write.csv(freqs, "Data/freqs/SS_md7000_shifts_persite.csv")



#window based AF
for (i in c("PWS","TB")){
    window_size <- 500000
    step_size <- 50000
    
    df<-read.csv(paste0("Data/freqs/",i, "_md7000_shifts_persite.csv",row.names=1))
    freqs_win <- winScan(x = df, 
                         groups = "chr", 
                         position = "pos", 
                         values = c("zt01","zt02","zt03","zt12","zt13","zt23"), 
                         win_size = window_size,
                         win_step = step_size,
                         funs = c("mean"))
    freqs_win <- na.omit(freqs_win)
    freqs_win <- freqs_win %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
    freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))
    
    write.csv(freqs_win, paste0("Data/freqs/", i, "_shifts_500k_50k.csv"))
    
}

#SS
window_size <- 500000
step_size <- 50000
df<-read.csv("Data/freqs//SS_md7000_shifts_persite.csv", row.names = 1)
freqs_win <- winScan(x = df, 
                     groups = "chr", 
                     position = "pos", 
                     values = c("zt12","zt13","zt23"), 
                     win_size = window_size,
                     win_step = step_size,
                     funs = c("mean"))
freqs_win <- na.omit(freqs_win)
freqs_win <- freqs_win %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

write.csv(freqs_win, "Data/freqs/MD7000/SS_shifts_500k_50k.txt")



######  Compare VCF maf vs. ANGSD estimated maf
#compare allele frequencies estimated from angsd and extracted from vcf files
col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF")
v1<-read.table("Data/new_vcf/AF/PWS91_maf05_af.frq",stringsAsFactors = FALSE, sep="\t", header = F, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF") )
v1$maf<-as.numeric(substr(v1$MAF, 3,10))
v1$Maj<-as.numeric(substr(v1$MajorAF, 3,10))

#Read maf01 .maf file
#t0$Maj<-1-t0$knownEM

af<-t0_2
af$vcf.af<-v1$maf

t91<-af[,c("chromo","position", "knownEM","vcf.af")]
t91m<-melt(t91, id.vars=c("chromo","position"))

t91m$chromo<-factor(t91m$chromo, levels=paste0("chr",1:26))

ggplot(t91m, aes(x=position, y=value, color=variable))+
    facet_wrap(~chromo, ncol=6)+
    geom_point(size=0.2, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')+theme(legend.title = element_blank())

ggsave("Output/VCF/AF_comparison_pws91_angsd.vs.vcf.png", width=15, height=10, dpi=150)


#chr1
ggplot(t91m[t91m$chromo=="chr1",], aes(x=position, y=value, color=variable))+
    geom_point(size=0.2, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')+theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 2, alpha=0.8)))+
    scale_color_discrete(labels=c('ANGSD', 'VCF'))

ggsave("Output/VCF/AF_comparison_pws91_angsd.vs.vcf_chr1.png", width=7, height=3.5, dpi=150)
#VCF -cutoff is maf0.05 for all 14 samples together -> within populations may have freq=0
#Angsd-estiamted -> maf does not go down to 0?

min(t91$knownEM) #0
min(t91$vcf.af) #0

### Plot allele frequency <0.05
t91_lowf<-t91[t91$knownEM<=0.05,]
t91_lowf$index<-1:nrow(t91_lowf)
ggplot()+
    geom_point(data=t91_lowf, aes(x=index, y=knownEM), color="#F8766D",size=0.2, alpha=0.5)+
    #facet_wrap(~chromo, ncol=6)+
    geom_point(data=t91_lowf, aes(x=index, y=vcf.af), color="gray",size=0.2, alpha=0.5)+
    theme_minimal()+
    xlab('')+theme(legend.title = element_blank())
ggsave("Output/VCF/AFbelow0.05_angsd.knownEM.png", width=15, height=8, dpi=100)

ggplot()+
    geom_point(data=t91_lowf[t91_lowf$chromo=="chr1",], aes(x=position, y=vcf.af, color="VCF"),color="#00BFC4",size=0.3, alpha=0.5)+
    geom_point(data=t91_lowf[t91_lowf$chromo=="chr1",], aes(x=position, y=-knownEM, color="ANGSD"),color="#F8766D",size=0.3, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')+
    annotate('text', x=0, y=-0.005, label="ANGSD", hjust=0)+
    annotate('text', x=0, y=0.005,  label="vcftools", hjust=0)
ggsave("Output/VCF/AFbelow0.05_angsd.knownEM_PWS91chr1.png", width=7, height=4, dpi=150)

#facet lables
methods <- c(`knownEM` = "ANGSD", `vcf.af` = "vcftools")

ggplot(t91m[t91m$chromo=="chr1",], aes(x=position, y=value, color=variable))+
    facet_wrap(~variable, ncol=1, labeller = as_labeller(methods))+
    geom_point(size=0.2, alpha=0.5)+
    theme_light()+
    ylab("MAF")+xlab('')+theme(legend.title = element_blank())+
    scale_color_discrete(guide="none")
ggsave("Output/VCF/AF.comapare.angsd.vs.vcftools.pws91_facet.pdf", width = 6, height = 6)

length(t91$knownEM[t91$knownEM==0]) #14
length(t91$vcf.af[t91$vcf.af==0]) #11730

#Plot only af<0.05
ggplot(t91m[t91m$chromo=="chr1"&t91m$value<0.05,], aes(x=position, y=value, color=variable))+
    facet_wrap(~variable, ncol=1, labeller = as_labeller(methods))+
    geom_point(size=0.2, alpha=0.5)+
    theme_light()+
    ylab("MAF")+xlab('')+theme(legend.title = element_blank())+
    scale_color_discrete(guide="none")
ggsave("Output/VCF/AF.comapare.angsd.vs.vcftools.pws91.less0.05_facet.png", width = 5, height = 4, dpi=150)






#Use maf01 .maf file form ANGSD
t0_2<-read.table("Data/new_vcf/AF/maf01/PWS91_maf01.mafs",stringsAsFactors = FALSE, sep="\t", header = T)

ggplot()+
    geom_point(data=t91[t91$chromo=="chr1",], aes(x=position, y=vcf.af, color="#F8766D"),size=0.2, alpha=0.5)+
    geom_point(data=t0_2[t0_2$chromo=="chr1",], aes(x=position, y=knownEM, color="#7CAE00"),size=0.2, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')

ggplot()+
    geom_point(data=t91[t91$chromo=="chr1",], aes(x=position, y=vcf.af, color="VCF"),color="#7CAE00",size=0.2, alpha=0.5)+
    geom_point(data=t0_2[t0_2$chromo=="chr1",], aes(x=position, y=-knownEM, color="ANGSD"),color="#F8766D",size=0.2, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')












####
df<-read.csv("~/Projects/Pacherring_Vincent/Data/Af_angsd.csv", header=F)
df2<-1-df
df2<-t(df2)
write.table(df2, "~/Projects/Pacherring_Vincent/Data/Af_angsd_major_transposed.csv", row.names = F, col.names = F, , sep = ",")


t91cut<-t91[t91$knownEM>=0.05,]
#260225 loci (as opposed to 330482)
t91cutm<-melt(t91cut, id.vars=c("chromo","position"))
ggplot(t91cutm[t91cutm$chromo=="chr1",], aes(x=position, y=value, color=variable))+
    facet_wrap(~chromo, ncol=6)+
    geom_point(size=0.2, alpha=0.5)+
    theme_minimal()+
    ylab("MAF")+xlab('')+theme(legend.title = element_blank())


freq<-data.frame(t0=t0$knownEM, t1=t1$knownEMnEM,t1=t1$knownEMnEM,t3=t3$knownEMnEM)


## convert the estiamted allele freq to maf matrix and remove the sites<0.05

#pws
for (i in 2:4){
    df<-get(paste0("t",i-1))
    df<-df[df$knownEM>=0.05,]
    df[,paste0("maf",i-1)]<-1-df$knownEM
    df$pos<-paste0(df$chromo,".",df$position)
    
    if(i==1) PWS<-df
    if (i!=1) {
        PWS<-merge(PWS, df[,c("pos",paste0("maf",i-1))], by="pos")
    }
}

#ss
s1<-read.table("Data/new_vcf/AF/SS96.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
s2<-read.table("Data/new_vcf/AF/SS06.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
s3<-read.table("Data/new_vcf/AF/SS17.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
for (i in 1:3){
    df<-get(paste0("s",i))
    df[,paste0("maf",i)]<-1-df$knownEM
    df<-df[df[,paste0("maf",i)]>=0.05,] #5571
    df$pos<-paste0(df$chromo,".",df$position)
    
    if(i==1) SS<-df
    if (i!=1) {
        SS<-merge(SS, df[,c("pos",paste0("maf",i))], by="pos")
    }
}

b0<-read.table("Data/new_vcf/AF/TB91.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
b1<-read.table("Data/new_vcf/AF/TB96.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
b2<-read.table("Data/new_vcf/AF/TB06.mafs",stringsAsFactors = FALSE, sep="\t", header = T)
b3<-read.table("Data/new_vcf/AF/TB17.mafs",stringsAsFactors = FALSE, sep="\t", header = T)

for (i in 1:4){
    df<-get(paste0("b",i-1))
    df[,paste0("maf",i-1)]<-1-df$knownEM
    df<-df[df[,paste0("maf",i-1)]>=0.05,] #5571
    df$pos<-paste0(df$chromo,".",df$position)
    
    if(i==1) TB<-df
    if (i!=1) {
        TB<-merge(TB, df[,c("pos",paste0("maf",i-1))], by="pos")
    }
}

nrow(PWS) #325139
nrow(SS) #325172
nrow(TB) #324284

#one matrix
AF<-merge(PWS[, c("pos","chromo","position","maf0","maf1","maf2","maf3")], SS[c("pos","maf1","maf2","maf3")], by="pos")
AF<-merge(AF, TB[, c("pos","maf0","maf1","maf2","maf3")], by="pos")
AF<-AF[order(AF$chromo, AF$position),]

write.table(AF[,2:12],"Data/AF_ansgd_Majaf.csv", sep=",",col.names =F, quote = F, row.names = F)

df<-t1
df$pos<-paste0(df$chromo,".",df$position)
removedsites<-df$pos[!(df$pos %in% AF$pos)]
