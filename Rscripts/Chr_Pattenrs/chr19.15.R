#### Plot which individuals have freq changes in this region
library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)
library(cowplot)

bla <- "#000000"
blu <- "#0072b2"
grb <- "#56b4e9"
lir <- "#cc79a7"
gre <- "#009e73"
red <- "#d55e00"
org <- "#e69f00"
yel <- "#f0e442"
gry<-  '#BBBBBB'

ch="chr19"
ch="chr15"

dir.create(paste0("Output/chr/DP7000/",ch,"/"))
vcf <- read.vcfR(paste0("Data/new_vcf/PH_MD7000_maf05_",ch,".vcf.gz"))

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[c("Sample","Population.Year")]
pops<-unique(pop_info$Population.Year)

#PWS populations
pws<-c("PWS91","PWS96","PWS07","PWS17")
tb<-c("TB91","TB96","TB06","TB17")
ss<-c("SS96","SS06","SS17")


Plots<-list()
#1
pop<-pws[1]
p_map<-pop_info[pop_info$Population.Year==pop,]

chr_plot <- genotype_plot(vcf_object  =  vcf, popmap = p_map,
                          snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),
                          plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot$genotypes<-updated    
Plots[[1]]<-combine_genotype_plot(chr_plot)

#2 & 3
for (i in 2:length(pws)){
    pop<-pws[i]
    p_map<-pop_info[pop_info$Population.Year==pop,]
    chr_plot <- genotype_plot(vcf_object  =  vcf, popmap = p_map,
                              snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),
                              plot_allele_frequency=F,invariant_filter = F)  
    if (i ==4){
    gt_plot<-chr_plot$genotypes   
    updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                       breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
    chr_plot$genotypes<-updated    
    Plots[[i]]<-chr_plot$genotypes
    }
    else {
        gt_plot<-chr_plot$genotypes   
        updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                           breaks=c("0","1","2"),guide="none")
        chr_plot$genotypes<-updated    
        Plots[[i]]<-chr_plot$genotypes
    }
}

pdf(paste0("Output/chr/DP7000/",ch,"/PWS_all_",ch,".pdf"),width=8,height=20)
ggdraw()+
    draw_plot(Plots[[1]],x=0,y=0.73,width=1 ,height=0.27)+
    draw_plot(Plots[[2]],x=0,y=0.5,width=1 ,height=0.23)+
    draw_plot(Plots[[3]],x=0,y=0.27,width=1 ,height=0.23)+
    draw_plot(Plots[[4]],x=0,y=0,width=1 ,height=0.27)
dev.off()


#SS populations
ss<-c("SS96","SS06","SS17")
#1.
p_map<-pop_info[pop_info$Population.Year==ss[1],]
chr_plot1 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  

gt_plot<-chr_plot1$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot1$genotypes<-updated
#2
p_map<-pop_info[pop_info$Population.Year==ss[2],]
chr_plot2 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot2$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot2$genotypes<-updated    
#3
p_map<-pop_info[pop_info$Population.Year==ss[3],]
chr_plot3 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot3$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
chr_plot3$genotypes<-updated    

pdf(paste0("Output/chr/DP7000/",ch,"/",ch,"_SS-all.pdf"), width =8, height = 15 )
ggdraw()+
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.65,width=1,height=0.35)+
    draw_plot(chr_plot2$genotypes,0,0.35,1,0.3)+
    draw_plot(chr_plot3$genotypes,0,0,1,0.35)
dev.off()

######TB
tb<-c("TB91","TB96","TB06","TB17")
#1.
p_map<-pop_info[pop_info$Population.Year==tb[1],]
chr_plot1 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  

gt_plot<-chr_plot1$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot1$genotypes<-updated
#2
p_map<-pop_info[pop_info$Population.Year==tb[2],]
chr_plot2 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot2$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot2$genotypes<-updated    
#3
p_map<-pop_info[pop_info$Population.Year==tb[3],]
chr_plot3 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot3$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot3$genotypes<-updated    

#4
p_map<-pop_info[pop_info$Population.Year==tb[4],]
chr_plot4 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot4$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
chr_plot4$genotypes<-updated    


pdf(paste0("Output/chr/DP7000/",ch,"/",ch,"_TB-all.pdf"), width =8, height = 20 )
ggdraw()+
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.73,width=1,height=0.27)+
    draw_plot(chr_plot2$genotypes,0,0.5,1,0.23)+
    draw_plot(chr_plot3$genotypes,0,0.27,1,0.23)+
    draw_plot(chr_plot4$genotypes,0,0,1,0.27)
dev.off()

### CA,BC,WA

ot<-c("BC17","WA17","CA17")
#1.
p_map<-pop_info[pop_info$Population.Year==ot[1],]
chr_plot1 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  

gt_plot<-chr_plot1$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot1$genotypes<-updated


p_map<-pop_info[pop_info$Population.Year==ot[2],]
chr_plot2 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot2$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot2$genotypes<-updated    
#3
p_map<-pop_info[pop_info$Population.Year==ot[3],]
chr_plot3 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot3$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
chr_plot3$genotypes<-updated    

pdf(paste0("Output/chr/DP7000/",ch,"/",ch,"_CA.WA.BC.pdf"), width =8, height = 15 )
ggdraw()+
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.65,width=1,height=0.35)+
    draw_plot(chr_plot2$genotypes,0,0.35,1,0.3)+
    draw_plot(chr_plot3$genotypes,0,0,1,0.35)
dev.off()


###### FInd breakpoints

ca<-pop_info[pop_info$Population.Year=="CA17",]

chr16_plot <- genotype_plot(vcf_object  =  vcf,   
                           popmap = ca,                           
                           snp_label_size = 1000000,                     
                           colour_scheme=c(yel,"red","blue"),
                           plot_allele_frequency=TRUE,                    
                           invariant_filter = F)   
pdf("Output/chr/DP7000/chr16/AF_chr16_CA17.pdf", width = 10, height = 5)
combine_genotype_plot(chr16_plot)
dev.off()


chr16_plot2 <- genotype_plot(vcf_object  =  vcf,   
                            popmap = ca,                           
                            snp_label_size = 1000000,                     
                            colour_scheme=c(yel,"red","blue"),
                            plot_allele_frequency=F,                    
                            invariant_filter = F)    


gt<-chr16_plot2$genotypes[["data"]]
positions<-data.frame(getFIX(vcf))

gt$pos<-rep(as.integer(positions$POS), each=70)
gt_prp<-aggregate(gt$GT, by=list(gt$pos), table)

gtlist<-gt_prp$x
names(gtlist)<-gt_prp$Group.1

genotypes<-data.frame(pos=as.integer(names(gtlist)))
genotypes$homo=0
genotypes$het=0
genotypes$homo.alt=0

for (i in 1: length(gtlist)){
    df<-data.frame(gtlist[[i]])
    if (length(df$Freq[df$Var1==0])!=0) genotypes$homo[i]<-df$Freq[df$Var1==0]
    if (length(df$Freq[df$Var1==1])!=0) genotypes$het[i]<-df$Freq[df$Var1==1]
    if (length(df$Freq[df$Var1==2])!=0) genotypes$homo.alt[i]<-df$Freq[df$Var1==2]
}

genotypes$pos<-as.integer(genotypes$pos)

ge2<-genotypes[genotypes$pos<=5000000,]

library(reshape2)
gty<-melt(ge2, id.vars="pos")
ggplot(gty, aes(x=pos, y=value, color=variable, fill=variable))+
    geom_bar(stat="identity", position="fill")

gty2<-gty[gty$pos>3500000&gty$pos<4200000,]
ggplot(gty2, aes(x=pos, y=value, color=variable, fill=variable))+
    geom_bar(stat="identity", position="fill")+
    scale_x_continuous(breaks=seq(3500000,4200000, 100000))+
    theme(axis.text.x = element_text(angle=90))

