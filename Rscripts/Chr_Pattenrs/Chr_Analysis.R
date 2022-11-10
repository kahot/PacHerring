#### Plot which individuals have freq changes in this region
library(ggplot2)
library(tidyverse)
library(vcfR)
#install.packages("remotes")
#remotes::install_github("JimWhiting91/genotype_plot")
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


#dir.create("Output/chr/")
#dir.create("Output/chr/DP7000")


pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[c("Sample","Population.Year")]
pops<-unique(pop_info$Population.Year)


#Specify the chromosome
ch<-"chr20"

dir.create(paste0("Output/chr/DP7000/",ch))

vcf <- read.vcfR(paste0("Data/new_vcf/PH_",ch,".vcf.gz"))

vcf <- read.vcfR(paste0("Data/new_vcf/PWSonly/PWSonly_chr8/PWSonly_ch8_maf05.vcf"))
vcf <- read.vcfR(paste0("Data/new_vcf/PWS/PWS_ch8_maf05.vcf"))
vcf <- read.vcfR(paste0("Data/new_vcf/PWSonly/PWSonly_chr4/PWSonly_ch4_maf05.vcf"))
vcf <- read.vcfR(paste0("Data/new_vcf/PWSonly/PWSonly_chr20/PWSonly_ch20_maf05.vcf"))

#PWS population 
pws<-c("PWS91","PWS96","PWS07","PWS17")
#1.
p_map<-pop_info[pop_info$Population.Year==pws[1],]
chr_plot1 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                              snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                              plot_allele_frequency=F,invariant_filter = F)  
    
gt_plot<-chr_plot1$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                       breaks=c("0","1","2"),guide="none")
chr_plot1$genotypes<-updated
combine_genotype_plot(chr_plot1)

#2
p_map<-pop_info[pop_info$Population.Year==pws[2],]
chr_plot2 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot2$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot2$genotypes<-updated
combine_genotype_plot(chr_plot2)
#3
p_map<-pop_info[pop_info$Population.Year==pws[3],]
chr_plot3 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot3$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot3$genotypes<-updated    
combine_genotype_plot(chr_plot3)

#4
p_map<-pop_info[pop_info$Population.Year==pws[4],]
chr_plot4 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot4$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
chr_plot4$genotypes<-updated    
combine_genotype_plot(chr_plot4)


pdf(paste0("Output/chr/DP7000/",ch,"/",ch,"_PWS-all.pdf"), width =8, height = 20 )
ggdraw()+
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.73,width=1,height=0.27)+
    draw_plot(chr_plot2$genotypes,0,0.5,1,0.23)+
    draw_plot(chr_plot3$genotypes,0,0.27,1,0.23)+
    draw_plot(chr_plot4$genotypes,0,0,1,0.27)
dev.off()


combine_genotype_plot(chr_plot3)
combine_genotype_plot(chr_plot4)



### chr20 groups ###
pw20<-read.csv("Output/PCA/PWSonly_with.groups_pca_chr20.csv", row.names = 1)
pw20<-pw20[,c("Sample","group")]

# Group1.
p_map<-pw20[pw20$group=="Group1",]
chr_plot1 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  

gt_plot<-chr_plot1$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot1$genotypes<-updated
combine_genotype_plot(chr_plot1)

# Group2
p_map<-pw20[pw20$group=="Group2",]
chr_plot2 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, 
                           snp_label_size = 1000000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot2$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),guide="none")
chr_plot2$genotypes<-updated
combine_genotype_plot(chr_plot2)


#Group3
p_map<-pw20[pw20$group=="Group3",]
chr_plot3 <- genotype_plot(vcf_object  =  vcf,    popmap = p_map, start=20000000, end=28000000,
                           snp_label_size = 200000,colour_scheme=c(yel,"red","blue"),                  
                           plot_allele_frequency=F,invariant_filter = F)  
gt_plot<-chr_plot3$genotypes   
updated<-gt_plot+scale_fill_manual(values=c(yel,"red","blue"),name="Genotype",na.value = "gray90",
                                   breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))
chr_plot3$genotypes<-updated    
combine_genotype_plot(chr_plot3)

png(paste0("Output/PCA/chr20/PWSonlt_ch20_groups.png"), width =8, height = 20, res=300)
ggdraw()+
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.66,width=1,height=0.34)+
    draw_plot(chr_plot2$genotypes,0,0.33,1,0.31)+
    draw_plot(chr_plot3$genotypes,0,0,1,0.35)
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
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.63,width=1,height=0.37)+
    draw_plot(chr_plot2$genotypes,0,0.37,1,0.28)+
    draw_plot(chr_plot3$genotypes,0,0,1,0.37)
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
    draw_plot(combine_genotype_plot(chr_plot1),x=0,y=0.63,width=1,height=0.37)+
    draw_plot(chr_plot2$genotypes,0,0.37,1,0.28)+
    draw_plot(chr_plot3$genotypes,0,0,1,0.37)
dev.off()


#### Genotype info per individual
for (i in 1:length(pops)){
    pop<-pops[i]
    p_map<-pop_info[pop_info$Population.Year==pop,]
    
    chr_plot <- genotype_plot(vcf_object  =  vcf,   
                              popmap = p_map,                           
                              #start=14000000,
                              #end=30000000,
                              snp_label_size = 1000000,                     
                              colour_scheme=c(yel,"red","blue"),
                              plot_allele_frequency=F,                    
                              invariant_filter = F)  
    
    gt<-chr_plot$genotypes[["data"]]
    gt$Sample<-p_map$Sample
    write.csv(gt, paste0("Output/chr/DP7000/",ch,"/", pops[i], "_gtinfo.csv"))
    
}
    
    

## Chr8  What are the 3 groups?
pca<-read.csv("Output/PCA/noTB_pca_chr8.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1<=-0.003, ]
gp3.1<-pca[pca$PC1>0.07&pca$pop!="CA", ]
gp3.2<-pca[pca$PC1>0.13&pca$pop=="CA", ]
gp3<-rbind(gp3.1, gp3.2) #29
sampl<-c(gp1$Sample, gp3$Sample)
gp2<-pca[!(pca$Sample %in% sampl),  ] #176


table(gp1$yr.pop)
# BC17  CA17 PWS07 PWS17 PWS91 PWS96  SS06  SS17  SS96  WA17 
#   51    27    27    25    37    51    30    46    58    55 

table(gp2$yr.pop)
#BC17  CA17 PWS07 PWS17 PWS91 PWS96  SS06  SS17  SS96  WA17 
#  12    31    17    28    19    19    11    17    17    14 

table(gp3$yr.pop)
#BC17  CA17 PWS07 PWS17 PWS91 PWS96  SS17  SS96  WA17 
#   1    12     2     3     2     2     1     3     3 


gp1$Group<-1
gp2$Group<-2
gp3$Group<-3

gp<-rbind(gp1[,c("Group","Sample","pop","year","yr.pop")],gp2[,c("Group","Sample","pop","year","yr.pop")],gp3[,c("Group","Sample","pop","year","yr.pop")])
write.csv(gp, "Output/chr/DP7000/chr8/chr8_PCAgroups_updated.csv")
#pws17 and CA17 has higher proportions of group 1 and 2


gp1_tb<-data.frame(table(gp1$yr.pop))
gp2_tb<-data.frame(table(gp2$yr.pop))
gp3_tb<-data.frame(table(gp3$yr.pop))

gp1_tb$Group<-1
gp2_tb$Group<-2
gp3_tb$Group<-3


group<-rbind(gp1_tb,gp2_tb, gp3_tb)

group$Var1<-factor(group$Var1, levels = c("PWS91","PWS96","PWS07","PWS17","SS96","SS06", "SS17","BC17","WA17","CA17"  ))
ggplot(group, aes(x=Var1, y=Freq, fill=factor(Group)))+
    geom_bar(stat="identity", position="fill")+
    scale_fill_manual(values=c(red, gre, blu), name="Group")+
    xlab("")+ylab("Proportion ")+ggtitle("Chr8")
ggsave("Output/chr/DP7000/chr8/Chr8_3group_proportion.pdf", width = 8, height = 5.5)


#individulas in group 1
group1<-gp[gp$Group==1, ]
ch8.group<-gp1

#chr12
pca<-read.csv("Output/PCA/noTB_pca_chr12.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1>0.1, ] #45
gp3<-pca[pca$PC1<0.02, ] #505
sampl<-c(gp1$Sample, gp3$Sample)
gp2<-pca[!(pca$Sample %in% sampl),  ] #71

gp1$Sample[gp1$pop=="CA"] #44
ch8.group$Sample[ch8.group$pop=="CA"] #12
intersect(gp1$Sample[gp1$pop=="CA"],ch8.group$Sample[ch8.group$pop=="CA"])
#"0551_CA17" "0681_CA17" "0772_CA17" "0773_CA17" "0832_CA17" "1246_CA17"
#Only 6 individuals intersect
gp1$Group<-1
gp2$Group<-2
gp3$Group<-3
gp<-rbind(gp1[,c("Group","Sample","pop","year","yr.pop")],gp2[,c("Group","Sample","pop","year","yr.pop")],gp3[,c("Group","Sample","pop","year","yr.pop")])

write.csv(gp, "Output/chr/DP7000/chr12_PCAgroups.csv")


## chr16
pca<-read.csv("Output/PCA/noTB_pca_chr16.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1>0.15, ] #15
gp3<-pca[pca$PC1<0, ] #523
sampl<-c(gp1$Sample, gp3$Sample)
gp2<-pca[!(pca$Sample %in% sampl),  ] #83

intersect(gp1$Sample[gp1$pop=="CA"],ch8.group$Sample[ch8.group$pop=="CA"])
#2 individuals
#[1] "0551_CA17" "0844_CA17"

gp1$Group<-1
gp2$Group<-2
gp3$Group<-3
gp<-rbind(gp1[,c("Group","Sample","pop","year","yr.pop")],gp2[,c("Group","Sample","pop","year","yr.pop")],gp3[,c("Group","Sample","pop","year","yr.pop")])

write.csv(gp, "Output/chr/DP7000/chr16_PCAgroups.csv")

gp16<-read.csv("Output/chr/DP7000/chr16_PCAgroups.csv", row.names = 1)
gp12<-read.csv("Output/chr/DP7000/chr12_PCAgroups.csv", row.names = 1)

intersect(gp16$Sample[gp16$Group==1], gp12$Sample[gp12$Group==1])
#"0551_CA17" "0771_CA17" "0775_CA17" "0846_CA17" "0976_CA17" "1233_CA17" "1234_CA17" "1242_CA17"
#"1247_CA17"

#chr7
pca<-read.csv("Output/PCA/noTB_pca_chr7.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1>0.1, ] #40
gp3<-pca[pca$PC1<0, ] #518
sampl<-c(gp1$Sample, gp3$Sample)
gp2<-pca[!(pca$Sample %in% sampl),  ] #63
gp1$Group<-1
gp2$Group<-2
gp3$Group<-3
gp<-rbind(gp1[,c("Group","Sample","pop","year","yr.pop")],gp2[,c("Group","Sample","pop","year","yr.pop")],gp3[,c("Group","Sample","pop","year","yr.pop")])
write.csv(gp, "Output/chr/DP7000/chr7_PCAgroups.csv")
gp7<-gp

intersect(gp7$Sample[gp7$Group==1], gp12$Sample[gp12$Group==1])
#[1] "0548_CA17" "0549_CA17" "0551_CA17" "0683_CA17" "0685_CA17" "0686_CA17" "0769_CA17" "0773_CA17"
#[9] "0828_CA17" "0829_CA17" "0830_CA17" "0832_CA17" "0842_CA17" "0845_CA17" "0847_CA17" "0848_CA17"
#[17] "0913_CA17" "0917_CA17" "0969_CA17" "0970_CA17" "0972_CA17" "0973_CA17" "0976_CA17" "1247_CA17"

intersect(gp16$Sample[gp16$Group==1], gp7$Sample[gp7$Group==1])


#samples with the inversion
inv7<-read.csv("Output/chr7/CA_samples_with_inversions_chr7.csv", row.names = 1)
intersect(gp7$Sample[gp7$Group==1], inv7$Sample)
#all samples with inversions (CA 28) +  "0435_PWS07" "1217_WA17" are in Group1

sink("Output/chr/DP7000/Samples_with.chr72.inversions.txt")
cat(gp7$Sample[gp7$Group==1])
sink(NULL)
# [1] "0545_CA17" "0547_CA17" "0548_CA17" "0549_CA17" "0551_CA17" "0683_CA17" "0684_CA17" "0685_CA17"
#[9] "0686_CA17" "0687_CA17" "0688_CA17" "0769_CA17" "0773_CA17" "0828_CA17" "0829_CA17" "0830_CA17"
#[17] "0832_CA17" "0841_CA17" "0842_CA17" "0845_CA17" "0847_CA17" "0848_CA17" "0913_CA17" "0914_CA17"
#[25] "0917_CA17" "0918_CA17" "0919_CA17" "0920_CA17" "0969_CA17" "0970_CA17" "0971_CA17" "0972_CA17"
#[33] "0973_CA17" "0975_CA17" "0976_CA17" "1236_CA17" "1239_CA17" "1247_CA17"


inv12<-read.csv("Output/AF/chr12/CA_samples_with_inversions.csv", row.names = 1)
gp12$Sample[gp12$Group==1]

intersect(inv12$Sample, gp12$Sample[gp12$Group==1])
sink("Output/chr/DP7000/Samples_with.chr12.inversions.txt")
cat(gp12$Sample[gp12$Group==1])
sink(NULL)

sink("Output/chr/DP7000/Samples_with.chr12.chr7.inversions.txt")
cat(intersect(gp7$Sample[gp7$Group==1], gp12$Sample[gp12$Group==1]))
sink(NULL)
#24 have both inversions

#look at chr16 inversions
inv16<-read.csv("Output/chr/DP7000/chr16/CA_samples_with_inversions_chr16.csv", row.names = 1)
intersect(inv16$Sample,gp16$Sample[gp16$Group==1] ) #13
# "0551_CA17" "0771_CA17" "0775_CA17" "0844_CA17" "0846_CA17" "0971_CA17" "0976_CA17" "1233_CA17"
# "1234_CA17" "1236_CA17" "1242_CA17" "1247_CA17" "1248_CA17"
sink("Output/chr/DP7000/Samples_with.chr16.inversion.txt")
cat(gp16$Sample[gp16$Group==1])
sink(NULL)


sink("Output/chr12/Samples_with.chr12.chr7.inversions.txt")
cat(intersect(gp7$Sample[gp7$Group==1], intersect(gp12$Sample[gp12$Group==1], gp16$Sample[gp16$Group==1])))
sink(NULL)
#3 individuals have all 3 inversions



## chr20
pca<-read.csv("Output/PCA/noTB_pca_chr20.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1>0.1&pca$PC2>-0.05, ] #32
gp3.1<-pca[pca$PC1<0.025&pca$PC2<0.045, ] #506
#remove 3 samples
pca[pca$yr.pop=="SS06"&pca$PC1>0.01&pca$PC1<0.02&pca$PC2<0,] #0511_SS06
pca[pca$yr.pop=="PWS07"&pca$PC1>0.015&pca$PC1<0.025&pca$PC2>0&pca$PC2<0.023,] # 1146_PWS07
pca[pca$yr.pop=="WA17"&pca$PC1>0.015&pca$PC1<0.025&pca$PC2>=-0.01&pca$PC2<0.01,]  # 0796_WA17
gp3<-gp3.1[!(gp3.1$Sample %in% c("0511_SS06","1146_PWS07","0796_WA17")),] #503

sampl<-c(gp1$Sample, gp3$Sample)
gp2<-pca[!(pca$Sample %in% sampl),  ] #86

gp1$Group<-1
gp2$Group<-2
gp3$Group<-3
gp<-rbind(gp1[,c("Group","Sample","pop","year","yr.pop")],gp2[,c("Group","Sample","pop","year","yr.pop")],gp3[,c("Group","Sample","pop","year","yr.pop")])

write.csv(gp, "Output/chr/DP7000/chr20_PCAgroups.csv")

## ch15
pca<-read.csv("Output/PCA/noTB_pca_chr15.csv")
pca$yr.pop<-paste0(pca$pop,substr(pca$year, 3,4))
gp1<-pca[pca$PC1>0, ] #32
gp3.1<-pca[pca$PC1< (-0.035)&pca$PC2>0, ] #506
gp3.2<-pca[pca$PC1< (-0.075)&pca$PC2<0, ] 
all<-c(gp1$Sample, gp3.1$Sample,gp3.2$Sample)
gp2<-pca[!(pca$Sample %in% all),]

gp1$Group="Group1"
gp1$Group.sub="Group1"
gp3.1$Group="Group3"
gp3.1$Group.sub="Group3.1"
gp3.2$Group="Group3"
gp3.2$Group.sub="Group3.2"
gp2$Group="Group2"
gp2$Group.sub="Group2"

ch15<-rbind(gp1,gp2, gp3.1,gp3.2)
write.csv(ch15, "Output/PCA/chr15_PCAgroups.csv")
write.table(ch15$Sample[ch15$Group=="Group1"],"Data/ch15_group1.txt", quote = F, row.names = F, col.names = F)
write.table(ch15$Sample[ch15$Group=="Group2"],"Data/ch15_group2.txt", quote = F, row.names = F, col.names = F)
write.table(ch15$Sample[ch15$Group=="Group3"],"Data/ch15_group3.txt", quote = F, row.names = F, col.names = F)

#####
gp20<-gp
inv20<-read.csv("Output/chr/DP7000/chr20/CA_samples_with_inversions_chr20.csv", row.names=1)
intersect(inv20$Sample, gp20$Sample) #all 30 in group1
sink("Output/chr/DP7000/Samples_with.chr20.inversions.txt")
cat(gp20$Sample[gp20$Group==1])
sink(NULL)


gp16<-read.csv("Output/chr/DP7000/chr16_PCAgroups.csv", row.names = 1)
gp12<-read.csv("Output/chr/DP7000/chr12_PCAgroups.csv", row.names = 1)
gp7<-read.csv("Output/chr/DP7000/chr7_PCAgroups.csv", row.names = 1)


intersect(gp20$Sample[gp16$Group==1], gp12$Sample[gp12$Group==1])
# [1] "0546_CA17" "0552_CA17" "0681_CA17" "0685_CA17" "0772_CA17" "0773_CA17" "0774_CA17" "0829_CA17"
# [9] "0830_CA17" "0843_CA17"

inv3<-intersect(gp7$Sample[gp7$Group==1], intersect(gp12$Sample[gp12$Group==1], gp16$Sample[gp16$Group==1]))
intersect(inv3,gp20$Sample[gp16$Group==1]) #0
intersect(gp7$Sample[gp7$Group==1], intersect(gp12$Sample[gp12$Group==1], gp20$Sample[gp16$Group==1]))
#[1] "0685_CA17" "0773_CA17" "0829_CA17" "0830_CA17"

#create a venn diagram

library(ggvenn)
chr7<-gp7$Sample[gp7$Group==1]
chr12<-gp12$Sample[gp12$Group==1]
chr16<-gp16$Sample[gp16$Group==1]
chr20<-gp20$Sample[gp20$Group==1]
x<-list(chr7=chr7,chr12=chr12,chr16=chr16,chr20=chr20)

library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

ggvenn(
    x, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
)

library("ggVennDiagram")
ggVennDiagram(x, label_alpha = 0)+ ggplot2::scale_fill_gradient(low="white",high = "#0073C2FF")
ggsave("Output/chr/DP7000/venn.diagram_inversion_overlap.pdf", width = 6.5, height = 5)

i1<-intersect(chr7,chr20)
i2<-intersect(chr12,chr16)
intersect(i1,i2)
#"0976_CA17" "1247_CA17"

gp8<-read.csv("Output/chr/DP7000/chr8_PCAgroups.csv", row.names = 1)
gp8$Sample[gp8$Group==1]



### create proportion of grouping in each pop
library(dplyr)

#1. chr12
ch="chr12"

chs<-c("chr7","chr12","chr16","chr20")
for (i in 1:4){
    df<-read.csv(paste0("Output/chr/DP7000/",chs[i],"_PCAgroups.csv"), row.names = 1)
    pop.sum<-df %>% count(Group, yr.pop)
    pop.sum$yr.pop<-factor(pop.sum$yr.pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))
    ggplot(data=pop.sum, aes(x=yr.pop, y=n, fill=factor(Group)))+
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=c("red", org, blu), labels=c("Group1","Group2","Group3"))+
        xlab("")+ylab("Proportion")+theme(legend.title = element_blank())+
        ggtitle(chs[i])
    ggsave(paste0("Output/chr/DP7000/",chs[i],"_group_barplot.pdf"), width = 8, height = 5.5)
    
}

pop.sum<-gp12 %>% count(Group, yr.pop)

pop.sum$yr.pop<-factor(pop.sum$yr.pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))
ggplot(data=pop.sum, aes(x=yr.pop, y=n, fill=factor(Group)))+
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=c("red", org, blu), labels=c("Group1","Group2","Group3"))+
    xlab("")+ylab("Proportion")+theme(legend.title = element_blank())+
    ggtitle(ch)
ggsave("Output/chr/DP7000/chr8_group_barplot.pdf", width = 8, height = 5.5)



