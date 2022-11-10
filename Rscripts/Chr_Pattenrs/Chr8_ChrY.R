library(ggplot2)
library(tidyverse)
library(vcfR)
library(GenotypePlot)
library(ape)
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


# Visualize variation on chromosome 

#my_vcf <- read.vcfR(paste("Data/vcfs/chrs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_",chr,".vcf.gz", sep = ""))
#my_vcf <- read.vcfR("Data/vcfs/chr8_sex.vcf")

my_vcf <- read.vcfR("Data/chrY/chr8_sex2.vcf")


meta <- read.csv("Data/Sample_metadata_892pops.csv")
my_popmap <- meta[c("Sample","Population.Year")]

#ad <- extract.gt(my_vcf, element = "AD")
#AD_frequency(ad, allele = 1L, sum_type = 0L, decreasing = 1L)


#read the variant position info for chrY of Atlantic Herring
chrY_vars<-read.csv("Data/chrY/pnas.2009925117.sd01.csv")
#overlapping positions from Atlantic Herring
poss<-getPOS(my_vcf)
sex.snps<-intersect(chrY_vars$Pos, poss) #therea re 26 snps that fall in the sex chromosome region identified in Atlantic Herring

sink("Data/chrY/dAF_pos.txt") #if not using for loop, there will be a space at the beginning from the second row on. Don't know why.
for (i in 1: length(sex.snps)){
    cat(paste0("chr8","\t",sex.snps[i],"\n"))
}    
sink(NULL)


# plot variation for individuals grouped by population
new_plot <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = my_popmap,                              
                          cluster        = TRUE,                           
                          snp_label_size = 50000,                          
                          colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
                          #missingness=0.8
                        )

pdf(paste0("Output/chr8/chr8_sexregion_snps_cluster_67SNPs.pdf"),width=20,height=10)
combine_genotype_plot(new_plot)
dev.off()

plot2 <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = my_popmap,                              
                          cluster        = F,                           
                          snp_label_size = 50000,                          
                          colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/chr8_sexregion_snps_byPopGroup_67SNPs.pdf"),width=12,height=8)
combine_genotype_plot(plot2)
dev.off()


#selected 26 positions using dAF_pos.txt
vcf <- read.vcfR("Data/chrY/chr8_sex_selected.recode.vcf")

plot3 <- genotype_plot(vcf_object  =  vcf,
                       popmap = my_popmap,                              
                       cluster        = T,                           
                       snp_label_size = 50000,                          
                       colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Selected_chr8_sexregion_cluster_26loci.pdf"),width=12,height=8)
combine_genotype_plot(plot3)
dev.off()

### group individuals based on GT (homo/hetero)
#Start with PWS
pws<-my_popmap[grep("PWS", my_popmap$Sample),]
pws1 <- genotype_plot(vcf_object  =  vcf,
                       popmap = pws,                              
                       cluster        = T,                           
                       snp_label_size = 50000, colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6') )
#combine_genotype_plot(pws1)

#extract the genotype data
pws_gt<-pws1$genotypes[["data"]]

#Add the sample id
pws_gt$Sample<-pws$Sample

#first pull out all homozygous individuals

pwsGT<-data.frame(Sample=pws$Sample)
for (i in 1:length(pws$Sample)){
    id<-pws$Sample[i]
    df<-pws_gt[pws_gt$Sample==id,]
    #remove missing locus
    df<-df[!is.na(df$GT),]
    #count homozygous sites
    pwsGT$homo.loci[i]<-nrow(df[df$GT==0,])
    pwsGT$hetero.loci[i]<-nrow(df[df$GT==2,])
    pwsGT$homo_percent[i]<-pwsGT$homo.loci[i]/nrow(df)*100
    pwsGT$hetero_percent[i]<-pwsGT$hetero.loci[i]/nrow(df)*100
}

#proportion of all homozygous individuals per population
pwsGT$homo<-FALSE
pwsGT$homo[pwsGT$homo_percent==100]<-TRUE
pwsGT$pop<-substr(pwsGT$Sample, 6,10)
prop.table(table(pwsGT$homo, pwsGT$pop), 2)
#            PWS07     PWS17     PWS91     PWS96
#FALSE 0.7826087 0.5357143 0.7586207 0.7777778
#TRUE  0.2173913 0.4642857 0.2413793 0.2222222


### Expand to the entire populations
all_gt<-plot3$genotypes[["data"]]
#all_gt<-new_plot$genotypes[["data"]]

#Add the sample id
all_gt$Sample<-my_popmap$Sample

#first pull out all homozygous individuals

allGT<-data.frame(Sample=my_popmap$Sample)
for (i in 1:length(my_popmap$Sample)){
    id<-my_popmap$Sample[i]
    df<-all_gt[all_gt$Sample==id,]
    #remove missing locus
    df<-df[!is.na(df$GT),]
    #count homozygous sites
    allGT$homo.loci[i]<-nrow(df[df$GT==0,])
    allGT$alt.loci[i]<-nrow(df[df$GT==2,])
    allGT$homo_percent[i]  <-allGT$homo.loci[i]/nrow(df)*100
    allGT$alt_percent[i]<-allGT$alt.loci[i]/nrow(df)*100
    allGT$hetero_percent[i]<-nrow(df[df$GT==1,])/nrow(df)*100
}

allGT$alt.and.hetero<-allGT$alt_percent+allGT$hetero_percent

barplot(allGT$alt.and.hetero)

plot(allGT$homo_percent, pch=".")
plot(allGT$alt_percent, pch=".")
plot(allGT$alt.and.hetero, pch=".")


#how many have all homozygous?
nrow(allGT[allGT$homo_percent==100,])
#261
nrow(allGT[allGT$homo_percent>=80,])
#431

nrow(allGT[allGT$alt.and.hetero>=50,])
#161
nrow(allGT[allGT$alt.and.hetero>=30,])
#376

#Let's call >-80% as female
allGT$sex<-NA
#likely female
allGT$sex[allGT$homo_percent>=80]<-"F"
#pretty sure female
allGT$sex[allGT$homo_percent==100]<-"FF"

#likely male
allGT$sex[allGT$alt.and.hetero>=30]<-"M"

allGT$sex[is.na(allGT$sex)]<-"M_likely"

#Last 2 loci seems unimportant in determination
positions<-getPOS(vcf)
chrY_vars[chrY_vars$Pos==positions[26],]


library(showtext)
showtext_auto()

female = intToUtf8(9792)
male = intToUtf8(9794)

#Atlantic Herring positions and dAF of chrY. Add the SNP locations from PacHerring data
ggplot(chrY_vars, aes(x=Pos, y=dAF))+
    geom_vline(xintercept=positions, color="pink", alpha=0.7)+
    geom_point(color="blue", size=1)+
    theme_classic()+xlab('')+ylab(paste0("dAF between ",female," & ",male))
ggsave("Output/chr8/Atlnatic_herring_chrY-dAF_positions+ourSNPloci.pdf", width = 6, height = 4.5)



#our marker locations
daf<-chrY_vars[chrY_vars$Pos %in% positions,]
mean(daf$dAF)
#0.4753846
write.csv(allGT,"Output/chr8/Individual_sex_from_genotype.csv")

#now group the populations in to the sex categories
newpop<-allGT[,c("Sample", "sex")]
newpop$sex<-factor(newpop$sex, levels=c("FF","F","M","M_likely"))
newpop<-newpop[order(newpop$sex),]
out1 <- genotype_plot(vcf_object  =  vcf,
                       popmap = newpop,                              
                       cluster        = F,                           
                       snp_label_size = 50000,                          
                       colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Selected_chr8_sexregion_assignedSex_26loci.pdf"),width=12,height=8)
combine_genotype_plot(out1)
dev.off()

#plot this on 67 loci data
out1_67 <- genotype_plot(vcf_object  =  my_vcf,
                      popmap = newpop,                              
                      cluster        = F,                           
                      snp_label_size = 50000,                          
                      colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Chr8_sexregion_assignedSexfrom26loci_67loci.pdf"),width=12,height=8)
combine_genotype_plot(out1_67)
dev.off()



out2 <- genotype_plot(vcf_object  =  vcf,
                      popmap = newpop,                              
                      cluster        = T,                           
                      snp_label_size = 50000,                          
                      colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Selected_chr8_sexregion_assignedSex_cluster_26loci.pdf"),width=18,height=8)
combine_genotype_plot(out2)
dev.off()


dendrolabel<-data.frame(Sample=out2$dendro_labels)
dendrolabel$order<-1:892
dendrolabel<-merge(dendrolabel, newpop, by="Sample" )
dendrolabel<-dendrolabel[order(dendrolabel$order),]
dendrolabel$color<-apply(dendrolabel["sex"],1, function(x) {if (x["sex"]=="FF") "red"
                        else if  (x["sex"]=="F") "pink"
                        else if (x["sex"]=="M") "blue"
                        else if (x["sex"]=="M_likely") "lightblue"})


dendro_with_tips <- out2$dendrogram +
        geom_text(aes(x=1:length(out2$dendro_labels),
                  y=-0.5,
                  label=dendrolabel$sex), size=2, color=dendrolabel$color)
out2$dendrogram<-dendro_with_tips

pdf(paste0("Output/chr8/Selected_chr8_sexregion_snps_assignedSex_label_26loci.pdf"),width=28,height=16)
combine_genotype_plot(out2)
dev.off()

######### do the same for all 67 loci 
#remove the last few positions?
pos67<-getPOS(my_vcf)
gt<-new_plot$genotypes[["data"]]
#Add the sample id
gt$Sample<-my_popmap$Sample

allGT67<-data.frame(Sample=my_popmap$Sample)
for (i in 1:length(my_popmap$Sample)){
    id<-my_popmap$Sample[i]
    df<-gt[gt$Sample==id,]
    #remove missing locus
    df<-df[!is.na(df$GT),]
    #count homozygous sites
    allGT67$homo.loci[i]<-nrow(df[df$GT==0,])
    allGT67$alt.loci[i]<-nrow(df[df$GT==2,])
    allGT67$homo_percent[i]  <-allGT67$homo.loci[i]/nrow(df)*100
    allGT67$alt_percent[i]<-allGT67$alt.loci[i]/nrow(df)*100
    allGT67$hetero_percent[i]<-nrow(df[df$GT==1,])/nrow(df)*100
}

allGT67$alt.and.hetero<-allGT67$alt_percent+allGT67$hetero_percent

plot(allGT67$homo_percent,   pch=16, col=blu)
plot(allGT67$alt_percent,    pch=16, col=blu)
plot(allGT67$alt.and.hetero, pch=16, col=blu)

#how many have all homozygous?
nrow(allGT67[allGT67$homo_percent==100,])
#62
nrow(allGT67[allGT67$homo_percent>=80,])
#475
nrow(allGT67[allGT67$homo_percent>=90,])
#316

nrow(allGT67[allGT67$alt.and.hetero>=50,])
#25
nrow(allGT67[allGT67$alt.and.hetero>=30,])
#241
nrow(allGT67[allGT67$alt_percent>=20,])
#99


#assing the same sex categories
#Let's call >-80% as female
allGT67$sex<-NA
#likely female
allGT67$sex[allGT67$homo_percent>=80]<-"F"
#pretty sure female
allGT67$sex[allGT67$homo_percent==100]<-"FF"

#likely male
allGT67$sex[allGT67$alt.and.hetero>=30]<-"M"

allGT67$sex[is.na(allGT67$sex)]<-"M_likely"

#now group the populations in to the sex categories
newpop2<-allGT67[,c("Sample", "sex")]
newpop2$sex<-factor(newpop2$sex, levels=c("FF","F","M","M_likely"))
newpop2<-newpop2[order(newpop2$sex),]
out1_67 <- genotype_plot(vcf_object  =  my_vcf,
                      popmap = newpop2,                              
                      cluster        = F,                           
                      snp_label_size = 50000,                          
                      colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Sexregion_assignedSex_67loci.pdf"),width=12,height=8)
combine_genotype_plot(out1_67)
dev.off()

#plot this on 67 loci data
out1_orig <- genotype_plot(vcf_object  =  my_vcf,
                         popmap = newpop,                              
                         cluster        = F,                           
                         snp_label_size = 50000,                          
                         colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Chr8_sexregion_assignedSexfrom26loci_67loci.pdf"),width=12,height=8)
combine_genotype_plot(out1_orig)
dev.off()



out2 <- genotype_plot(vcf_object  =  vcf,
                      popmap = newpop,                              
                      cluster        = T,                           
                      snp_label_size = 50000,                          
                      colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/Selected_chr8_sexregion_assignedSex_cluster_26loci.pdf"),width=18,height=8)
combine_genotype_plot(out2)
dev.off()






#Plot PWS only?
pwspop<-newpop[grep("PWS", newpop$Sample),]

plot5 <- genotype_plot(vcf_object  =  vcf,
                       popmap = pwspop,                              
                       cluster        = F,                           
                       snp_label_size = 50000,                          
                       colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/PWS_Selected_chr8_sexregion_snps_assignedSex_26loci.pdf"),width=18,height=10)
combine_genotype_plot(plot5)
dev.off()



# plot allele frequencies for populations
afplot <- genotype_plot(vcf_object  = vcf,   
                          popmap = newpop,                           
                          snp_label_size = 50000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=TRUE,                    
                          invariant_filter = TRUE)                      

pdf("Output/chr8/Sex_group_AF_26loci.pdf",width=12,height=9)
combine_genotype_plot(afplot)
dev.off()

# plot allele frequencies for populations
afplot2 <- genotype_plot(vcf_object  = my_vcf,   
                        popmap = newpop2,                           
                        snp_label_size = 50000,                     
                        colour_scheme=c(yel,org,red),
                        plot_allele_frequency=TRUE,                    
                        invariant_filter = TRUE)                      

pdf("Output/chr8/Sex_group_AF_67loci.pdf",width=12,height=9)
combine_genotype_plot(afplot2)
dev.off()

#### look within Populations

pwspop<-merge(pwspop, meta[,c("Sample","pop","Year.Collected")], by="Sample")

sexratio<-table(pwspop$sex, pwspop$Year.Collected)
sexr<-prop.table(sexratio, 2)

sexdf<-data.frame(table(pwspop$sex, pwspop$Year.Collected))
colnames(sexdf)<-c("Sex","Year","Freq")

ggplot(sexdf,aes(x=Year, y=Freq, fill=Sex))+
    geom_bar(stat="identity", position = 'fill', alpha=0.8)+
    scale_fill_manual(values=c("red","pink", "blue","lightblue"))+
    theme_classic()+
    ggtitle("PWS Population Sex Ratio")+
    ylab('')
    
ggsave("Output/chr8/PWS_sex_ratio.pdf", width = 5.5, height = 4)

# do the same for using 67 loci
pwspop2<-newpop2[grep("PWS", newpop2$Sample),]
pwspop2<-merge(pwspop2, meta[,c("Sample","pop","Year.Collected")], by="Sample")
sexratio<-table(pwspop2$sex, pwspop2$Year.Collected)
sexr<-prop.table(sexratio, 2)

sexdf2<-data.frame(table(pwspop2$sex, pwspop2$Year.Collected))
colnames(sexdf2)<-c("Sex","Year","Freq")

ggplot(sexdf2,aes(x=Year, y=Freq, fill=Sex))+
    geom_bar(stat="identity", position = 'fill', alpha=0.8)+
    scale_fill_manual(values=c("red","pink", "blue","lightblue"))+
    theme_classic()+
    ggtitle("PWS Population Sex Ratio")+
    ylab('')
ggsave("Output/chr8/PWS_sex_ratio_67loci.pdf", width = 5.5, height = 4)

