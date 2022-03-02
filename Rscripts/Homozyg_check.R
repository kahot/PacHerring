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
my_vcf <- read.vcfR("Data/vcfs/chr8_sex.vcf")
#dna<-read.dna("Data/chrY/rev_comp_chrY_v1.0.fa",format = "fasta")
#gff <- read.table("Data/chrY/chrY.gff", sep="\t", quote="")

my_vcf <- read.vcfR("Data/chrY/chr8_sex2.vcf")


meta <- read.csv("Data/Sample_metadata_892pops.csv")

my_popmap <- meta[c("Sample","Population.Year")]

#chrom<-create.chromR(name='chrY', vcf=my_vcf, seq=dna, ann=gff)
#plot(chrom)

#allelic depth
#dp <- extract.gt(my_vcf, element='DP', as.numeric=TRUE)
#boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
#abline(h=seq(0,1e4, by=100), col="#C0C0C088")
#
##violin plot
#
#dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
#dpf <- dpf[ dpf$Depth > 0,]

#ggplot(dpf, aes(x=Sample, y=Depth)) + 
#    geom_violin(fill="#C0C0C0", adjust=1.0,scale = "count", trim=TRUE)+
#    theme_bw()+
#    theme(axis.title.x = element_blank(), 
#        axis.text.x = element_text(angle = 90, hjust = 1, size=6))+
#    stat_summary(fun.y=median, geom="point", shape=23, size=2)


ad <- extract.gt(my_vcf, element = "AD")

AD_frequency(ad, allele = 1L, sum_type = 0L, decreasing = 1L)

#read the variant informatio nfor chrY
chrY_vars<-read.csv("Data/chrY/pnas.2009925117.sd01.csv")

poss<-getPOS(my_vcf)

sex.snps<-intersect(chrY_vars$Pos, poss)

# plot variation for individuals grouped by population
new_plot <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = my_popmap,                              
                          cluster        = TRUE,                           
                          snp_label_size = 50000,                          
                          colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
                          #missingness=0.8
                        )

pdf(paste0("Output/chr8/chr8_sexregion_snps_cluster_68SNPs.pdf"),width=20,height=10)
combine_genotype_plot(new_plot)
dev.off()

plot2 <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = my_popmap,                              
                          cluster        = F,                           
                          snp_label_size = 50000,                          
                          colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

#pdf(paste0("Output/chr8/chr8_sexregion_snps.pdf"),width=12,height=8)
combine_genotype_plot(plot2)
#dev.off()


#selected 26 positions
vcf <- read.vcfR("Data/chrY/chr8_sex_selected.recode.vcf")

plot3 <- genotype_plot(vcf_object  =  vcf,
                       popmap = my_popmap,                              
                       cluster        = T,                           
                       snp_label_size = 50000,                          
                       colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/chr8_sexregion_snps_cluster_selected.pdf"),width=12,height=8)
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
    pwsGT$homo_percetn[i]<-pwsGT$homo.loci[i]/nrow(df)*100
    pwsGT$hetero_percetn[i]<-pwsGT$hetero.loci[i]/nrow(df)*100
}


## Expand to the entire populations
all_gt<-plot3$genotypes[["data"]]
all_gt<-new_plot$genotypes[["data"]]





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

barplot(allGT$homo_percent)
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


#positions and dAF from chrY
ggplot(chrY_vars, aes(x=Pos, y=dAF))+
    geom_vline(xintercept=positions, color="lightpink", alpha=0.5)+
    geom_point(color="gray30")+
    theme_classic()+xlab('')
ggsave("Output/chr8/chrY-dAF_positions+ourSNPloci.pdf", width = 6, height = 4.5)

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

pdf(paste0("Output/chr8/chr8_sexregion_snps_assignedSex_68loci.pdf"),width=12,height=8)
combine_genotype_plot(out1)
dev.off()


out2 <- genotype_plot(vcf_object  =  vcf,
                      popmap = newpop,                              
                      cluster        = T,                           
                      snp_label_size = 50000,                          
                      colour_scheme=c("#EEF4ED","#94B0E0",'#4F48C6'), 
)

pdf(paste0("Output/chr8/chr8_sexregion_snps_assignedSex_cluster_68loci.pdf"),width=18,height=8)
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

pdf(paste0("Output/chr8/chr8_sexregion_snps_assignedSex_label_68loci.pdf"),width=24,height=16)
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
combine_genotype_plot(plot5)






# plot allele frequencies for populations
afplot <- genotype_plot(vcf_object  = vcf,   
                          popmap = newpop,                           
                          snp_label_size = 50000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=TRUE,                    
                          invariant_filter = TRUE)                      

pdf("Output/chr8/Sex_group_AF_68loci.pdf",width=12,height=9)
combine_genotype_plot(afplot)
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


# Visualize allele frequency on chromosome 
chr <- "chr12"

my_vcf <- read.vcfR(paste("Data/vcfs/chrs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_",chr,".vcf.gz", sep = ""))

meta <- read.csv('Data/vcfs/sample_metadata.txt', comment.char = '#', sep = "\t")
my_popmap <- meta[c("Sample","Population_Year")]

#my_popmap <- head(my_popmap,100)
pops <- c("PWS17","CA17")

pops <- c("TB17","PWS17","SS17","BC17","WA17","CA17")

my_popmap <-  my_popmap[my_popmap$Population_Year %in% pops,]
my_popmap <-  sample_n(my_popmap[my_popmap$Population_Year %in% pops,],50)












# plot subset of chromosome
new_plot <- genotype_plot(vcf_object  =  my_vcf,   
                          #chr    = 12,                                 
                          #start  = 11700000,                           
                          #end    = 11800000,                           
                          popmap = my_popmap,                           
                          #cluster        = TRUE,                        
                          snp_label_size = 5000000,                     
                          colour_scheme=c(yel,org,red),
                          plot_allele_frequency=TRUE,                    
                          invariant_filter = TRUE)                      

#dendro_with_tips <- new_plot$dendrogram +
#                    geom_text(aes(x=1:length(new_plot$dendro_labels),y=-2.5,label=new_plot$dendro_labels))
#dendro_with_tips

pdf("Output/inversions/test_plot.pdf",width=10,height=8)
pdf("Output/inversions/t2017_plot.pdf",width=10,height=8)

combine_genotype_plot(new_plot)

dev.off()

