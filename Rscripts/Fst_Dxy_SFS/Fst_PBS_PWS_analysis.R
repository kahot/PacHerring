#Pacific Herring Fst all comparisons chromosome view 
library(ggplot2)
library(tidyverse)
#library(reticulate)
library(reshape2)
#library(plyranges)
library(seqinr)
library(gridExtra)
library(vcfR)

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

#pop_names = c("PWS91","PWS96","PWS07","PWS17","TB91","TB96","TB06","TB17","SS96","SS06","SS17","BC17","WA17","CA17")

# read the pbs calculated for 50k windows by Joe
pb1 <- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS07_PWS17.txt")
pb2<-read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS07.txt")
pb3<- read.delim("Data/fst_pbs/maf05/fst_pbs_50kb_win_10kb_step_folded_PWS91_PWS96_PWS17.txt")

#combine the pbs values for PWS96, 07, 17 against PWS91 in one data frame

colnames(pb1)[8:10]<-c("PBS91","PBS07","PBS17")
colnames(pb2)[8:10]<-c("PBS91","PBS96","PBS07")
colnames(pb3)[8:10]<-c("PBS91","PBS96","PBS17")

#PBS values are different depending on the combinations of 3 populations
dif<-pb3$PBS96-pb2$PBS96
mean(dif) #0.0004655003
median(dif) # 0.000416
plot(dif,pch="." ) #different values 

dif2<-pb1$PBS07-pb2$PBS96
mean(dif2) # 0.003671074  -an order of magnitude larger 
median(dif2)#0.002961
par(mfrow=c(2,1))
plot(dif, pch=".")
plot(dif2, pch=".")

dev.off()

#Use the PBS96 analyzed with PBS 07
pbs<-merge(pb1, pb2[,c("region", "PBS96")], by="region")
pbs$ch<-as.integer(gsub("chr",'',pbs$chr))

pbs<-pbs[order(pbs$ch, pbs$midPos),]
pbs$loc<-1:nrow(pbs)

#odds<-paste0("chr",seq(1,26, by=2))
evens<-paste0("chr",seq(2,26, by=2))
pbs$color<-"col1"
pbs$color[pbs$chr %in% evens]<-"col2"


ggplot(pbs, aes(x=loc, y=PBS96, color=color))+
    geom_point(size=0.2)+
    scale_color_manual(values=c("gray30",grb))+
    theme_classic()+
    ylab("PBS")+xlab('')+
    theme(legend.position = "none")+
    ggtitle("PWS96")


pbsm<-melt(pbs[c("chr","loc","color","PBS96","PBS07","PBS17")], id.vars=c("chr","loc","color"))
ggplot(pbsm, aes(x=loc, y=value, color=color))+
    facet_wrap(~variable, ncol = 1, strip.position="right")+
    geom_point(size=0.2)+
    scale_color_manual(values=c("gray30",grb))+
    theme_bw()+
    ylab("PBS")+xlab('Genome position')+
    theme(legend.position = "none")+
    theme(axis.text.x = element_blank())
ggsave("Output/fst_pbs/PWS_pbs_3yrs_stacked.pdf", width = 10, height = 7)


# top 1% outliers of PBSs (Fst/PBS from ANGSD were based on windows)
pb.out<-pbsm[order(abs(pbsm$value), decreasing = T),] #210222 windows
pb.out<-pb.out[1:1913,]

pb.out<-merge(pb.out, pbs[,c("loc","region","midPos")], by="loc", all.x=T)
table(pb.out$variable)
#PBS96 PBS07 PBS17 
#268  1084   561  

length(unique(pb.out$loc))
#1913/2102

p17<-pb.out[pb.out$variable=="PBS17",]
p07<-pb.out[pb.out$variable=="PBS07",]
p96<-pb.out[pb.out$variable=="PBS96",]

sum<-data.frame(table(pb.out$variable,pb.out$chr))
#       1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
#PBS96 15  0 15 14  4 20 20 12  7 14  8 28  5 15 15  9  4 18 15  7  8 11 26 11  1  0
#PBS07 66 38 58 54 27 59 45 39 19 87 45 65 35 57 48 59 46 47 61 38 58 28 40 35 11 12
#PBS17 23 28 31 25 26 32 21 19  9 34 28 38 24 33  5 42 32 20 37 27 15  8 40  8  8 10

chr<-aggregate(sum$Freq, by=list(sum$Var2), sum)
chr[order(chr$x),]

#etract the regions and find the annotations
#create the bed file
pops<-c("PWS96","PWS07","PWS17")
files<-c("p96","p07","p17")
for (i in 1: length(pops)){
    loc<-data.frame()
    n=1
    pbs<-get(files[i])
    while (n <=nrow(pbs)){
        x<-pbs$loc[n]
        
        if (!is.na(pbs$loc[n+1])& pbs$loc[n+1]!=(x+1)) {
            newrow=c(pbs$loc[n], pbs$chr[n],pbs$midPos[n]-25000, pbs$midPos[n]+25000 )
            loc<-rbind(loc, newrow)
            n=n+1
        }
        else if (is.na(pbs$loc[n+1])) n=n+1
        else if (pbs$loc[n+1]==(x+1)){
            k=0
            while (!is.na(pbs$loc[n+k]) & pbs$loc[n+k]==(x+k)) k=k+1
            newrow=c(pbs$loc[n], pbs$ch[n],pbs$midPos[n]-25000, pbs$midPos[n+k-1]+25000)
            loc<-rbind(loc, newrow)
            n=n+k
        }
    }        
    
    loc<-loc[,-1]
    #convert the numbers to non-scientific
    loc[,2]<-as.integer(loc[,2])
    loc[,3]<-as.integer(loc[,3])
    write.table(loc, paste0("Output/fst_pbs/PWS/", pops[i],"_OutlierPBS_loci.bed"), quote=F, sep="\t",row.names=F, col.names = F)
}





# Create a new vcf file containing the loci in p07_loc.
 ## at terminal
vcftools --gzvcf Data/vcfs/population_PWS07_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS07_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS07_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS17_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS17_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS96_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS96_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS96_pbsoutlier --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS07_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS17_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS96_PBS_outlier.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS07_PBS_outlier.vcf > PWS07_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS17_PBS_outlier.vcf > PWS17_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS96_PBS_outlier.vcf > PWS96_annotation


#snpEff gene results
pops<-c("PWS96","PWS07","PWS17")
for (f in 1:3){
    #read the annotation info
    df<-read.table(paste0("Output/fst_pbs/PWS/",pops[f],"_annotation"), header = F)
    annotations<-data.frame()
    for (i in 1: nrow(df)){
        anns<-unlist(strsplit(df$V4[i], "\\|"))
        anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
        annotations<-rbind(annotations, anns)
    }     

    colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
    Ano<-cbind(df[,1:3], annotations)
    colnames(Ano)[1:3]<-c("chr","pos","AF")
    #remove the duplicated annotations for deeper digging
    remove<-!duplicated(annotations)
    Ano2<-Ano[remove,]
    write.csv(Ano2, paste0("Output/fst_pbs/PWS/", pops[f],"_highPBS_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-unique(geneids)
    
    
    sink(paste0("Output/fst_pbs/PWS/", pops[f],"_geneid_list.txt"))
    cat(paste0(geneids,"; "))
    sink(NULL)
}



annotations<-data.frame()
for (i in 1: nrow(df)){
    anns<-unlist(strsplit(df$V4[i], "\\|"))
    anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
    annotations<-rbind(annotations, anns)
}      

colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
Ano<-cbind(df[,1:3], annotations)
colnames(Ano)[1:3]<-c("chr","pos","AF")


#remove the duplicated annotations for deeper digging
remove<-!duplicated(annotations)
Ano2<-Ano[remove,]
write.csv(Ano2, "Output/fst_pbs/PWS/PWS07_highPBS_genelist.csv")            

geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
geneids<-unique(geneids)
geneids<-geneids[nchar(geneids)<=18]
geneids<-unique(geneids)


sink("Output/fst_pbs/PWS/PWS07_geneid_list.txt")
cat(paste0(geneids,"; "))
sink(NULL)



#####  from snpEff results ###
# Read the gene file

p07<-read.delim("Output/fst_pbs/PWS/PWS07.genes.txt",skip=1 )

#unique genes

#remove the no-annotated genes
p07a<-p07[grep("ENSCHA", p07$X.GeneName,invert = T),]
length(unique(p07a$X.GeneName)) #570 genes

write.table(unique(p07a$X.GeneName), "Output/fst_pbs/PWS/pws07_unique_genes.txt", quote = F, row.names = F)


## 
Ano2[Ano2$Gene_name=="FHIT",]
Ano2[Ano2$Gene_name=="styx",]

findGene<-function(gene){
    df<-Ano2[Ano2$Gene_name==gene,]
    
    return(df)
}

findGene("DDAH1")
findGene("vaspa")
findGene("hoxa2b")
Ano2[Ano2$Gene_name2=="hoxa2b",]
