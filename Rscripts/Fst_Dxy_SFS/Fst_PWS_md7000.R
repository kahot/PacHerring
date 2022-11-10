#Pacific Herring Fst all comparisons chromosome view 
source("Rscripts/BaseScripts.R")
#library(seqinr)
library(gridExtra)
library(vcfR)

# read the pbs calculated for 50k windows by Joe
pb1 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_PWS91_PWS07_PWS17_maf00_50kWindow")
pb2 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_PWS91_PWS96_PWS17_maf00_50kWindow")
pb3 <- read.delim("Data/new_vcf/angsd/fromVCF/3D/fst_PWS96_PWS07_PWS17_maf00_50kWindow")

#combine the pbs values for PWS96, 07, 17 against PWS91 in one data frame
colnames(pb1)[8:10]<-c("PBS91","PBS07","PBS17")
colnames(pb2)[8:10]<-c("PBS91","PBS96","PBS17")
colnames(pb3)[8:10]<-c("PBS96","PBS07","PBS17")

colnames(pb1)[5:7]<-c("Fst91.07","Fst91.17","Fst07.17")
colnames(pb2)[5:7]<-c("Fst91.96","Fst91.17","Fst96.17")
colnames(pb3)[5:7]<-c("Fst96.07","Fst96.17","Fst07.17")

#pb2 is missing one region. Align all 3 data frames
pb1<-merge(pb1, pb2[,c("chr","midPos")], by=c("chr","midPos"))
pb3<-merge(pb3, pb2[,c("chr","midPos")], by=c("chr","midPos"))
pb1<-pb1[order(pb1$chr,pb1$midPos),]
pb2<-pb2[order(pb2$chr,pb2$midPos),]
pb3<-pb3[order(pb3$chr,pb3$midPos),]


# Are PBS values different depending on the combinations of 3 populations?
mean(pb3$PBS17-pb1$PBS17, na.rm=T)
#0.0003372794
mean(pb2$PBS17-pb1$PBS17, na.rm=T)
#-6.546866e-05
dif<-pb3$PBS07-pb1$PBS07
mean(dif) # -0.0003189494
mean(pb1$PBS91-pb2$PBS91)
# -7.008358e-05
mean(pb3$PBS96-pb2$PBS96)
#-0.0003748692
mean(pb3$PBS17-pb2$PBS17)
#0.0004027481
plot(dif,pch="." ) #different values 
mean(pb1$PBS17)

# differences in Fst values
dif2<-pb1$Fst91.17-pb2$Fst91.17
mean(dif2) # -4.500781e-06 
mean(pb1$Fst07.17-pb3$Fst07.17, na.rm=T)
#-1.80595e-05



#Merge Fst values
fst<-merge(pb1[,1:7], pb2[,c(1,5,7)], by="region")
fst<-merge(fst, pb3[,c(1,7)], by="region")
fst$ch<-as.integer(gsub("chr","",fst$chr))

fst<-fst[order(fst$ch, fst$midPos),]
fst$loc<-1:nrow(fst)

#odds<-paste0("chr",seq(1,26, by=2))
evens<-paste0("chr",seq(2,26, by=2))
fst$color<-"col1"
fst$color[fst$chr %in% evens]<-"col2"


ggplot(fst, aes(x=loc, y=Fst07.91, color=color))+
    geom_point(size=0.2)+
    scale_color_manual(values=c("gray30",grb))+
    theme_classic()+
    ylab("Fst")+xlab('')+
    theme(legend.position = "none")+
    ggtitle("PWS07 vs. PWS91")



####  PBS outliers comparison between files
pb1$pos.id<-paste0(pb1$chr,"_",pb1$midPos)
pb2$pos.id<-paste0(pb2$chr,"_",pb2$midPos)
pb3$pos.id<-paste0(pb3$chr,"_",pb3$midPos)

pb1.91.outlier<-pb1[order(pb1$PBS91, decreasing = T),]
pb1.91.outlier<-pb1.91.outlier[1:(nrow(pb1)/100),]

pb2.91.outlier<-pb2[order(pb2$PBS91, decreasing = T),]
pb2.91.outlier<-pb2.91.outlier[1:(nrow(pb2)/100),]

length(intersect(pb1.91.outlier$pos.id, pb2.91.outlier$pos.id))
#250 out of 723 are overlapping


pb1.07.outlier<-pb1[order(pb1$PBS07, decreasing = T),]
pb1.07.outlier<-pb1.07.outlier[1:(nrow(pb1)/100),]

pb3.07.outlier<-pb3[order(pb3$PBS07, decreasing = T),]
pb3.07.outlier<-pb3.07.outlier[1:(nrow(pb3)/100),]

length(intersect(pb1.07.outlier$pos.id, pb3.07.outlier$pos.id))
#258 out of 723 are overlapping

pb1.17.outlier<-pb1[order(pb1$PBS17, decreasing = T),]
pb1.17.outlier<-pb1.17.outlier[1:(nrow(pb1)/100),]

pb2.17.outlier<-pb2[order(pb2$PBS17, decreasing = T),]
pb2.17.outlier<-pb2.17.outlier[1:(nrow(pb2)/100),]

pb3.17.outlier<-pb3[order(pb3$PBS17, decreasing = T),]
pb3.17.outlier<-pb3.17.outlier[1:(nrow(pb3)/100),]

length(intersect(pb1.17.outlier$pos.id, pb3.17.outlier$pos.id))
#218 out of 723 are overlapping

length(intersect(pb2.17.outlier$pos.id, pb3.17.outlier$pos.id))
#2328 out of 723 are overlapping
length(Reduce(intersect, list(pb1.17.outlier$pos.id, pb2.17.outlier$pos.id, pb3.17.outlier$pos.id)))
# 137
Reduce(intersect, list(pb1.17.outlier$pos.id, pb2.17.outlier$pos.id, pb3.17.outlier$pos.id))
#  [1] "chr9_30235000"  "chr9_30245000"  "chr9_30225000"  "chr20_23485000" "chr9_30385000"  "chr9_30215000"  "chr9_30405000" 
#[8] "chr9_30395000"  "chr14_7135000"  "chr19_1995000"  "chr4_4685000"   "chr26_1615000"  "chr7_29495000"  "chr26_1635000" 
#[15] "chr19_5635000"  "chr12_28945000" "chr21_24415000" "chr19_5625000"  "chr16_5745000"  "chr10_19175000" "chr21_24585000"
#[22] "chr26_1625000"  "chr9_30255000"  "chr19_1985000"  "chr16_5735000"  "chr10_2095000"  "chr19_3045000"  "chr4_25275000" 
#[29] "chr19_5645000"  "chr12_28785000" "chr21_6125000"  "chr10_2105000"  "chr7_29505000"  "chr21_24595000" "chr10_19185000"
#[36] "chr6_1425000"   "chr19_2005000"  "chr23_18205000" "chr3_13935000"  "chr3_26445000"  "chr6_1415000"   "chr16_25815000"
#[43] "chr16_26825000" "chr17_16945000" "chr10_19195000" "chr11_7675000"  "chr2_7165000"   "chr10_2115000"  "chr26_11555000"
#[50] "chr2_7155000"   "chr13_25475000" "chr13_25485000" "chr5_5995000"   "chr17_19645000" "chr4_25355000"  "chr19_1555000" 
#[57] "chr14_2805000"  "chr20_24905000" "chr19_1545000"  "chr19_1535000"  "chr23_22165000" "chr2_7145000"   "chr5_5985000"  
#[64] "chr13_27735000" "chr19_7075000"  "chr19_7055000"  "chr21_26255000" "chr8_28495000"  "chr15_14765000" "chr2_3315000"  
#[71] "chr12_25205000" "chr10_3245000"  "chr5_6005000"   "chr8_28485000"  "chr20_20325000" "chr15_14755000" "chr19_1525000" 
#[78] "chr11_7665000"  "chr21_26245000" "chr19_7085000"  "chr2_7325000"   "chr12_25195000" "chr11_7685000"  "chr7_29735000" 
#[85] "chr2_7175000"   "chr3_13925000"  "chr2_4905000"   "chr13_25465000" "chr23_18215000" "chr11_7695000"  "chr5_13425000" 
#[92] "chr8_28475000"  "chr5_13415000"  "chr21_26265000" "chr18_12315000" "chr19_3075000"  "chr21_26225000" "chr5_3285000"  
#[99] "chr15_14775000" "chr10_21545000" "chr25_4995000"  "chr21_26235000" "chr5_13395000"  "chr26_11545000" "chr12_25185000"
#[106] "chr25_4985000"  "chr5_3265000"   "chr8_28465000"  "chr12_25165000" "chr12_28955000" "chr17_14465000" "chr3_21375000" 
#[113] "chr3_21405000"  "chr4_26265000"  "chr12_25175000" "chr3_26375000"  "chr2_3325000"   "chr25_5005000"  "chr20_24915000"
#[120] "chr17_14485000" "chr17_14475000" "chr19_35000"    "chr5_13405000"  "chr12_15385000" "chr2_7335000"   "chr23_1825000" 
#[127] "chr19_245000"   "chr12_8985000"  "chr19_3085000"  "chr25_2565000"  "chr5_3275000"   "chr15_25625000" "chr11_385000"  
#[134] "chr3_21365000"  "chr15_25635000" "chr8_25755000"  "chr19_3145000"






#Use the PBS96 analyzed with PBS 07
pbs<-merge(pb1, pb2[,c("region", "PBS17")], by="region")
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
ggsave("Output/fst_pbs/MD7000/PWS_pbs_3yrs_stacked.pdf", width = 10, height = 7)


#compare folded SFS vs. unfolded SFS results for PWS07 vs. PWS17

fold<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS17_folded_50kWindow")
unf<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS17_50kWindow")

mean(fold$Nsites) #[1] 0.009639717
mean(unf$Nsites) # [1] 0.008852109

fold<-read.delim("Data/new_vcf/angsd/Fst/folded_fst_PWS07_PWS91_50kWindow")
unf<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS91_50kWindow")

mean(fold$Nsites) # 0.008773152
mean(unf$Nsites)  # 0.008241249


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

#extract the regions and find the annotations
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
