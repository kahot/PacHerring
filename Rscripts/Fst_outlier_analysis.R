# Fst Outlier analysis

library(OutFLANK)  # outflank package
library(vcfR)
library(adegenet)
library(dartR)
library(qvalue)
library(stringr)
source("Rscripts/data_subset.r")

vcf <- read.vcfR("Data/vcfs/PWS_SS/PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz")

#convert vcf to a genind object
vcf2 = vcfR2genind(vcf)

vcf2$pop = as.factor(substr(indNames(vcf2), 6, 10))

#subset the population 
popNames(vcf2)
pws1<-popsub(vcf2, sublist=c("PWS91", "PWS96"))
pws2<-popsub(vcf2, sublist=c("PWS91", "PWS07"))
pws3<-popsub(vcf2, sublist=c("PWS91", "PWS17"))
pws4<-popsub(vcf2, sublist=c("PWS96", "PWS17"))
pws5<-popsub(vcf2, sublist=c("PWS96", "PWS07"))
pws<-popsub(vcf2, sublist=c("PWS91","PWS96","PWS07","PWS17"))



# Run OutFLANK using dartR wrapper script
outflnk = gl.outflank(pws1, qthreshold = 0.05, plot = FALSE) #all are above q=0.01

outflnk = gl.outflank(pws5, qthreshold = 0.05, NumberOfSamples=2)
outflnk = gl.outflank(pws3, qthreshold = 0.05, NumberOfSamples=2)
outflnk = gl.outflank(pws, qthreshold = 0.05, NumberOfSamples=4)

#parameters to consider
#OutFLANK(fst,LeftTrimFraction=0.01,RightTrimFraction=0.01,
#               Hmin=0.05,NumberOfSamples=2,qthreshold=0.01)


# Extract OutFLANK results
outflnk.df = outflnk$outflank$results

# Remove duplicated rows for each SNP locus
rowsToRemove = seq(1, nrow(outflnk.df), by = 2)
outflnk.df = outflnk.df[-rowsToRemove, ]

# Print number of outliers (TRUE)
outflnk.df$OutlierFlag %>% summary
#   Mode   FALSE    TRUE    NA's 
#logical  202134       3    1008 

outlier_indexes = which(outflnk.df$OutlierFlag == TRUE)
outlierID = locNames(pws1)[outlier_indexes]
outlierID
#for pws91 vs. pws96 [1] "chr15_50932"   "chr16_2409773" "chr21_3520066"

#for pws91 vs. pws07 [1] "chr12_22754114"
#for pws91 vs. pws17 "   "chr1_756203"   "chr7_30659499"
#for pws96 vs. pws17 [1] "chr5_9799"      "chr11_7419995"  "chr11_24774233" "chr13_26955074" "chr16_2352301" 
# "chr17_22575376" "chr24_17373226"
#for pws96 vs. pws07 34 sties
#[1] "chr1_3123484"   "chr2_14324749"  "chr3_24222025"  "chr4_9083533"   "chr4_9084957"  
#[6] "chr4_22238158"  "chr5_8511831"   "chr5_23775393"  "chr5_26676625"  "chr6_8935875"  
#[11] "chr7_8102471"   "chr9_2071851"   "chr9_12429929"  "chr10_20397116" "chr10_20397477"
#[16] "chr13_11912894" "chr13_16816385" "chr13_16816423" "chr13_23781998" "chr16_802380"  
#[21] "chr16_1342977"  "chr17_7235506"  "chr17_17951088" "chr18_3453437"  "chr19_12448329"
#[26] "chr19_13586833" "chr20_20265999" "chr22_8035135"  "chr22_22029976" "chr22_22030337"
#[31] "chr23_6768314"  "chr23_15825559" "chr23_18895171" "chr24_2461074" 


# Convert Fsts <0 to zero
outflnk.df$FST[outflnk.df$FST < 0] = 0 

# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])

# Plot He versus Fst
ggplot(data = outflnk.df)+
    geom_point(aes(x = He, y = FST, colour = OutlierFlag), size=1)+
    scale_colour_manual(values = c("#80808066","red"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
    ggtitle("OutFLANK outlier test")+
    theme_bw()+
    xlab(hetlab)+
    ylab(fstlab)+
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("Output/Fst/PWS96vsPWS07.outflank_result.pdf", width = 6, height=4)

write.csv(outflnk.df,"Output/Fst/PWS96vsPWS07.outflank_result_q0.05.csv")



###########
# run comparison against Sitka pop

#subset the population 
#popNames(vcf2)
#Compare the same time

comp1<-popsub(vcf2, sublist=c("PWS96", "SS96"))
comp2<-popsub(vcf2, sublist=c("PWS07", "SS06"))
comp3<-popsub(vcf2, sublist=c("PWS17", "SS17"))

for (i in 2:3){
    comp<-get(paste0("comp",i))
    pops<- popNames(comp)
    print(pops)
    outflnk = gl.outflank(comp, qthreshold = 0.05, NumberOfSamples=2)
    outflnk.df = outflnk$outflank$results
    
    # Remove duplicated rows for each SNP locus
    rowsToRemove = seq(1, nrow(outflnk.df), by = 2)
    outflnk.df = outflnk.df[-rowsToRemove, ]
    outlier_indexes = which(outflnk.df$OutlierFlag == TRUE)
    outlierID = locNames(pws1)[outlier_indexes]
    # Print number of outliers (TRUE)
    print(outflnk.df$OutlierFlag %>% summary)
    print(outlierID)
    
    outflnk.df$FST[outflnk.df$FST < 0] = 0 
    
    # Italic labels
    fstlab = expression(italic("F")[ST])
    hetlab = expression(italic("H")[e])
    
    ggplot(data = outflnk.df)+
        geom_point(aes(x = He, y = FST, colour = OutlierFlag), size=1)+
        scale_colour_manual(values = c("#80808066","red"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
        ggtitle("OutFLANK outlier test")+
        theme_bw()+
        xlab(hetlab)+
        ylab(fstlab)+
        theme(legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
        ggtitle(paste0(pops[1]," vs.",pops[2]))
    ggsave(paste0("Output/Fst/",pops[1],"vs",pops[2],".outflank_result.pdf"), width = 6, height=4)
    write.csv(outflnk.df,paste0("Output/Fst/",pops[1],"vs",pops[2],".outflank_result_q0.05.csv"))
}

#Compare the 3 time points
Outlier<-data.frame()
for (i in 1:3){
    comp<-get(paste0("comp",i))
    pops<- popNames(comp)
    print(pops)
    df<-read.csv(paste0("Output/Fst/",pops[1],"vs",pops[2],".outflank_result_q0.05.csv"), row.names = 1)
    df<-df[df$OutlierFlag==T,]
    df<-df[!is.na(df$FST),]
    if (nrow(df)==0) next
    else {
        df$Comparison<-paste0(pops[1],"vs",pops[2])
        Outlier<-rbind(Outlier, df)
    }
}

Outlier$chr<-substr(Outlier$LocusName, 1, 5)
Outlier$chr<-gsub("_", "", Outlier$chr)

Outlier$pos<-str_extract(Outlier$LocusName, ("_"),("\\d.+"))

