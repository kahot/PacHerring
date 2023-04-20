
meta<-read.csv("Data/EVOS_MasterSheet_JoeMcGirr_April2020.csv")

table(meta$Year.Collected)
#1991 1996 2006 2007 2017 
#192  288  149   71  576 

table(meta$Location, meta$Year.Collected)
#                                      1991 1996 2006 2007 2017
#British Columbia                        0    0    0    0   96
#Cordova, PWS                            0    0    0    0   96
#Gravina Point                           0    0    0   71    0
#Green Island                           96    0    0    0    0
#Nunavachak                              0    0   75    0    0
#Port Orchard, Manchester, Washington    0    0    0    0   96
#Rocky Bay                               0   96    0    0    0
#San Francisco Bay, California           0    0    0    0   96
#Sitka Sound                             0   96   74    0   96
#Togiak Bay                             96   96    0    0   96

#Togiak Bay: Western end of the sampling region
#Sitka Sound: Eastern end of the sampling region


#Create master metadata sheet

pop_info <- read.table("Data/familiarize/EVOS_MasterSheet_JoeMcGirr_April2020.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
aligned <- read.table("Data/familiarize/aligned_samples.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
pop_info <- pop_info[pop_info$Sample %in% aligned$sample,]

samples <- read.table("Data/plink/plates_1_through_5_rm.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(samples)[names(samples)=="V1"] <- "Sample"
sample_info <- inner_join(samples,pop_info, by = "Sample")

sample_info$Pop<-sample_info$Population.Year
sample_info <- sample_info %>% separate(Pop, into = c("pop", "year"), sep = "(?<=[A-Za-z])(?=[0-9])")


write.csv(sample_info,"Data/Sample_metadata_892pops.csv", row.names = F)

sam$


sam<-read.csv("Data/Sample_metadata_892pops.csv")
sam$newID<-paste0(sam$Population.Year,"_", substr(sam$Sample, 1,4))
sink("Data/newID.txt")
cat(paste0(sam$newID,"\n"))
sink(NULL)

fam<-read.table("Data/vcfs/ph_maf05_id.fam")
fam$V2<-sam$newID
write.table(fam,"Data/vcfs/ph_id_modified.fam",  quote=F, sep=" ", col.names = F, row.names = F)



## Calculate the read count numbers for the samples used in the experiment 9.6.22

pops<-read.csv("Data/Sample_metadata_892pops.csv")
## Read picard metrics output files from Joe
wgs<-data.frame()
for (i in 1: nrow(pops)){
    id<-pops$Sample[i]
    df<-read.table(paste0("/Volumes/Kaho_Data/PacHerring/Data/familiarize/wgsMetrics/",id,".collect_wgs_metrics.txt"),header = TRUE, stringsAsFactors = FALSE, nrow= 1)
    df$Sample<-id
    wgs<-rbind(wgs, df)
}

popinfo<-pops[,c("Sample","Population.Year","pop","Year.Collected")]

bam_summary<-merge(popinfo, wgs, by="Sample")

byPopYr<-aggregate(bam_summary$MEAN_COVERAGE, by=list(bam_summary$Population.Year), mean)
write.csv(byPopYr, "Output/QC/Bam_metrics_summary.csv")

bam_summary$Population.Year<-factor(bam_summary$Population.Year, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "BC17","WA17","CA17"))
write.csv(bam_summary, "Output/QC/Read_coverage_bam_file_summary.csv")

library(ggplot2)
source("Rscripts/BaseScripts.R")

bam_summary<-read.csv("Output/QC/Read_coverage_bam_file_summary.csv", row.names = 1)
byPopYr<-read.csv("Output/QC/Bam_metrics_summary.csv", row.names = 1)

bam_summary$pop<-factor(bam_summary$pop, levels=c("TB","PWS","SS", "BC","WA","CA"))
byPopYr$Group.1<-factor(byPopYr$Group.1,  levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))

bam_summary$Population.Year<-factor(bam_summary$Population.Year, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))

ggplot()+
    geom_boxplot(data=bam_summary, aes(x=Population.Year, y=MEAN_COVERAGE, color=pop))+
    theme_classic()+ylab("Read coverage")+xlab('')+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    geom_point(data=byPopYr, aes(x=Group.1, y=x))+
    scale_color_manual(values=cols, guide="none")+
    annotate(geom="text", x = 12, y=2.45,label="Mean coverage = 1.061")
ggsave("Output/QC/Mean_coverage_plot_per_pop.year.png", width = 6, height = 4, dpi=300)
mean(bam_summary$MEAN_COVERAGE)
#1.060702 the dataset-wide coverage

byPopYr
#   Group.1         x
#1     BC17 1.1433385
#2     CA17 1.3067545
#3    PWS07 1.2245732
#4    PWS17 0.9596604
#5    PWS91 1.0283849
#6    PWS96 1.0858443
#7     SS06 0.9830740
#8     SS17 1.0059452
#9     SS96 1.1034367
#10    TB06 1.1069350
#11    TB17 0.9412129
#12    TB91 1.1359431
#13    TB96 0.9357896
#14    WA17 0.9048108

mean(byPopYr$x) #1.061836
mean(bam_summary$SD_COVERAGE)







reads<-popinfo
for (i in 1: nrow(pops)){
    id<-pops$Sample[i]
    df<-read.table(paste0("Data/QC/",id,".stats.txt"),header = F, fill =T )
    reads$Mapped_read_count[i]<-df[1,1]
    
}

### Check if there were any sequencing biases reflected in higher Fst observed between 1991 and 2006/07 samples

t1<-data.frame(table(pops$Sequence.Plate[pops$Year.Collected=="1991"]))
# 6  7  9 10 11 13 14 15 17 18 19 
# 3  8 15  7 22 16 16  8  8  5 24 
colnames(t1)<-c("Plate","1991")
t2<-data.frame(table(pops$Sequence.Plate[pops$Population.Year=="TB06"|pops$Population.Year=="PWS07"]))
# 6  7  8  9 11 12 14 15 17 18 20 
#15 23  8  8  8  9  3  8  8  8 
 
colnames(t2)<-c("Plate","2006")

plates<-merge(t1, t2, by="Plate", all=T)
plates[is.na(plates)]<-0
rownames(plates)<-plates$Plate

#chi square test
chisq.test(plates[,2:3])
#Pearson's Chi-squared test
#
#data:  plates[, 2:3]
#X-squared = 94.549, df = 12, p-value = 6.461e-15

table(pops$Sequence.Plate, pops$Year.Collected)

## TB06 and PWS07 have the higher pi out of all years, Not SS

table(pops$Sequence.Plate, pops$Population.Year)
#    BC17 CA17 PWS07 PWS17 PWS91 PWS96 SS06 SS17 SS96 TB06 TB17 TB91 TB96 WA17
#6    24    0    15     8     0     0    0    0    0    0    0    3    0    8
#7     8    0    15     0     0     0    8    0    8    8    8    8    0    0
#8     8    0     0     8     0     8   15    8    8    8    0    0    0    0
#9     8    8     8     0     7     0    0    8    0    0    0    8   16    0
#10    0    0     0     8     7     8    0    8    8    0    8    0    8    8
#11    0    8     0     8    15     0    0    0    0    8    0    7    0   16
#12    8    0     0     0     0     8    0    8    8    9    0    0    8    8
#13    0   14     0     0    16     0    0    8   15    0    0    0    0    8
#14    0    8     0     0     0     0    0    0   16    3    0   16    8    8
#15    0    8     0     8     0     8    8   16    0    8    0    8    0    0
#16    0    8     0     8     0    24    0    0    8    0    0    0   14    0
#17    0    0     0     8     0     8    8    0    0    8   16    8    3    0
#18    8    0     8     0     5     8    0    0    7    0   24    0    0    0
#19    0    0     0     0     8     0    0    8    0    0   16   16    8    8
#20    0   16     0     0     0     0    2    0    0    0    0    0    8    8

# No consistent plate contributing high pi (PWS07, CA17 PWS96 are the three highest pi)








### total reads
# use samtools to extract bam statistics (NoReads.sh)

files<-list.files("Data/QC/", pattern=".stats.txt")
reads<-data.frame(file=gsub(".stats.txt",'',files))
for (i in 1: length(files)){
    df<-read.table(paste0("Data/QC/", files[i]), fill = T)
    reads$Total[i]<-df[1,1]
    reads$Mapped[i]<-df[7,1]
}

write.csv(reads, "Output/QC/Bam_TotalReads_summary.csv")

reads$Total<-as.integer(reads$Total)
sum(reads$Total)
#9,065,600,620

