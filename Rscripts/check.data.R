
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




sam<-read.csv("Data/Sample_metadata_892pops.csv")
sam$newID<-paste0(sam$Population.Year,"_", substr(sam$Sample, 1,4))
sink("Data/newID.txt")
cat(paste0(sam$newID,"\n"))
sink(NULL)

fam<-read.table("Data/vcfs/ph_maf05_id.fam")
fam$V2<-sam$newID
write.table(fam,"Data/vcfs/ph_id_modified.fam",  quote=F, sep=" ", col.names = F, row.names = F)
