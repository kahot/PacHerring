#Run detectRUNS to identify ROH regions
library(detectRUNS)
library(ggplot2)
library(gridExtra)
#source("Rscripts/detectRuns_customizedPlots.R")

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected","Population.Year")]
colnames(pop_info)[3]<-"year"

#ph data
#plink ped is tab-delimited format. Need to convert it to space-delimited format;
#run   "  sed 's/\t/ /g' in.ped > out.ped   "

pops<-unique(pop_info$Population.Year)
i=4
pops[i]
gtfilepath<-paste0("Output/ROH/ped/", pops[i],".ped")
mappath<-paste0("Output/ROH/ped/", pops[i],"_maf00.map")


#maxMissWindow needs to be >=3 to have reasonable numbers for our data
slidingRuns <- slidingRUNS.run(
    genotypeFile = gtfilepath,
    mapFile = mappath, 
    windowSize = 100, 
    threshold = 0.05,
    minSNP = 50, 
    ROHet = FALSE, 
    maxOppWindow = 1, 
    maxMissWindow = 3,
    maxGap = 500000, 
    #minLengthBps = 1000,  #minimum length of run in bps 
    minDensity = 1/100000, # SNP/kbps
    maxOppRun = NULL,
    maxMissRun = NULL
) 


consecutiveRuns<-consecutiveRUNS.run(
    genotypeFile =gtfilepath,
    mapFile = mappath,
    minSNP = 20,
    ROHet = F,
    maxGap = 10^6,
    minLengthBps = 100000,
    maxOppRun = 1,
    maxMissRun = 2
)

#Export the results as csv files
write.csv(slidingRuns, "Output/ROH/PWS17_slidingruns_100-3-500k_100kdensity.csv")
#write.csv(consecutiveRuns, "Output/ROH/consecutiveruns_15-2-100k.csv")

#slidingRuns<-read.csv("Output/ROH/PWS17_slidingruns_50-3-500k.csv")


#Create summary
summaryList <- summaryRuns(
    runs = slidingRuns, mapFile = mappath, genotypeFile = gtfilepath, 
    Class = 1, snpInRuns = T)

summaryList2 <- summaryRuns(
    runs = consecutiveRuns, mapFile = mappath, genotypeFile = gtfilepath, 
    Class = 1, snpInRuns = F)


#ROH count and size
count1<-summaryList$summary_ROH_count
#write.csv(count.t,"Output/ROH/ROHcount_per_sample_sliding_15-3-100k.csv")

#count2<-summaryList2$summary_ROH_count
#count.t2<-data.frame(t(count2))
#colnames(count.t2)[1]<-"size0-6"
#write.csv(count.t2,"Output/ROH/ROHcount_per_sample2_consecutive_15_2_100k.csv")

count<-cbind(count.t,count.t2)
colnames(count)<-c("sliding","consecutive")
count$pop<-rownames(count)
countm<-melt(count, id.vars="pop")
ggplot(countm,aes(x=pop, y=value, fill=variable))+
    geom_bar(stat="identity",position=position_dodge(width = 0.8))+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90), legend.title = element_blank())+
    ylab("Total N of ROH ")
ggsave("Output/ROH/Total_number_of_ROH_perPopGroup.pdf", width = 7, height = 4.5)

#the average number of ROH per chromosome 
count_chr<-summaryList$summary_ROH_mean_chr
write.csv(count_chr,"Output/ROH/ROHcount_per_chr_sliding_50-3-500k.csv")

#%ROH per chromosome
pchr<-summaryList$summary_ROH_percentage_chr
#write.csv(pchr,"Output/ROH/ROHpercent_per_chr_per_sample_sliding_15-3-100k.csv")

pops<-unique(pop_info$Population.Year)

#pchr$pop<-rownames(pchr)
#pchr[,1:26]<-apply(pchr[,1:26], 1, as.numeric)

#colnames(pchr)[1:26]<-paste0("chr",colnames(pchr)[1:26])
#pchrm<-melt(pchr, id.vars="pop")

ggplot(pchr,aes(x=chrom, y=PWS17))+
    geom_bar(stat="identity")



pchrm$pop<-factor(pchrm$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","TB91","TB96","TB06","TB17","WA17","CA17","BC17"))
plots<-list()
for (i in 1:26){
    ch<-paste0("chr",i)
    df<-pchrm[pchrm$variable==ch,]
    plots[[i]]<-ggplot(df,aes(x=pop, y=value, fill=pop))+
        geom_bar(stat="identity")+
        theme_classic()+
        theme(legend.position = "none")+
        theme(axis.text.x = element_text(angle=90))+
        ylab("% ROH per chromosome")+xlab('')+
        ggtitle(ch)
}
#pdf("Output/ROH/PercentROH_perChromosome_perPop.pdf", width = 18, height=26)
#do.call(grid.arrange, c(plots, ncol=3))
#dev.off()



#####

pchr.sum<-pchr[,1:2]
pchr.sum$Mean<-rowMeans(pchr[,3:28])

ggplot(pchr.sum,aes(x=Population.Year, y=Mean))+
    geom_point(color="blue", position=position_jitter(width = .1))+
    theme_bw()
ggsave("Output/ROH/AveROHpercent_perPopulation.pdf", width = 8, height = 5)

#percentage of ROH class in each individual
per<-summaryList$summary_ROH_percentage

mean(pchr$PWS17)
### Froh
# This function calculates the individual inbreeding coefficients based on runs of
#' homozygosity (ROH), either per-chromosome (chromosome-wide) or based on the
#' entire genome (genome-wide). See details of calculations below

#' This function calculates the individual inbreeding coefficients based on runs of
#' homozygosity (ROH) using only ROH of specific size classes.
df1<-summaryList$result_Froh_class

pops<-unique(pop_info$Population.Year)
fclass<-data.frame()
for (i in 1: length(pops)){
    df<-df1[df1$id==pops[i],]
    fclass<-rbind(fclass,df[1,])
}

write.csv(fclass, "Output/ROH/Froh_class2.csv")

df2<-summaryList$result_Froh_genome_wide

#pop
gwide<-data.frame(pop=pops)
for (i in 1: length(pops)){
    df<-df2[df2$group==pops[i],]
    gwide$sum[i]<-mean(df$sum, na.rm = T)
    gwide$Froh_genome[i]<-mean(df$Froh_genome, na.rm=T)
}

write.csv(gwide, "Output/ROH/Froh_genome_wide_sliding_15-3-100k.csv")

df3<-summaryList$result_Froh_chromosome_wide
cwide<-data.frame(pop=pops)
for (i in 1: length(pops)){
    df<-df3[df3$group==pops[i],]
    cwide[i,2:27]<-colMeans(df[,3:28], na.rm = T)
}

write.csv(cwide, "Output/ROH/Froh_perChromosome_sliding_15-3-100k.csv")


foh<-df2[,c(2,4)]
ggplot(df2,aes(x=group, y=Froh_genome))+
    geom_point(color="pink",position=position_jitter(width = .1))+
    theme_bw()+xlab('')

pchrm$pop<-factor(pchrm$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","TB91","TB96","TB06","TB17","WA17","CA17","BC17"))
plots<-list()
for (i in 1:26){
    ch<-paste0("chr",i)
    df<-pchrm[pchrm$variable==ch,]
    plots[[i]]<-ggplot(df,aes(x=pop, y=value, fill=pop))+
        geom_bar(stat="identity")+
        theme_classic()+
        theme(legend.position = "none")+
        theme(axis.text.x = element_text(angle=90))+
        ylab("% ROH per chromosome")+xlab('')+
        ggtitle(ch)
}
pdf("Output/ROH/PercentROH_perChromosome_perPop.pdf", width = 18, height=26)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()





#The dataframe “SNPinRun” contains, for each SNP, the proportion of times it falls inside a run in any given population/group
head(summaryList$SNPinRun)


### PLOTS ###

# Plot directly all runs detected in each individual against their position along the chromosome
# Separate plots per chromosome are produced

#plot_Runs <- function(runs, suppressInds=FALSE, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {

plot_Runs(runs = slidingRuns, suppressInds =F, savePlots=T, separatePlots=F, outputName="Output/ROH/PWS17_Slidingrun_plot_")

#runs can be plotted against their position along the chromosome
plot_StackedRuns(runs=slidingRuns, savePlots=T, separatePlots=FALSE, outputName="Output/ROH/PWS17_SlidingRuns_stacked_100-3-500k")

plot_StackedRuns(runs=slidingRuns[slidingRuns$group=="PWS91",])


#plot_StackedRuns_ph(runs=consecutiveRuns, savePlots=TRUE, separatePlots=FALSE, outputName=NULL)


#the proportion of times each SNP falls inside a run in any given population/group
plot_SnpsInRuns(
    runs = slidingRuns[slidingRuns$chrom==1,], genotypeFile = gtfilepath, 
    mapFile = mappath)


        
            
            
            
# identify the position of a runs (ROH) peak
# Thershold =  % of individuals in that population/groupo 
topRuns <- tableRuns(
    runs =  slidingRuns,  genotypeFile = gtfilepath, mapFile = mappath, 
    threshold = 0.7)




# the proportion of times each SNP falls inside a run plotted against SNP positions in all chromosomes together
plot_manhattanRuns(
    runs = slidingRuns[slidingRuns$group=="PWS96",], 
    genotypeFile = gtfilepath, mapFile = mappath)
