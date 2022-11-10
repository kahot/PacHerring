#from https://github.com/pinskylab/codEvol

## examine results from site reshuffled FST from ANGSD genotypes
## run after angsd_fst_siteshuffle_null.sh

#################
# parameters
#################
minloci <- 10 # should match angsd_fst_siteshuffle_null.r
winsz <- 50000 # window size
winstp <- 10000 # window step

###########################
# load functions
###########################
require(data.table)
#require(plyr)
require(ggplot2)
require(RColorBrewer)

calcp <- function(fst, null) return((sum(null > fst)+1)/(length(null)+1)) # equation from North et al. 2002 Am J Hum Gen

#####################
# read in and prep data
#####################

years<-c("PWS91","PWS96","PWS07","PWS17")
comb<-t(combn(years, 2))

# continuous nucleotide position for the whole genome
chrmax <- fread('Data/new_vcf/chr_sizes.bed')
chrmax<-chrmax[,-2]
colnames(chrmax)<-c("chr", "len")
chrmax$start<-c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

#setkey(chrmax, chr)

prunedFst_gw<-data.frame(pop1=comb[,1],pop2=comb[,2])
#for (p in 1: nrow(comb)){
    for (p in 1:5){
    #genome-wide max Fst from reshuffling (unlinked sites)
    null<-fread(paste0('Data/shuffle/', comb[p,1],".",comb[p,2],'_fst_siteshuuffle.csv.gz'))
    #window-based Fst from angsd            
    fstwin <- fread(paste0('Data/new_vcf/angsd/fromVCF/2D/fst_',comb[p,1],"_",comb[p,2],'_50kWindow_maf00'), header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) # all sites
    #persite fst (AB) -pruned sites
    AB <- fread(paste0('Data/shuffle/fst_',comb[p,1],"_",comb[p,2],'_persite_maf00_ldtrim.csv.gz'))
    
    # merge nucleotide position into the frequency files
    setkey(AB, CHROM)
    AB <- AB[chrmax[, .(CHROM = chr, start)], ]
    AB[, posgen := POS + start]
    AB[,start := NULL]
    
    # Calc genome-wide FST
    prunedFst_gw$gwFst[p]<-AB[!is.na(A), sum(A)/sum(B)]
    
    # Calc windowed FST
    # create new columns as indices for windows
    for(j in 1:(winsz/winstp)){
        AB[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    }
    
    # calc fst and # snps per window
    for(j in 1:(winsz/winstp)){
        if(j ==1){
            fstwin <- AB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
        } 
        if(j > 1){
            fstwin <- rbind(fstwin, AB[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
        } 
    }
    
    ## null model stats
    fstats<-null[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
    prunedFst_gw$null.Fstmax[p]<-fstats$max[1]
    prunedFst_gw$null.Fst.95q[p]<-fstats$u95[1]
    
    #Calculate p-values
    pval<-fstwin[, p := calcp(fst, null$x), by = .(CHROM, midPos)]
    
    
    #combined the datasets
    pair<-paste0(comb[p,1],"_",comb[p,2])
    fstwin[, pop := pair ]
    write.csv(fstwin, paste0("Output/Fst/fst_siteshuffle.pairwise_", pair,".csv"))
    #data <- rbind(data, fstwin)
}

write.csv(prunedFst_gw, "Output/Fst/fst_siteshuffle_pws_summary.csv")

data<-data.frame()
for(i in 1: nrow(comb)){
    pair<-paste0(comb[i,1],"_",comb[i,2])
    df<-fread(paste0("Output/Fst/fst_siteshuffle.pairwise_", pair,".csv"))
    df<-df[!is.na(fst) & midPos > 0, ]
    df[,ch:= as.integer(gsub("chr",'', CHROM))]
    setkey(df,CHROM)
    df<-df[chrmax[, .(CHROM = chr, start)], ]
    df[,pos.gw:= midPos+start]
    #setkey(df, ch, midPos)
    data<-rbind(data, df)
}

write.csv(data, file = gzfile('Output/Fst/fst_siteshuffle.angsd.PWS.csv.gz'), row.names = FALSE)
    

# plots
#pop as factor 
data[, pop := factor(pop, levels = c("PWS91_PWS96","PWS91_PWS07","PWS91_PWS17","PWS96_PWS07","PWS96_PWS17","PWS07_PWS17"))]

# plot p-value vs. position
cols <- brewer.pal(4, 'Paired')[rep(1:2,13)]
#isplay.brewer.all(n=4)
cols<-rep(c("steelblue","lightblue"), times=13)

ggplot(data[nloci >= minloci, ], aes(pos.gw, -log10(p), color = factor(ch))) + 
    geom_point(size = 0.2, alpha = 0.3) +
    facet_wrap(~pop, ncol = 1) +
    scale_color_manual(values = cols, guide='none') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')+
    theme_bw()+
    theme(axis.text.x = element_blank())+
    xlab("Genome")
ggsave("Output/Fst/fst.siteshuffle.p_vs_pos.PWS.png", width = 7.5, height = 8.5,dpi = 300)


# plot p-value vs. nloci (gatk loci)
ggplot(data, aes(nloci, -log10(p), color = pop)) + 
    geom_point(size = 0.2, alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') + 
    scale_x_log10()+theme_bw()
ggsave("Output/Fst/fst.siteshuffle.p_vs_nloci.PWS.png", width = 6, height =4,dpi = 300)



#################
# print outliers
#################


outliers<-data[p < 0.05, ]
write.csv(outliers, "Output/Fst/Fst_nullshuffling_outliers_PWS.csv")
