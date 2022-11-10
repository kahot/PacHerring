# shuffle ANGSD FST A and B values across sites and calculate windowed FST to get a null distribution of max genome-wide FST
# run last part of angsd_fst.sh first to get the *.fst.AB.gz files
# if using GATK nodam2 loci, groups into linkage blocks based on ngsLD output. Need to run ngsLD_find_blocks.sh/.r first

# to run on saga

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) != 1) stop("Have to specify whether all loci (0) or GATK loci (1)", call.=FALSE)

gatkflag <- as.numeric(args[1])

if(gatkflag == 1){
    print('Using only GATK nodam2 loci. Will trim to unlinked blocks.')
} else if(gatkflag == 0){
    print('Using all loci. Note that this option does not trim to unlinked groups yet.')
} else {
    stop(paste(gatkflag, ' is not 0 or 1. Please only specify 0 (all loci) or 1 (gatk loci).'))
}


# parameters
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles
minloci <- 2 # minimum number of loci per window to consider


outfile1<-'Output/Analysis/PWS07_PWS17_siteshuuffle.csv.gz'
outfilecan <- 'analysis/Can_40.Can_14.fst.siteshuffle.csv.gz' # used if all loci are used
outfilelof0711 <- 'analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz'
outfilelof0714 <- 'analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz'
outfilelof1114 <- 'analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz'

# load functions
require(data.table)


#############
# Prep data
#############
# load fst A/B data
#can <- fread('analysis/Can_40.Can_14.fst.AB.gz')
#setnames(can, c('CHROM', 'POS', 'A', 'B'))

pws0717 <- fread('Data/new_vcf/angsd/fromVCF/2D/fst_PWS07_PWS17_persite_maf00.txt.gz')
setnames(pws0717, c('CHROM', 'POS', 'A', 'B'))

#lof0711 <- fread('analysis/Lof_07.Lof_11.fst.AB.gz')
#setnames(lof0711, c('CHROM', 'POS', 'A', 'B'))
#
#lof0714 <- fread('analysis/Lof_07.Lof_14.fst.AB.gz')
#setnames(lof0714, c('CHROM', 'POS', 'A', 'B'))
#
#lof1114 <- fread('analysis/Lof_11.Lof_14.fst.AB.gz')
#setnames(lof1114, c('CHROM', 'POS', 'A', 'B'))


# trim to gatk loci if requested
#if(gatkflag == 1){
    # list of loci to use
    gatk <- fread('data_2020.05.07/GATK_filtered_SNP_no_dam2.tab')
    setnames(gatk, c('CHROM', 'POS', 'REF', 'ALT'))
    
    # trim
    can <- can[gatk, on = c('CHROM', 'POS')]
    lof0711 <- lof0711[gatk, on = c('CHROM', 'POS')]
    lof0714 <- lof0714[gatk, on = c('CHROM', 'POS')]
    lof1114 <- lof1114[gatk, on = c('CHROM', 'POS')]
    
    # set new outfile names
    outfilecan <- 'analysis/Can_40.Can_14.gatk.nodam.fst.siteshuffle.csv.gz'
    outfilelof0711 <- 'analysis/Lof_07.Lof_11.gatk.nodam.fst.siteshuffle.csv.gz'
    outfilelof0714 <- 'analysis/Lof_07.Lof_14.gatk.nodam.fst.siteshuffle.csv.gz'
    outfilelof1114 <- 'analysis/Lof_11.Lof_14.gatk.nodam.fst.siteshuffle.csv.gz'
}

???
# remove unplaced
#can <- can[!(CHROM %in% 'Unplaced'), ]
#lof0711 <- lof0711[!(CHROM %in% 'Unplaced'), ]
#lof0714 <- lof0714[!(CHROM %in% 'Unplaced'), ]
#lof1114 <- lof1114[!(CHROM %in% 'Unplaced'), ]


# trim to unlinked loci if GATK loci were requested
if(gatkflag == 1){
    ld <- fread('analysis/ld.blocks.gatk.nodam.csv.gz') # linkage blocks from ngsLD_find_blocks.r
    
    can <- merge(can, ld[, .(CHROM, POS, cluster = cluster_can)], all.x = TRUE)
    lof0711 <- merge(lof0711, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)
    lof0714 <- merge(lof0714, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)
    lof1114 <- merge(lof1114, ld[, .(CHROM, POS, cluster = cluster_lof)], all.x = TRUE)
    
    # average Fst in linkage blocks and find a locus near the middle of each block to keep
    findmid <- function(POS){ # function to return a value near the middle of a vector of positions
        mn <- mean(range(POS))
        return(POS[which.min(abs(POS - mn))])
    }
    
    can[, keep := 1] # column for marking which loci to keep
    can[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
    canclust <- can[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
    can <- merge(can, canclust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
    can[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
    can <- can[keep == 1, ] # drop all loci in clusters that aren't midpoints
    can[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns
    
    lof0711[, keep := 1] # column for marking which loci to keep
    lof0711[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
    lof0711clust <- lof0711[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
    lof0711 <- merge(lof0711, lof0711clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
    lof0711[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
    lof0711 <- lof0711[keep == 1, ] # drop all loci in clusters that aren't midpoints
    lof0711[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns
    
    lof0714[, keep := 1] # column for marking which loci to keep
    lof0714[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
    lof0714clust <- lof0714[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
    lof0714 <- merge(lof0714, lof0714clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
    lof0714[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
    lof0714 <- lof0714[keep == 1, ] # drop all loci in clusters that aren't midpoints
    lof0714[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns
    
    lof1114[, keep := 1] # column for marking which loci to keep
    lof1114[!is.na(cluster), keep := 0] # in general, drop loci that are in a cluster
    lof1114clust <- lof1114[!is.na(cluster), .(CHROM = unique(CHROM), midPOS = findmid(POS), nloci = length(POS), sumA = sum(A), sumB = sum(B)), by = cluster] # find the cluster midpoints and sum of the FST components
    lof1114 <- merge(lof1114, lof1114clust[, .(CHROM, POS = midPOS, sumA, sumB, nloci, clustermid = 1)], all.x = TRUE) # label the cluster midpoints
    lof1114[clustermid == 1, ':='(keep = 1, A = sumA, B = sumB)] # move A and B from cluster mids over to the right column for FST calculations and mark them to keep
    lof1114 <- lof1114[keep == 1, ] # drop all loci in clusters that aren't midpoints
    lof1114[, ':='(keep = NULL, sumA = NULL, sumB = NULL, clustermid = NULL)] # drop extra columns
    
    # write out fst trimmed to unlinked blocks
    write.csv(can, gzfile('analysis/Can_40.Can_14.gatk.nodam.ldtrim.fst.AB.csv.gz'))
    write.csv(lof0711, gzfile('analysis/Lof_07.Lof_11.gatk.nodam.ldtrim.fst.AB.csv.gz'))
    write.csv(lof0714, gzfile('analysis/Lof_07.Lof_14.gatk.nodam.ldtrim.fst.AB.csv.gz'))
    write.csv(lof1114, gzfile('analysis/Lof_11.Lof_14.gatk.nodam.ldtrim.fst.AB.csv.gz'))
    
}

pw<-pws0717
# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
    pw[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    #lof0711[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    #lof0714[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    #lof1114[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# mark windows with < minloci for removal
rem <- 0 # number of windows removed for each of the 4 comparisons
for(j in 1:(winsz/winstp)){
    pwwin <- pw[, .(nsnps = length(POS)), by = .(win = get(paste0('win', j)))] # calc num snps per window
    rem[1] <- rem[1] + pwwin[, sum(nsnps < minloci)] # record number to be removed
    pwwin[, (paste0('win', j, 'keep')) := 1] # create col to mark which windows to keep
    pwwin[nsnps < minloci, (paste0('win', j, 'keep')) := 0] # mark windows to remove
    pwwin[, nsnps := NULL] # drop column
    setnames(pwwin, "win", paste0('win', j)) # change col name
    pw <- merge(pw, pwwin, by = paste0('win', j), all.x = TRUE) # merge keeper col back to full dataset
    
}

rem # number of windows removed for each comparison




####################################
# shuffle and recalc windowed FST
####################################
colnms <- c('CHROM', 'POS', paste0('win', 1:(winsz/winstp)), paste0('win', 1:(winsz/winstp), 'keep')) # list of column names we want out of the base data.table

# PWS07-17
print('Starting PWS07-17')

for(i in 1:nrep){
    cat(i); cat(' ')
    # create new dataset
    inds <- sample(1:nrow(pw), nrow(pw), replace = FALSE)
    temp <- cbind(pw[, ..colnms], pw[inds, .(A, B)]) # shuffle FSTs across positions
    
    # calc fst for each window to keep
    for(j in 1:(winsz/winstp)){
        temp2 <- temp[get(paste0('win', j, 'keep')) == 1, ] # trim to windows to keep. can't combine with next line for some reason.
        if(j ==1) tempfsts <- temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))]
        if(j > 1) tempfsts <- rbind(tempfsts, temp2[, .(fst = sum(A)/sum(B)), by = .(CHROM, POS = get(paste0('win', j)))])
    }
    
    # save the max windowed fst
    # exclude windows with negative midpoints
    if(i == 1) maxfst <- tempfsts[POS > 0, max(fst, na.rm = TRUE)]	
    if(i > 1) maxfst <- c(maxfst, tempfsts[POS > 0, max(fst, na.rm = TRUE)])
}

print(paste('Max:', max(maxfst, na.rm = TRUE), '; 95th:', quantile(maxfst, prob = 0.95, na.rm = TRUE)))

write.csv(maxfst, gzfile(outfile1), row.names = FALSE)

rm(maxfst)








## examine results from site reshuffled FST from ANGSD genotypes
## run after angsd_fst_siteshuffle_null.sh

#################
# parameters
#################
minloci <- 2 # should match angsd_fst_siteshuffle_null.r
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

# max FST per genome from reshuffling (all sites)
# nullcan <- fread('analysis/Can_40.Can_14.fst.siteshuffle.csv.gz')
# nulllof0711 <- fread('analysis/Lof_07.Lof_11.fst.siteshuffle.csv.gz')
# nulllof0714 <- fread('analysis/Lof_07.Lof_14.fst.siteshuffle.csv.gz')
# nulllof1114 <- fread('analysis/Lof_11.Lof_14.fst.siteshuffle.csv.gz')

# max FST per genome from reshuffling (GATK nodam2 sites, combined in linkage blocks, >1 SNP per window)
nullcangatk <- fread('analysis/Can_40.Can_14.gatk.nodam.fst.siteshuffle.csv.gz')
nulllof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.nodam.fst.siteshuffle.csv.gz')
nulllof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.nodam.fst.siteshuffle.csv.gz')
nulllof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.nodam.fst.siteshuffle.csv.gz')


# sliding window FST A and B components from ANGSD, after collapsing to unlinked loci
# header is missing the fst column, so have to skip and make our own
# need to make the all loci AB files
# can <- fread('analysis/Can_40.Can_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) # all sites
# lof0711 <- fread('analysis/Lof_07.Lof_11.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
# lof0714 <- fread('analysis/Lof_07.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 
# lof1114 <- fread('analysis/Lof_11.Lof_14.slide', skip = 1, header = FALSE, col.names = c('region', 'chr', 'midPos', 'Nsites', 'fst')) 

cangatk <- fread('analysis/Can_40.Can_14.gatk.nodam.ldtrim.fst.AB.csv.gz', drop = 1) # gatk sites
lof0711gatk <- fread('analysis/Lof_07.Lof_11.gatk.nodam.ldtrim.fst.AB.csv.gz', drop = 1) 
lof0714gatk <- fread('analysis/Lof_07.Lof_14.gatk.nodam.ldtrim.fst.AB.csv.gz', drop = 1) 
lof1114gatk <- fread('analysis/Lof_11.Lof_14.gatk.nodam.ldtrim.fst.AB.csv.gz', drop = 1) 


# nucleotide position for the whole genome (start position for each chr)
chrmax <- fread('data/lg_length.csv')
chrmax$start=c(0,cumsum(chrmax$len)[1:(nrow(chrmax)-1)])
setkey(chrmax, chr)

# merge nucleotide position into the frequency files
# setkey(can, chr)
# can <- can[chrmax[, .(chr, start)], ]
# can[, posgen := midPos + start]
# can[,start := NULL]
# 
# setkey(lof0711, chr)
# lof0711 <- lof0711[chrmax[, .(chr, start)], ]
# lof0711[, posgen := midPos + start]
# lof0711[,start := NULL]
# 
# setkey(lof0714, chr)
# lof0714 <- lof0714[chrmax[, .(chr, start)], ]
# lof0714[, posgen := midPos + start]
# lof0714[,start := NULL]
# 
# setkey(lof1114, chr)
# lof1114 <- lof1114[chrmax[, .(chr, start)], ]
# lof1114[, posgen := midPos + start]
# lof1114[,start := NULL]

setkey(cangatk, CHROM)
cangatk <- cangatk[chrmax[, .(CHROM = chr, start)], ]
cangatk[, posgen := POS + start]
cangatk[,start := NULL]

setkey(lof0711gatk, CHROM)
lof0711gatk <- lof0711gatk[chrmax[, .(CHROM = chr, start)], ]
lof0711gatk[, posgen := POS + start]
lof0711gatk[,start := NULL]

setkey(lof0714gatk, CHROM)
lof0714gatk <- lof0714gatk[chrmax[, .(CHROM = chr, start)], ]
lof0714gatk[, posgen := POS + start]
lof0714gatk[,start := NULL]

setkey(lof1114gatk, CHROM)
lof1114gatk <- lof1114gatk[chrmax[, .(CHROM = chr, start)], ]
lof1114gatk[, posgen := POS + start]
lof1114gatk[,start := NULL]

######################
# Calc genome-wide FST
######################

# already removes unplaced
cangatk[!is.na(A), sum(A)/sum(B)]
lof0711gatk[!is.na(A), sum(A)/sum(B)]
lof0714gatk[!is.na(A), sum(A)/sum(B)]
lof1114gatk[!is.na(A), sum(A)/sum(B)]


######################
# Calc windowed FST
######################

# create new columns as indices for windows
for(j in 1:(winsz/winstp)){
    cangatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    lof0711gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    lof0714gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    lof1114gatk[, (paste0('win', j)) := floor((POS - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
}

# calc fst and # snps per window
for(j in 1:(winsz/winstp)){
    if(j ==1){
        canfstsgatk <- cangatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
        lof0711fstsgatk <- lof0711gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
        lof0714fstsgatk <- lof0714gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
        lof1114fstsgatk <- lof1114gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))]
    } 
    if(j > 1){
        canfstsgatk <- rbind(canfstsgatk, cangatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
        lof0711fstsgatk <- rbind(lof0711fstsgatk, lof0711gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
        lof0714fstsgatk <- rbind(lof0714fstsgatk, lof0714gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
        lof1114fstsgatk <- rbind(lof1114fstsgatk, lof1114gatk[, .(fst = sum(A)/sum(B), nloci = length(POS)), by = .(CHROM, midPos = get(paste0('win', j)))])
    } 
}


#######################
## null model stats
#######################

# nullcan[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof0711[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof0714[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
# nulllof1114[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

nullcangatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0711gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof0714gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]
nulllof1114gatk[, .(max = max(x), u95 = quantile(x, probs = 0.95))]

###########################
# calc p-values per window
###########################

# can[, p := calcp(fst, nullcan$x), by = .(chr, midPos)]
# lof0711[, p := calcp(fst, nulllof0711$x), by = .(chr, midPos)]
# lof0714[, p := calcp(fst, nulllof0714$x), by = .(chr, midPos)]
# lof1114[, p := calcp(fst, nulllof1114$x), by = .(chr, midPos)]

canfstsgatk[, p := calcp(fst, nullcangatk$x), by = .(CHROM, midPos)]
lof0711fstsgatk[, p := calcp(fst, nulllof0711gatk$x), by = .(CHROM, midPos)]
lof0714fstsgatk[, p := calcp(fst, nulllof0714gatk$x), by = .(CHROM, midPos)]
lof1114fstsgatk[, p := calcp(fst, nulllof1114gatk$x), by = .(CHROM, midPos)]


######################
# combine the datasets
######################
# all loci
# lof0711[, pop := 'lof0711']
# lof0714[, pop := 'lof0714']
# lof1114[, pop := 'lof1114']
# can[, pop := 'can']
# 
# dat <- rbind(lof0711, lof0714, lof1114, can)
# nrow(dat)
# dat
# 
# dat[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]


# gatk loci
canfstsgatk[, pop := 'can']
lof0711fstsgatk[, pop := 'lof0711']
lof0714fstsgatk[, pop := 'lof0714']
lof1114fstsgatk[, pop := 'lof1114']

datgatk <- rbind(canfstsgatk, lof0711fstsgatk, lof0714fstsgatk, lof1114fstsgatk)
nrow(datgatk)
datgatk

datgatk[, pop := factor(pop, levels = c('can', 'lof0711', 'lof0714', 'lof1114'))]

# remove NAs and negative windows
datgatk <- datgatk[!is.na(fst) & midPos > 0, ]

# sort by window
setkey(datgatk, pop, CHROM, midPos)

##############
# Write out
##############

# write.csv(dat, file = gzfile('output/fst_siteshuffle.angsd.csv.gz'), row.names = FALSE)
write.csv(datgatk, file = gzfile('output/fst_siteshuffle.angsd.gatk.csv.gz'), row.names = FALSE)

##############
# plots
##############

# plot p-value vs. position (all loci)
# cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
# p1 <- ggplot(dat, aes(posgen, -log10(p), color = chr)) + 
#   geom_point(size = 0.2, alpha = 0.3) +
#   facet_wrap(~pop, ncol = 1) +
#   scale_color_manual(values = cols) +
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
# p1
# 
# ggsave(plot = p1, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.png', width = 7.5, height = 6, units = 'in', dpi = 300)
# 


# plot p-value vs. position (gatk loci)
# only where nloci >= minloci
cols <- brewer.pal(4, 'Paired')[rep(1:2,12)]
p2 <- ggplot(datgatk[nloci >= minloci, ], aes(midPos, -log10(p), color = CHROM)) + 
    geom_point(size = 0.2, alpha = 0.3) +
    facet_wrap(~pop, ncol = 1) +
    scale_color_manual(values = cols) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey')
p2

ggsave(plot = p2, device = 'png', filename = 'figures/fst.siteshuffle.p_vs_pos.gatk.png', width = 7.5, height = 6, units = 'in', dpi = 300)


# plot p-value vs. nloci (gatk loci)
ggplot(datgatk, aes(nloci, -log10(p), color = pop)) + 
    geom_point(size = 0.2, alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') + 
    scale_x_log10()

#################
# print outliers
#################


datgatk[p < 0.05, ]

