#https://github.com/pinskylab/codEvol 

# shuffle theta values across sites and calculate windowed averages to get a null distribution of max genome-wide thetas

#################################
# parameters
#################################
winsz <- 50000 # window size
winstp <- 10000 # window step
nrep <- 1000 # number of reshuffles

####################
# load functions
require(data.table)

# after https://github.com/ANGSD/angsd/blob/master/misc/stats.cpp
a1f <- function(nsam) return(sum(1/seq(1, nsam-1)))
a2f <- function(nsam) return(sum(1/(seq(1, nsam-1)*seq(1, nsam-1))))
b1f <- function(nsam) return((nsam + 1)/(3*(nsam-1)))
b2f <- function(nsam) return((2*(nsam*nsam + nsam + 3))/(9*nsam*(nsam - 1)))
c1f <- function(a1, b1) return(b1 - (1/a1))
c2f <- function(nsam, a1, a2, b2) return(b2 - ((nsam + 2)/(a1*nsam)) + (a2/(a1 * a1)))
e1f <- function(a1, c1) return(c1/a1)
e2f <- function(a1, a2, c2) return(c2/((a1*a1) + a2))

# Tajima's D calculation
# nsam: sample size
# thetaW: Watterson's theta (# segregating sites / a1)
# sumk: theta pi (average number of SNPs in pairwise comparisons)
# after https://github.com/ANGSD/angsd/blob/master/misc/stats.cpp
tajd <- function(nsam, thetaW, sumk){
    a1 <- a1f(nsam)
    segsites <- thetaW * a1
    if(segsites == 0) return(0)
    a2 <- a2f(nsam)
    b1 <- b1f(nsam)
    b2 <- b2f(nsam)
    c1 <- c1f(a1, b1)
    c2 <- c2f(nsam, a1, a2, b2)
    e1 <- e1f(a1, c1)
    e2 <- e2f(a1, a2, c2)
    res <- (sumk - (thetaW))/sqrt((e1*segsites) + ((e2*segsites)*(segsites-1)))
    return(res)
}

# calculate thetas in windows from data
calcstat <- function(temp1, temp2, nchr1, nchr2){
    # create new columns as indices for windows
    # need multiple columns because overlapping windows
    for(j in 1:(winsz/winstp)){
        temp1[, (paste0('win_', j)) := floor((Pos - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
        temp2[, (paste0('win_', j)) := floor((Pos - (j-1)*winstp)/winsz)*winsz + winsz/2 + (j-1)*winstp]
    }
    
    # calc ave theta for each window
    for(j in 1:(winsz/winstp)){
        nm1 <- paste0('win_', j)
        if(j ==1){
            bins1 <- temp1[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW1 = sum(exp(Watterson), na.rm = TRUE),
                                                                                       tP1 = sum(exp(Pairwise), na.rm = TRUE)), 
                                                                                   by = .(Chromo, WinCenter)]
            bins2 <- temp2[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW2 = sum(exp(Watterson), na.rm = TRUE),
                                                                                       tP2 = sum(exp(Pairwise), na.rm = TRUE)), 
                                                                                   by = .(Chromo, WinCenter)]
        }
        if(j > 1){
            bins1 <- rbind(bins1, temp1[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW1 = sum(exp(Watterson), na.rm = TRUE),
                                                                                                    tP1 = sum(exp(Pairwise), na.rm = TRUE)), 
                                                                                                by = .(Chromo, WinCenter)])
            bins2 <- rbind(bins2, temp2[, .(Chromo, Watterson, Pairwise, WinCenter = get(nm1))][, .(tW2 = sum(exp(Watterson), na.rm = TRUE),
                                                                                                    tP2 = sum(exp(Pairwise), na.rm = TRUE)), 
                                                                                                by = .(Chromo, WinCenter)])
        }
    }
    
    # calculate Tajima's D
    bins1[, tD1 := tajd(nchr1, tW1, tP1), by = 1:nrow(bins1)]
    bins2[, tD2 := tajd(nchr2, tW2, tP2), by = 1:nrow(bins2)]
    
    # merge by window and calculate change
    bins <- merge(bins1, bins2, by = c("Chromo", "WinCenter"))
    bins[, ':='(tWd = tW2 - tW1, tPd = tP2 - tP1, tDd = tD2 - tD1)]
    
    return(bins)
}

# calculate min and max change from reshuffled data
shufflestat <- function(dat1, dat2, nchr1, nchr2){
    inds1 <- sample(1:nrow(dat1), nrow(dat1), replace = FALSE)
    temp1 <- cbind(dat1[, .(Chromo, Pos)], dat1[inds1, .(Watterson, Pairwise)]) # shuffle theta across position pairs
    inds2 <- sample(1:nrow(dat2), nrow(dat2), replace = FALSE)
    temp2 <- cbind(dat2[, .(Chromo, Pos)], dat2[inds2, .(Watterson, Pairwise)]) # second time-point
    
    bins <- calcstat(temp1, temp2, nchr1, nchr2)
    
    # return the max and min windowed change
    # exclude windows with negative midpoints
    return(bins[WinCenter >= winsz/2, .(mintWd = min(tWd, na.rm = TRUE), 
                                        maxtWd = max(tWd, na.rm = TRUE),
                                        mintPd = min(tPd, na.rm = TRUE), 
                                        maxtPd = max(tPd, na.rm = TRUE),
                                        mintDd = min(tDd, na.rm = TRUE), 
                                        maxtDd = max(tDd, na.rm = TRUE))])
    
}


######################
# Load data
######################

years<-c("PWS91","PWS96","PWS07","PWS17")
comb<-t(combn(years, 2))
popsize<-read.csv("Data/popinfo/popsize.csv")

for (p in 1: nrow(comb)){
    dat1<-fread(paste0('Data/new_vcf/angsd/fromVCF/theta/',comb[p, 1],'.theta.mini'))
    dat2<-fread(paste0('Data/new_vcf/angsd/fromVCF/theta/',comb[p, 2],'.theta.mini'))
    
    nchr1 <- popsize$Freq[popsize$Var1==comb[p, 1]]
    nchr2 <- popsize$Freq[popsize$Var1==comb[p, 2]]
    outfilesuff <- paste0(comb[p, 1],"_" ,comb[p, 2])

    setnames(dat1, '#Chromo', 'Chromo',skip_absent=T)
    setnames(dat2, '#Chromo', 'Chromo',skip_absent=T)
    
    
    #######################################
    # Calculate change on unshuffled data
    #######################################
    
    # full dataset (include all loci)
    out <- calcstat(dat1, dat2, nchr1, nchr2)
    outfile <- paste0('/home/ktist/ph/data/angsd/theta/theta_change_region_', winsz, '.', outfilesuff, '.csv.gz')
    write.csv(out, file = gzfile(outfile))
    
    ################################
    # Run reshuffle calculations
    # shuffle and recalc windowed thetas
    ################################
    
    for(i in 1:nrep){
        cat(i); cat(' ')
        tempminmax <- shufflestat(dat1, dat2, nchr1, nchr2)
        
        if(i == 1) minmax <- tempminmax # save the results from this iteration
        if(i > 1) minmax <- rbind(minmax, tempminmax)
    }
    
    write.csv(minmax, gzfile(paste0('/home/ktist/ph/data/angsd/theta/theta.siteshuffle.', outfile, '.csv.gz')), row.names = FALSE) # write out all iterations
    
    rm(minmax)
    
}


