# From https://github.com/pinskylab/codEvol calcNe_ANGSD.r


### Jorde & Ryman/NeEstimator approach
# Jorde & Ryman 2007

# Ne in # diploid individuals
# based on NeEstimator manual v2.1 
jrNe2 <- function(maf1, maf2, n1, n2, gen){
    Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles
    
    z <- (maf1+maf2)/2 # for the first allele
    z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
    Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
    Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual
    
    sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals
    
    S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
    S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
    Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
    return(gen/(2*Fsprime)) # calculation of Ne in # diploid individuals
}

# bootstrap over loci to get CIs
# Jorde & Ryman Ne estimator, for boot() to use
jrNe2boot <- function(data, gen, indices){
    maf1 <- data$freq1[indices]
    maf2 <- data$freq2[indices]
    n1 <- data$nInd1[indices]
    n2 <- data$nInd2[indices]
    
    Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles
    
    z <- (maf1+maf2)/2 # for the first allele
    z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
    Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
    Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual
    
    sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals
    
    S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
    S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
    Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
    Ne <- gen/(2*Fsprime)
    if(Ne < 0) Ne <- Inf
    return(Ne) # calculation of Ne in # diploid individuals
}

# for block bootstrapping across LGs
# frin https://stackoverflow.com/questions/11919808/block-bootstrap-from-subject-list
jrNe2block <- function(lgs, gen, alldata, indices){
    mydata <- do.call("rbind", lapply(indices, function(n) subset(alldata, chromo==lgs[n])))
    return(jrNe2boot(mydata, gen, indices = 1:nrow(mydata)))
}
