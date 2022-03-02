#AF-vapeR https://github.com/JimWhiting91/afvaper
# remotes::install_github("JimWhiting91/afvaper")
library(afvaper,verbose = F)
library(vcfR)
library(stringr)
#each of the plots is a ggplot object that can be extracted and edited however you like, for e.g. we can remove the title and change the colour using standard ggplot syntax:
library(ggplot2)


#Each chromosome need to be in a separate file
#Divide them up into chromosomes

#to Zip VCF file, use bgzip
#bgzip -c PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf > PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz

# Index the .vcf.gz
#tabix -p vcf PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz

#create a bash script to divide them into 26 chromosomes
sink("Data/vcfs/PWS_SS/Subset_by_chromosome.sh")
cat("#!/bin/bash\n")
#26 chromosomes
for (i in 1: 26){
    ch<-paste0("chr",i)
    phrase<-paste0('bcftools view PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz -r ',ch,' -Oz --threads 10 > PWS_SS_',ch,".vcf.gz")
    cat(phrase)
    cat("\n")
}
sink(NULL)


#Start with chr1

vcf_in <- read.vcfR("Data/vcfs/PWS_SS/PWS_SS_chr1.vcf.gz")
vcf_in

popmap<-read.csv("Data/Sample_metadata_892pops.csv")
popmap<-popmap[popmap$pop=="PWS"|popmap$pop=="SS",]
popmap<-popmap[,c("Sample","Population.Year")]

#populations to compare
vector_list<-list(c("PWS91","PWS96"),
              c("PWS91","PWS07"),
              c("PWS91","PWS17"),
              c("PWS91","SS06"))
names(vector_list)<-c("PWS96","PWS07","PWS17","SS06")

#Calculate Allele Frequency Change Vectors
# Set our window size
window_snps = 200

# Calculate Allele Frequency Change Vector Matrices
AF_input <- calc_AF_vectors(vcf = vcf_in,
                            window_size = window_snps,
                            popmap = popmap,
                            vectors = vector_list,
                            n_cores = 10)

# Show features of input...
print(paste0("Number of windows = ",length(AF_input)))
print(paste0("Number of SNPs per window = ",ncol(AF_input[[1]])))
print(paste0("Number of vectors per window = ",nrow(AF_input[[1]])))

# How many permutations to run
null_perm_N = 1000

# Calculate Allele Frequency Change Vector Matrices
null_input <- calc_AF_vectors(vcf = vcf_in,
                              window_size = window_snps,
                              popmap = popmap,
                              vectors = vector_list,
                              n_cores = 4,
                              null_perms = null_perm_N)

print(paste0("Number of null windows = ",length(null_input)))
print(paste0("Number of SNPs per window = ",ncol(null_input[[1]])))
print(paste0("Number of vectors per window = ",nrow(null_input[[1]])))


## Now run the all 26 chromosomes
#Divide the desired total of 10,000 (or more) up between chromosomes based on the relative sizes (this info is available in a genome fasta index for e.g.).
# calculate the approximate chromosome lengths
map<-read.table("Data/vcfs/PWS_SS/PWS_SS_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.map")

chrsize<-data.frame(chr=paste0("chr",1:26))
for (i in 1:26){
    df<-map[map$V1==i,]
    begin<-df$V4[1]-1
    end<-df$V4[nrow(df)]
    length<-end-begin
    chrsize$length[i]<-length
}


# How many permutations do we want in total?
total_perms <- 100000

# Fetch proportional size of all chromosomes
chr_props <- chrsize$length/sum(chrsize$length)
chr_perms <- data.frame(chr=chrsize$chr,
                        perms=round(chr_props * total_perms))

# This gives us approximately 10000 null perms in total, distributed across the genome according to relative size of chromosomes...
print(chr_perms)

# We can now loop over our chromosome VCFs.
vfiles<-list.files("Data/vcfs/PWS_SS/", pattern="_chr.*.vcf.gz$")

vcf_list<-list()
for (i in 1:26){
    vcf_list[[i]]<-read.vcfR(paste0("Data/vcfs/PWS_SS/PWS_SS_chr",i,".vcf.gz"))
    names(vcf_list)[i]<-paste0("chr",i)
}    


#different comparisons
vector_list<-list(c("PWS91","PWS96"),
                  c("PWS91","PWS07"),
                  c("PWS91","SS96"),
                  c("PWS91","SS06"))
names(vector_list)<-c("PWS96","PWS07","SS96","SS06")




# Run the simulations of a total 100k devided up based on the chromosomal size  
# Increase the permutation later

all_chr_res <- lapply(1:26,function(i){
    
    # First read in the VCF
    #chr_vcf <- read.vcfR(paste0("Data/vcfs/PWS_SS/PWS_SS_chr",i,".vcf.gz"))
    chr_vcf <- vcf_list[[i]]
    
    # Calculate AFV
    chr_AF_input <- calc_AF_vectors(vcf = chr_vcf,
                                    window_size = window_snps,
                                    popmap = popmap,
                                    vectors = vector_list,
                                    n_cores = 4)
    
    # Calculate null AFV
    chr_null_input <- calc_AF_vectors(vcf = chr_vcf,
                                      window_size = window_snps,
                                      popmap = popmap,
                                      vectors = vector_list,
                                      n_cores = 4,
                                      null_perms = chr_perms$perms[i])
    
    ## We could save these to some temporary file, e.g.
    # saveRDS(list(chr_AF_input,chr_null_input),paste0("chr",i,"_AFV.rds"))
    
    # Return our results
    return(list(chr_AF_input,chr_null_input))
})

saveRDS(all_chr_res, file="Output/all_chr_res_AFVapeR.RData")

# To fetch all of  chr AFV, we take the first element of each list element
# Note: the merge_eigen_res() func is the same as unlist(,recursive=F)
AF_input <- merge_eigen_res(lapply(all_chr_res,'[[',1))

# All null, we take the second element of each list element
null_input <- merge_eigen_res(lapply(all_chr_res,'[[',2))

# We now have our whole genome's worth of AFV matrices and null matrices in a single input
c(head(names(AF_input)),tail(names(AF_input)))



#perform Eigen Analysis Over Allele Frequency Matrices
# Perform eigen analysis
eigen_res <- lapply(AF_input,eigen_analyse_vectors)

# View chromosomal regions:
head(names(eigen_res))

# View eigenvalue distribution of first matrix
eigen_res[[1]]$eigenvals
#Eigenvector_1 Eigenvector_2 Eigenvector_3 Eigenvector_4 
#    2.4505478     0.7298721     0.4741240     0.3454561 


# View eigenvector loadings of first matrix
eigen_res[[1]]$eigenvecs
#      Eigenvector_1 Eigenvector_2 Eigenvector_3 Eigenvector_4
#PWS96    -0.5211819     0.4081990    -0.3810525    -0.6454006
#PWS07    -0.5403778    -0.1282842    -0.5103588     0.6565583
#PWS17    -0.4345403    -0.8036379     0.2612301    -0.3116081
#SS06     -0.4975298     0.4136215     0.7253218     0.2351374


# View head of SNP scores
head(eigen_res[[1]]$A_matrix)

#Find Null Cutoffs
#Using our null_input, we can output a matrix containing the cutoff expectations from the null distribution for each eigenvector for various significance thresholds:
    
# Get cutoffs for 95%, 99% and 99.9%
null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.95,0.99,0.999))
null_cutoffs
#                   95%      99%    99.9%
#Eigenvector 1 2.693177 2.815598 2.981889
#Eigenvector 2 3.253051 3.334198 3.461600
#Eigenvector 3 3.689391 3.731468 3.801997
#Eigenvector 4 4.000000 4.000000 4.000000

#Calculate empirical p-values
# Calculate p-vals
pvals <- eigen_pvals(eigen_res,null_input)

# Show lowest pvals
head(pvals)
#                     Eigenvalue_1 Eigenvalue_2 Eigenvalue_3 Eigenvalue_4
#chr1:24227-590334       0.4149117    0.1603468   0.14553709   0.61653767
#chr1:592915-1303586     0.5901982    0.4535409   0.18963621   0.27592448
#chr1:1303603-1902327    0.6721466    0.6000480   0.39206216   0.61653767
#chr1:1903017-2548082    0.1599768    0.1737165   0.09907802   0.01901962
#chr1:2548509-2706025    0.5359993    0.5598188   0.52895942   0.27592448
#chr1:2711231-3280983    0.8693726    0.8386732   0.57957841   0.49612008


## Plot Eigenvalues Along Chromosomes ##

# Plot the raw eigenvalues, and visualise the cutoff of 99%
all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"99%"])

# Show the plots for eigenvalue 1
all_plots[[1]]
all_plots[[2]]


# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
all_plots_p <- eigenval_plot(eigen_res,null_vectors = null_input,plot.pvalues = T)

# Show the plots for eigenvalue 1
all_plots_p[[1]]
all_plots_p[[2]]

#We can also exploit that all windows are named in the format chr:start-end to plot specific chromosomes using grep():
#Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
chr1_windows <- grep("chr1:",names(eigen_res))
all_plots_p_chr1 <- eigenval_plot(eigen_res[chr1_windows],null_vectors = null_input,plot.pvalues = T)

# Show the plots for eigenvalue 1
all_plots_p_chr1[[1]]


chr7_windows <- grep("chr7",names(eigen_res))
all_plots_p_chr7 <- eigenval_plot(eigen_res[chr7_windows],null_vectors = null_input,plot.pvalues = T)

# Show the plots for eigenvalue 1
all_plots_p_chr7[[1]]


chr12_windows <- grep("chr12",names(eigen_res))
all_plots_p_chr12 <- eigenval_plot(eigen_res[chr12_windows],null_vectors = null_input,plot.pvalues = T)
# Show the plots for eigenvalue 1
all_plots_p_chr12[[1]]

#each of the plots is a ggplot object that can be extracted and edited however you like, for e.g. we can remove the title and change the colour using standard ggplot syntax:
library(ggplot2)

# Pull the figure 
eig1_pval_fig <- all_plots_p[[1]]

# Edit
eig1_pval_fig + theme(title = element_blank()) + geom_step(colour="red2")


#Pull Significant Windows
# Recall the use of find_null_cutoffs() to fetch a matrix of cutoffs...
# null_cutoffs

# Find significant windows above 99.9% null permutation
significant_windows <- signif_eigen_windows(eigen_res,null_cutoffs[,"99%"])

# Display 'outliers'
significant_windows

# 1. 99.9%
#$`Eigenvector 1`
#[1] "chr1:21044354-22086421" "chr2:2962050-3824120"   "chr10:5016085-5372849" 

#$`Eigenvector 2`
#[1] "chr22:850151-1309144"

# 2. 99%
#$`Eigenvector 1`
#[1] "chr1:21044354-22086421" "chr1:28966415-29992653" "chr2:2962050-3824120"   "chr4:31245579-32246236"
#[5] "chr7:17257596-17698794" "chr8:9511398-10310547"  "chr10:310335-582470"    "chr10:5016085-5372849" 
#[9] "chr13:3504785-3976096" 
#
#$`Eigenvector 2`
#[1] "chr1:13501178-14186611"  "chr1:17853682-18933195"  "chr7:1242406-1775673"    "chr12:18576865-18887232"
#[5] "chr15:5914449-7392600"   "chr15:12635237-13580922" "chr15:15054028-15918124" "chr16:22149558-23232529"
#[9] "chr22:850151-1309144"   
#
#$`Eigenvector 3`
#[1] "chr1:8617410-9178692"    "chr4:29579426-30239739"  "chr7:23076749-24138551"  "chr9:23056244-23546536" 
#[5] "chr12:22750862-23214025" "chr12:25314629-26276329" "chr14:12133360-12738993" "chr15:13584367-14293476"
#[9] "chr24:2799555-3495911"  
#
#$`Eigenvector 4`
#[1] "chr3:1708712-2088369"    "chr14:11418754-12133359"



### Summarise Outliers
# Summarise parallel evolution in windows that are significant on eigenvector 1
eig1_parallel <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.5,
                                              eigenvector = 1)
# Show results
head(eig1_parallel)
# 1. cutoff (correlation) =0.5 
#               window_id eigenvector eigenvalue parallel_lineages    parallel_pops antiparallel_pops
#1 chr1:21044354-22086421        Eig1   3.022109                 2       PWS07,SS06                  
#2 chr1:28966415-29992653        Eig1   2.835110                 3 PWS96,PWS07,SS06                  
#3   chr2:2962050-3824120        Eig1   3.045063                 2      PWS96,PWS17                  
#4 chr4:31245579-32246236        Eig1   2.931434                 1            PWS96                  
#5 chr7:17257596-17698794        Eig1   2.936297                 2      PWS96,PWS17                  
#6  chr8:9511398-10310547        Eig1   2.871694                 2      PWS96,PWS07      

# 1. cutoff (correlation) =0.3 
#               window_id eigenvector eigenvalue parallel_lineages          parallel_pops antiparallel_pops
#1 chr1:21044354-22086421        Eig1   3.022109                 4 PWS96,PWS07,PWS17,SS06                  
#2 chr1:28966415-29992653        Eig1   2.835110                 4 PWS96,PWS07,PWS17,SS06                  
#3   chr2:2962050-3824120        Eig1   3.045063                 4 PWS96,PWS07,PWS17,SS06                  
#4 chr4:31245579-32246236        Eig1   2.931434                 4 PWS96,PWS07,PWS17,SS06                  
#5 chr7:17257596-17698794        Eig1   2.936297                 4 PWS96,PWS07,PWS17,SS06                  
#6  chr8:9511398-10310547        Eig1   2.871694                 4 PWS96,PWS07,PWS17,SS06                  
#7    chr10:310335-582470        Eig1   2.843782                 4 PWS96,PWS07,PWS17,SS06                  
#8  chr10:5016085-5372849        Eig1   3.008630                 4 PWS96,PWS07,PWS17,SS06                  
#9  chr13:3504785-3976096        Eig1   2.820830                 4 PWS96,PWS07,PWS17,SS06                 




## Explore candidate regions ##
#Find a focal candidate region and whether the signal is localised.
#The per-SNP scores are stored within the 'A_matrix' slot of each window’s entry in eigen_res. This a matrix with a row per SNP and a column per eigenvector, given the per-SNP per-eigenvector score of association:

# Fetch an A matrix for chr1
A_mat <- eigen_res[[1]]$A_matrix
head(A_mat)

to_plot <- data.frame(snp=rownames(A_mat),
                      eig1_score=A_mat[,1])

to_plot <- to_plot %>% tidyr::separate("snp",into=c("chr","pos"),sep="_")
to_plot$pos <- as.integer(to_plot$pos)

ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
    geom_point(color="darkblue")+
    labs(y="Eig1 Score",x="Pos (bp)")+
    theme_bw()
ggsave("Output/AFvapeR_output/Eig1Score_alongGenome_chr1.pdf", height = 4, width = 7)
