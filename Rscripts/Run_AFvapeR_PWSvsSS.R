#AF-vapeR https://github.com/JimWhiting91/afvaper
library(afvaper,verbose = F)
library(vcfR)
library(stringr)
library(ggplot2)
library(gridExtra)


#Each chromosome need to be in a separate file
#Divide them up into chromosomes

popmap<-read.csv("Data/Sample_metadata_892pops.csv")
popmap<-popmap[popmap$pop=="PWS"|popmap$pop=="SS",]
popmap<-popmap[,c("Sample","Population.Year")]

#populations to compare
vector_list<-list(c("PWS91","PWS96"),
                  c("PWS91","PWS07"),
                  c("PWS91","SS96"),
                  c("PWS91","SS06"))
#just between 91-96
vector_list<-list(c("PWS91","PWS96"),
                  c("PWS91","SS96"))


names(vector_list)<-c("PWS96","SS96")



#Calculate Allele Frequency Change Vectors
# Set our window size
window_snps = 200


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

#saveRDS(all_chr_res, file="Output/all_chr_res_AFVapeR_PWSvsSS.RData")
saveRDS(all_chr_res, file="Output/all_chr_res_AFVapeR_PWSvsSS_96only.RData")

all_chr_res<-readRDS(file="Output/all_chr_res_AFVapeR_PWSvsSS.RData")


# To fetch all of  chr AFV, we take the first element of each list element
# Note: the merge_eigen_res() func is the same as unlist(,recursive=F)
AF_input <- merge_eigen_res(lapply(all_chr_res,'[[',1))

# All null, we take the second element of each list element
null_input <- merge_eigen_res(lapply(all_chr_res,'[[',2))

# We now have our whole genome's worth of AFV matrices and null matrices in a single input
c(head(names(AF_input)),tail(names(AF_input)))



#perform Eigen Analysis Over Allele Frequency Matrices
# The eigen_res output is a list containing per chromosome-region 1) eigenvalue distribution, 2) eigenvector loadings, and 
# 3) the projected A matrix that shows per-SNP scores for each eigenvector.


#In linear algebra, an eigenvector (characteristic vector of a linear transformation) is a nonzero vector that changes at most by a scalar factor 
#when that linear transformation is applied to it. The corresponding eigenvalue, often denoted by λ
# is the factor by which the eigenvector is scaled.


# Perform eigen analysis
eigen_res <- lapply(AF_input,eigen_analyse_vectors)

saveRDS(eigen_res, file="Output/eigeb_res_AFVapeR_PWSvsSS91-96.RData")


# View chromosomal regions: (1008)
head(names(eigen_res))

# View eigenvalue (λ) distribution of first matrix=first chromosomal period
eigen_res[[1]]$eigenvals
#Eigenvector_1 Eigenvector_2 Eigenvector_3 Eigenvector_4 
#    2.5561616     0.6352067     0.4335604     0.3750713 


# View eigenvector loadings of first matrix
eigen_res[[1]]$eigenvecs
#      Eigenvector_1 Eigenvector_2 Eigenvector_3 Eigenvector_4
#PWS96     0.5336672     0.1198154    0.01767346     0.8369774
#PWS07     0.5175362    -0.1225420   -0.79354127    -0.2956890
#SS96      0.4705731    -0.7006709    0.49341097    -0.2101592
#SS06      0.4753069     0.6925951    0.35570425    -0.4097193


# View head of SNP scores
head(eigen_res[[1]]$A_matrix)

#Find Null Cutoffs
#Using our null_input, we can output a matrix containing the cutoff expectations from the null distribution for each eigenvector for various significance thresholds:
    
# Get cutoffs for 95%, 99% and 99.9%  
# value for Eigenvector 2 is the sum of Eigenvalues 1 + 2, 
null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.95,0.99,0.999))
null_cutoffs
#                   95%      99%    99.9%
#Eigenvector 1 2.748449 2.866056 3.045489
#Eigenvector 2 3.293476 3.371688 3.526830
#Eigenvector 3 3.709700 3.750160 3.827790
#Eigenvector 4 4.000000 4.000000 4.000000

#Calculate empirical p-values
# calculate one-tailed p-values by comparing our observed eigenvalues to the null distribution for each window

# Calculate p-vals
pvals <- eigen_pvals(eigen_res,null_input)

# Show lowest pvals
head(pvals)
#                     Eigenvalue_1 Eigenvalue_2 Eigenvalue_3 Eigenvalue_4
#chr1:24227-590334      0.31060379   0.23981520   0.43134137  0.284584308
#chr1:592915-1303586    0.67170657   0.65977680   0.23108538  0.757994840
#chr1:1303603-1902327   0.92760145   0.83299334   0.66831663  0.624817504
#chr1:1903017-2548082   0.07900842   0.02574949   0.08583828  0.002929941
#chr1:2548509-2706025   0.54274915   0.41974161   0.04569909  0.757994840
#chr1:2711231-3280983   0.86676266   0.77327453   0.80156397  0.105317894



## Plot Eigenvalues Along Chromosomes ##

# Plot the raw eigenvalues, and visualise the cutoff of 99%
all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"99%"])

# Show the plots for eigenvalue 1
all_plots[[1]]
all_plots[[2]]
all_plots[[3]]



# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
# Empirical p-values' inclusion here is to aid in visualisation rather than to be used as explicit tests of significance, and changing the number of null permutations will change the empirical p-values:
all_plots_p <- eigenval_plot(eigen_res,null_vectors = null_input,plot.pvalues = T)

# Show the plots for eigenvalue 1
all_plots_p[[1]]
all_plots_p[[2]]

all_plots_p[[1]]+
    geom_hline(yintercept=1.3, color="dodgerblue1")+
    geom_hline(yintercept=2, color="red")


#export p-values
p.mat<-apply(pvals, 2, function(x) 10^(-x))
p.mat<-data.frame(p.mat)
write.csv(mat,"Output/AFvapeR_output/PWSvsSS_Pvalues.csv")


#All windows are named in the format chr:start-end: Can plot specific chromosomes using grep():
#Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc. for chr1
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

#99%
#$`Eigenvector 1`
#[1] "chr1:13501178-14186611"  "chr1:21044354-22086421"  "chr2:2962050-3824120"    "chr4:31245579-32246236" 
#[5] "chr5:11260732-12358256"  "chr7:23076749-24138551"  "chr8:9511398-10310547"   "chr10:310335-582470"    
#[9] "chr10:5016085-5372849"   "chr13:18335978-19518809" "chr14:13711187-14717832" "chr15:877103-1771586"   
#[13] "chr15:1772090-2425625"   "chr15:5914449-7392600"   "chr15:11752605-12633617" "chr15:12635237-13580922"
#[17] "chr15:13584367-14293476" "chr15:15054028-15918124" "chr15:15919271-16381740"
#
#$`Eigenvector 2`
#[1] "chr1:28966415-29992653"  "chr6:26726655-27359673"  "chr12:18576865-18887232" "chr12:22750862-23214025"
#[5] "chr15:9187945-9989916"   "chr20:1070169-2016781"   "chr22:850151-1309144"   
#
#$`Eigenvector 3`
#[1] "chr2:30915914-31498409"  "chr6:13825086-14681734"  "chr6:25139271-25505092"  "chr7:24713032-25385786" 
#[5] "chr12:21474149-22043017" "chr13:15800016-16219633" "chr15:29495-870222"      "chr15:4617672-5100185"  
#
#$`Eigenvector 4`
#[1] NA

# 2. 99.9%
##> significant_windows
#$`Eigenvector 1`
#[1] "chr1:21044354-22086421"  "chr4:31245579-32246236"  "chr10:5016085-5372849"   "chr13:18335978-19518809"
#[5] "chr15:5914449-7392600"   "chr15:15054028-15918124"
#
#$`Eigenvector 2`
#[1] "chr22:850151-1309144"




### Summarise Outliers
# which of our replicate vectors are associated with the relevant eigenvector?
# whether vectors are associated in the same (parallel) or different (antiparallel) direction?




# Summarise parallel evolution in windows that are significant on eigenvector 1
eig1_parallel <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.5,
                                              eigenvector = 1)
eig1_parallel1 <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 1)


# Show results
eig1_parallel
# 1. cutoff (correlation) =0.5. 99.9%
#                window_id eigenvector eigenvalue parallel_lineages parallel_pops antiparallel_pops
#1  chr1:21044354-22086421        Eig1   3.148683                 2    PWS96,SS96                  
#2  chr4:31245579-32246236        Eig1   3.091426                 2    PWS96,SS96                  
#3   chr10:5016085-5372849        Eig1   3.059516                 1         PWS07                  
#4 chr13:18335978-19518809        Eig1   3.097149                 2     SS96,SS06                  
#5   chr15:5914449-7392600        Eig1   3.135559                 2     SS96,SS06                  
#6 chr15:15054028-15918124        Eig1   3.111057                 2    PWS96,SS96          
            
# 1. cutoff (correlation) =0.3. 99.9%
#                window_id eigenvector eigenvalue parallel_lineages         parallel_pops antiparallel_pops
#1  chr1:21044354-22086421        Eig1   3.148683                 4 PWS96,PWS07,SS96,SS06                  
#2  chr4:31245579-32246236        Eig1   3.091426                 4 PWS96,PWS07,SS96,SS06                  
#3   chr10:5016085-5372849        Eig1   3.059516                 4 PWS96,PWS07,SS96,SS06                  
#4 chr13:18335978-19518809        Eig1   3.097149                 4 PWS96,PWS07,SS96,SS06                  
#5   chr15:5914449-7392600        Eig1   3.135559                 4 PWS96,PWS07,SS96,SS06                  
#6 chr15:15054028-15918124        Eig1   3.111057                 4 PWS96,PWS07,SS96,SS06                  


#cut off = 0.5, 99%
#                 window_id eigenvector eigenvalue parallel_lineages    parallel_pops antiparallel_pops
#1   chr1:13501178-14186611        Eig1   2.902749                 3 PWS96,PWS07,SS06                  
#2   chr1:21044354-22086421        Eig1   3.148683                 2       PWS96,SS96                  
#3     chr2:2962050-3824120        Eig1   3.029637                 2       PWS96,SS96                  
#4   chr4:31245579-32246236        Eig1   3.091426                 2       PWS96,SS96                  
#5   chr5:11260732-12358256        Eig1   2.866199                 2       PWS96,SS96                  
#6   chr7:23076749-24138551        Eig1   2.896856                 2       PWS96,SS96                  
#7    chr8:9511398-10310547        Eig1   2.886182                 2      PWS96,PWS07                  
#8      chr10:310335-582470        Eig1   2.939674                 3  PWS07,SS96,SS06                  
#9    chr10:5016085-5372849        Eig1   3.059516                 1            PWS07                  
#10 chr13:18335978-19518809        Eig1   3.097149                 2        SS96,SS06                  
#11 chr14:13711187-14717832        Eig1   2.883277                 2       PWS96,SS96                  
#12    chr15:877103-1771586        Eig1   2.920033                 2       PWS96,SS96                  
#13   chr15:1772090-2425625        Eig1   2.881083                 1             SS96                  
#14   chr15:5914449-7392600        Eig1   3.135559                 2        SS96,SS06                  
#15 chr15:11752605-12633617        Eig1   2.966010                 3  PWS96,SS96,SS06                  
#16 chr15:12635237-13580922        Eig1   2.996318                 2       PWS96,SS96                  
#17 chr15:13584367-14293476        Eig1   2.962886                 3  PWS96,SS96,SS06                  
#18 chr15:15054028-15918124        Eig1   3.111057                 2       PWS96,SS96                  
#19 chr15:15919271-16381740        Eig1   2.936024                 2       PWS96,SS96   

eig1_parallel1
eig1<-merge(eig1_parallel, eig1_parallel1[,c(1,5)], by="window_id")


eig2_parallel <- summarise_window_parallelism(window_id = significant_windows[[2]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.5,
                                              eigenvector = 2)

eig2_parallel2 <- summarise_window_parallelism(window_id = significant_windows[[2]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 2)

eig2_parallel

write.csv(eig2_parallel, "Output/AFvapeR_output/PWSvsSS_eig2_parallel_summary_0.5_99.csv")
write.csv(eig2_parallel2, "Output/AFvapeR_output/PWSvsSS_eig2_parallel_summary_0.3_99.csv")

colnames(eig1)[c(4:5,7)]<-c("parallel_lin_0.5","eig1.0.5", "eig1.0.3")

eig12<-merge(eig1, eig2_parallel[c(1,5:7)], by="window_id", all=T)

# Show results
head(eig2_parallel)

## Explore candidate regions ##
#Find a focal candidate region and whether the signal is localised.
#The per-SNP scores are stored within the 'A_matrix' slot of each window’s entry in eigen_res. This a matrix with a row per SNP and a column per eigenvector, given the per-SNP per-eigenvector score of association:


# Fetch the A matrix for significant windows,
# lookat eig2_parallel 2

windws<-eig2_parallel2$window_id[c(1,3,5,7,9,11,13)]
plots<-list()
for (i in windws){
    no<-which(names(eigen_res)==i)
    A_mat <- eigen_res[[no]]$A_matrix
    
    to_plot <- data.frame(snp=rownames(A_mat),
                          eig1_score=A_mat[,1])
    
    to_plot <- to_plot %>% tidyr::separate("snp",into=c("chr","pos"),sep="_")
    to_plot$pos <- as.integer(to_plot$pos)
    
    plots[[i]]<-ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
        geom_point(color="darkblue")+
        labs(y="Eig1 Score",x="Pos (bp)")+
        theme_bw()+
        ggtitle(i)
    write.csv(to_plot, paste0("Output/AFvapeR_output/Parallel_PWS96_SS96/PWSvsSS_eig12_pararell_region_",i,".csv"))
}
pdf("Output/AFvapeR_output/Parallel_PWS96_SS96/PWSvsSS_eig12_Pararell_regions.pdf", height = 8, width = 8)
do.call(grid.arrange, c(plots, ncol=2))
dev.off()


A_mat <- eigen_res[[12]]$A_matrix

to_plot <- data.frame(snp=rownames(A_mat),
                      eig1_score=A_mat[,1])

to_plot <- to_plot %>% tidyr::separate("snp",into=c("chr","pos"),sep="_")
to_plot$pos <- as.integer(to_plot$pos)

ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
    geom_point(color="darkblue")+
    labs(y="Eig1 Score",x="Pos (bp)")+
    theme_bw()+
    ggtitle("chr12")


### What gene/region are these high eig1 score snp belong to?

pwsss<-list.files("Output/AFvapeR_output/Parallel_PWS96_SS96/", pattern= "^PWSvsSS_eig12.+.csv$")
#chromosome 1,2,4,5,7,14,15
for (i in 1: length(pwsss)){
    eig<-read.csv(paste0("Output/AFvapeR_output/Parallel_PWS96_SS96/", pwsss[i]))
    ch<-eig$chr[1]
    ch<-gsub("chr",'',ch)
    #read the chromosome gff file
    gff<-read.table(paste0("Data/annotations/chr/",ch), sep = "\t")
    
    #select snp with Eig1 score >=0.4
    eig<-eig[abs(eig$eig1_score)>=0.4,]
    eig<-eig[order(abs(eig$eig1_score), decreasing = T),]
    geneinfo<-data.frame()
    for (j in 1: nrow(eig)){
        pos<-eig$pos[j]
        gene<-gff[which(gff$V4>=pos)[1]:(which(gff$V4>=pos)[1]+3),]
        gene$SNP<-pos
        geneinfo<-rbind(geneinfo, gene)
    }    
    fname<-gsub("PWSvsSS_eig1scores_region_eig1_","", pwsss[i])
    fname<-gsub(".csv","",fname)
    write.csv(geneinfo, paste0("Output/AFvapeR_output/Parallel_PWS96_SS96/Anno_PWSvsSS_eig12_SNPs_over0.4_region_",fname,".csv"))  
}

#Extract the positions
int.pos<-data.frame()
for (i in 1: length(pwsss)){
    eig<-read.csv(paste0("Output/AFvapeR_output/Parallel_PWS96_SS96/", pwsss[i]))
    ch<-eig$chr[1]
    ch<-gsub("chr",'',ch)
    #select snp with Eig1 score >=0.4
    eig<-eig[abs(eig$eig1_score)>=0.4,]
    
    int.pos<-rbind(int.pos,eig)
}

#create bed file form int.pos 
bed1<-int.pos[,c("chr","pos")]
bed1$end<-bed1$pos
write.table(bed1, "Output/AFvapeR_output/PWS_SS_eig12_Pararell_candidate_positions.bed", sep = "\t", quote = F,row.names = F, col.names = F)

#
#bed1<-int.pos[,c("chr","pos")]
#bed1$start<-bed1$pos-10
#bed1$end<-bed1$pos+10
#write.csv(int.pos, "Output/AFvapeR_output/PWS_SS_eig1_candidate_positions.csv")


### Look at eig1 loci
windws<-eig1_parallel$window_id
Plots<-list()
for (i in windws){
    no<-which(names(eigen_res)==i)
    A_mat <- eigen_res[[no]]$A_matrix
    
    to_plot <- data.frame(snp=rownames(A_mat),
                          eig1_score=A_mat[,1])
    
    to_plot <- to_plot %>% tidyr::separate("snp",into=c("chr","pos"),sep="_")
    to_plot$pos <- as.integer(to_plot$pos)
    
    Plots[[i]]<-ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
        geom_point(color="darkblue")+
        labs(y="Eig1 Score",x="Pos (bp)")+
        theme_bw()+
        ggtitle(i)
    write.csv(to_plot, paste0("Output/AFvapeR_output/PWSvsSS_eig1Pararell_region_",i,".csv"))
}

pdf("Output/AFvapeR_output/Parallel_PWS96_SS96/PWSvsSS_eig1Pararell_regions.pdf", height = 20, width = 12)
do.call(grid.arrange, c(Plots, ncol=3))
dev.off()

### What gene/region are these high eig1 score snp belong to?

par1<-list.files("Output/AFvapeR_output/", pattern= "^PWSvsSS_eig1Pararell.+.csv$")

for (i in 1: length(par1)){
    eig<-read.csv(paste0("Output/AFvapeR_output/", par1[i]))
    ch<-eig$chr[1]
    ch<-gsub("chr",'',ch)
    #read the chromosome gff file
    gff<-read.table(paste0("Data/annotations/chr/",ch), sep = "\t")
    
    #select snp with Eig1 score >=0.3
    eig<-eig[abs(eig$eig1_score)>=0.3,]
    eig<-eig[order(abs(eig$eig1_score), decreasing = T),]
    geneinfo<-data.frame()
    for (j in 1: nrow(eig)){
        pos<-eig$pos[j]
        gene<-gff[which(gff$V4>=pos)[1]:(which(gff$V4>=pos)[1]+3),]
        gene$SNP<-pos
        geneinfo<-rbind(geneinfo, gene)
    }    
    fname<-gsub("PWSvsSS_eig1scores_region_eig1_","", pwsss[i])
    fname<-gsub(".csv","",fname)
    write.csv(geneinfo, paste0("Output/AFvapeR_output/Parallel_PWS96_SS96/Anno_eig1_pararell_SNPs_over0.3_region_",fname,".csv"))  
}

#Extract the positions
int.pos<-data.frame()
for (i in 1: length(par1)){
    eig<-read.csv(paste0("Output/AFvapeR_output/", par1[i]))
    ch<-eig$chr[1]
    ch<-gsub("chr",'',ch)
    #select snp with Eig1 score >=0.3
    eig<-eig[abs(eig$eig1_score)>=0.3,]
    
    int.pos<-rbind(int.pos,eig)
}

#create bed file form int.pos 
bed1<-int.pos[,c("chr","pos")]
bed1$end<-bed1$pos
write.table(bed1, "Output/AFvapeR_output/PWS_SS_eig1_pararell_candidate_positions.bed", sep = "\t", quote = F,row.names = F, col.names = F)


## Read the AF of candidate loci 

afiles<-list.files("Output/AF/",pattern = "freq")
AFlists<-list()
afreqs<-data.frame()
for (i in 1: length(afiles)){
    af<-read.table(paste0("Output/AF/",afiles[i]), sep="\t", row.names = NULL)
    colnames(af)[1:4]<-colnames(af)[2:5]
    colnames(af)[5:6]<-c("Freq1","Freq2")
    df<-af[af$POS %in% int.pos$pos,]
    df$Pop<-pops[i]
    afreqs<-rbind(afreqs,df) 
    #AFlists[[i]]<-df
    #names(AFlists)[i]<-pops[i]
}

afreqs$freq<-substr(afreqs$Freq2, 3,11)
afreqs$freq<-as.numeric(afreqs$freq)


#look at PWS and SS
poplist1<-c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17")
af1<-afreqs[afreqs$Pop %in% poplist1,]

positions<-unique(af1$POS)
po1<-positions[1:30]

posC<-afreqs[afreqs$POS==496434,]
posC$location<-substr(posC$Pop, 1,2)
posC$Pop<-factor(posC$Pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))

ggplot(posC, aes(x=Pop, y=freq, fill=location))+
    geom_bar(stat="identity")

7355172

pos2<-afreqs[afreqs$POS==7355172,]
pos2$location<-substr(pos2$Pop, 1,2)
pos2$Pop<-factor(pos2$Pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))

ggplot(pos2, aes(x=Pop, y=freq, fill=location))+
    geom_bar(stat="identity")


#first 30 loci
ggplot(af1[af1$POS %in%po1, ], aes(x=factor(POS), y=freq, fill=Pop))+
    geom_bar(stat="identity", position=position_dodge())
