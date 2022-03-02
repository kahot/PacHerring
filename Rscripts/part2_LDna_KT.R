########################################################################
library(LDna)
library(dplyr)
library(tidyr)
library(tibble)

########################################################################
# Read in data
ld<-read.table("Data/plink/pops/population_PWS91_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld")
colnames(ld)<-ld[1,]
ld<-ld[-1,]

ld1<-ld[ld$CHR_A=="1",]
ld1$R2<-as.numeric(ld1$R2)

#######
#create a pairwise LD matrix

l<-ld1[1:50,c("SNP_A","SNP_B","R2")]
lm<-acast(l,formula = SNP_A ~ SNP_B, value.var = 'R2',drop = F, fill=0)
lm<-data.frame(lm)
colnames(lm)<-gsub("\\.","\\:", colnames(lm))

lm<-lm[order(rownames(lm)),]
lm<-lm[,order(colnames(lm))]

for (i in 1:nrow(lm)){
    rown<-rownames(lm)[i]
    coln<-colnames(lm)[i]
    if(rown==coln) next
    if(rown!=coln) {
        if (!(rown %in% colnames(lm))){
            lm<-add_column(lm, newcol=rep(0,times=nrow(lm)), .before=i)
            colnames(lm)[colnames(lm)=="newcol"]<-rown}
    }
}

lm<-lm[order(rownames(lm)),]
lm<-lm[,order(colnames(lm))]

#Loop through coluumns
for (i in 1:ncol(lm)){
    rown<-rownames(lm)[i]
    coln<-colnames(lm)[i]
    if (!(coln %in% rownames(lm))){
        lm<-rbind(lm,rep(0,times=ncol(lm)))
        rownames(lm)[nrow(lm)]<-coln
    }
}

i=1
write.csv(lm, paste0("Output/LD/LD_pairwise_matrix_chr",i,',csv'))


#lower triangle is needed for LDna so transpose the matrix
l1m<-t(lm)

#ld1<-as.data.frame(ld1m)
#fileName <- "Output/LD/" 
#fileName2 <- "locus_metadata.txt" 


#ld.df <- read.delim(fileName)
#meta_df <- read.delim(fileName2)
ld.mat<-l1m
class(ld.mat) # check to make sure the data is loaded as a dataframe
ld.mat[1:5, 1:5]
########################################################################
# Run LDNa

# LDnaRaw function takes a lower diagonal matrix of pairwise LD values 
# and produces data files for subsequent LD network analyses.

ldna <- LDnaRaw(ld.mat, digits = 4, )
ldna$clusterfile #check output
ldna$stats

# Explore clustering in the data set
par(mfcol = c(1,3))

plotLDnetwork(ldna, LDmat = ld.mat, option = 1, threshold = 0.25)
plotLDnetwork(ldna, LDmat = ld.mat, option = 1, threshold = 0.50)
plotLDnetwork(ldna, LDmat = ld.mat, option = 1, threshold = 0.75)


# Use the extractClusters function to identify outlier clusters
# and plot cluster merging by LD threshold
# parameters: 
# ldna = ldna object
# min.edges = Minimum number of edges for a cluster that is shown as a branch in a tree (min # connections)
# Groups of loci are connected by "edges" that represent LD values above a given threshold
# phi = Controls λ_{lim} which sets the threshold above which λ are considered as outliers. Default is two, values below this are not recommended.
par(mfcol = c(1, 2))
clusters1 <- extractClusters(ldna, LDmat = ld.mat, min.edges = 5, lambda.lim = 0.0001)

# 10,000 foot level Summary table of the results
summary1 <- summaryLDna(ldna, clusters1, ld.mat)
summary1
plotLDnetwork(ldna, LDmat = ld.mat, option = 2, summary=summary1, threshold = 0.25)



# of course this program gives you a nested list oflists as output. Grrr.
# Now we have to parse this output.

# The assignment of loci to clusters is here
my_clust <- clusters1$clusters

# This is some code that will collapse the list of lists into a dataframe.
# Do not be alarmed by the NAs they are just placeholders to make the dataframe equal in length
temp_df <- data.frame(lapply(my_clust, "length<-", max(lengths(my_clust))), check.names = FALSE)
head(temp_df)

# Turn this into tidy format
results_df <- pivot_longer(data = temp_df, 
                     cols  = 1:ncol(temp_df), 
                     names_to = "cluster", 
                     values_to = "locus")

# get rid of the NAs
results_df <- results_df %>%
  drop_na()

head(results_df)

results_df <- left_join(results_df, meta_df, by = "locus")

summary1$cluster.name <- rownames(summary1)




#Explanation of output values of summaryLDna:
# lambda = the change in LD when two clusters merge. High values of lambda indicate the merger of large clusters of strongly associated loci
# Type specifies if a cluster is a "single outlier cluster", SOC or a "compound oulier cluster"
# Take a peek at the loci that are within the outlier groups
#"Merge.at" specfies the LD threshold for cluster merger. 
#"Median.LD" gives the median LD of all pairwise values between loci in a cluster.
#"MAD.LD" gives their unscaled median absolute deviation.


# save the results for plottings

write.table(results_df, file = "results_LDNA_batch_1_firstsnp_GCA900700415.e20f7.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

write.table(summary1, file = "summary_LDNA_batch_1_firstsnp_GCA900700415.e20f7.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

