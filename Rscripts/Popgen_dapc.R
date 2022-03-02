
ind<-read.table("~/programs/gdc/ph_filtered_snps_ind")
ind$Pop<-substr(ind$V1,6,10)
ind<-ind[,c(1,2,4)]

write.table(ind, "~/programs/gdc/ph_filtered_snps_ind2", sep="\t", row.names = F)


library(adegenet)
library(vcfR)

ph<-read.vcfR("Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf")
save(ph, file="Output/ph_vcfR.RData")

ph_genind <- vcfR2genind(ph)
save(ph_genind, file="ph_genind.RData")




#phvcf<-read.vcf("Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf",from = 1, to =1e6)

#library(reconproGS)
#ph <- vcf2genind("Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf")
save(ph, file="ph_genind.RData")



#Find clusters
phgroups<-find.clusters(ph_genind, max.n.clust = 40)
table(phgroups$grp)
phgroups$grp[phgroups$grp==2] #CA
phgroups$grp[phgroups$grp==3] #TB
phgroups$grp[phgroups$grp==4]

#dapc
dapc1<-dapc(ph_genind,phgroups$grp )

#plot the result
scatter(dapc1)

scatter(dapc1,1,1, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)


library(pcadapt)
