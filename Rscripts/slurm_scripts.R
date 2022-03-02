# create slurm scripts to run evalAdmix
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)
    
    
qfiles<-list.files("Data/ngsadmix/qfiles/")
ffiles<-list.files("Data/ngsadmix/", pattern =".fopt.gz")

sink("Output/Slurmscripts/calculateAF_eachPop.sh")
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=calculateAF
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error calculateAF.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("modeul load vcftools\n\n")
for (i in 1: length(pops)){
    cat(paste0("vcftools --gzvcf /home/jamcgirr/ph/data/vcfs/split_pops/maf05/population_",pops[i], "_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz "))
    cat(paste0("--freq --out /home/ktist/ph/data/AF/",pops[i],"_freq \n"))
}
sink(NULL)
#

#index vcf maf00 files for bcf tools 
sink("Output/Slurmscripts/bcftools_index.sh")
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=bcftools_index
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error bcftools_index.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load bcftools\n\n")
for (i in 1: length(pops)){
    cat(paste0("bcftools index /home/ktist/ph/data/vcf/maf00/population_",pops[i], "_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5.vcf.gz \n"))
}
sink(NULL)



## Obtain population genetic statistiscs for each pop group
windws<-read.table("Data/vcfs/windows.txt", sep="\t")

#Too big to fit in 1 script. Divide them up to each pop
for (i in 1: length(pops)){
    sink(paste0("Output/Slurmscripts/het_perwindow_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat("#SBATCH --job-name=het_perwindow
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error het_perwindow.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load bcftools\n\n")
    cat(paste0("mkdir /home/ktist/ph/data/vcf/maf00/stats/", pops[i],"/ \n"))
    for (j in 1: nrow(windws)){
        reg<-paste0(windws$V1[j],":",windws$V2[j],"-",windws$V3[j])
        cat(paste0("bcftools stats -r ", reg," -s - /home/ktist/ph/data/vcf/maf00/population_",pops[i], "_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5.vcf.gz "))
        cat(paste0("| grep '^PSC' > /home/ktist/ph/data/vcf/maf00/stats/",pops[i],"/", pops[i], "_stats_",j,"\n"))        
    }
    sink(NULL)
}
#

#add file name and concatinate all files

for (i in 1: length(pops)){
    sink(paste0("Output/Slurmscripts/catFiles_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat("#SBATCH --job-name=catFiles_
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error het_perwindow.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
    cat("\n")
    cat(paste0("cd ",pops[i],"/\n"))
    cat(paste0("for f in *; do sed -i "))
    cat('"s/$/\\t$f/" $f; done \n')
    cat(paste0("cat $(ls -t) > ","/home/ktist/ph/data/vcf/maf00/stats/", pops[i], "_statsFile\n"))
    cat("cd ..\n")
    sink(NULL)
}

sink(paste0("Output/Slurmscripts/pop_stats_.sh"))
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=pop_stats_
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error het_perwindow.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
cat("\n")
for (i in 1: length(pops)){

    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/vcf/maf00/population_",pops[i], "_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5.vcf.gz --window-pi 100000 --out /home/ktist/ph/data/vcf/maf00/pi/",pops[i],"\n"))
}
sink(NULL)




sink(paste0("Output/Slurmscripts/ROH_plink.sh"))
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=ROH_plink
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error ROH_plink.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
cat("\n")
for (i in 1: length(pops)){
    cat(paste0("plink -vcf /home/ktist/ph/data/vcf/maf00/population_",pops[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5.vcf.gz --const-fid ",pops[i]," --make-bed --out /home/ktist/ph/data/vcf/maf00/",pops[i],"_maf00\n"))
}
sink(NULL)


sink(paste0("Output/Slurmscripts/ROH_plink3.sh"))
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=ROH_plink2
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error ROH_plink2.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
cat("\n")
cat("module load plink\n\n")
for (i in 1: length(pops)){
    cat(paste0("plink -bfile /home/ktist/ph/data/vcf/maf00/",pops[i],"_maf00 --homozyg group --homozyg-snp 100 --homozyg-gap 500 --homozyg-window-missing 15 --homozyg-density 75 --out /home/ktist/ph/data/vcf/maf00/",pops[i],"_miss15\n"))
} 

sink(NULL)


sink(paste0("Output/Slurmscripts/makePed.sh"))
#cat("#!/bin/bash -l\n\n")
#cat("#SBATCH --job-name=ROH_plink2
##SBATCH --mem=12G
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=12
##SBATCH --error ROH_plink2.err
##SBATCH --time=24:00:00
##SBATCH -p high \n")
#cat("\n")
#cat("module load plink\n\n")
for (i in 1: length(pops)){
    cat(paste0("plink --bfile /home/ktist/ph/data/vcf/maf00/",pops[i],"_maf00 --recode --tab --out /home/ktist/ph/data/vcf/maf00/",pops[i],"_maf00 \n"))
} 

sink(NULL)



sink(paste0("Output/Slurmscripts/zip_Ped.sh"))
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=zip_Ped
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error zip_Ped.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
cat("\n")
for (i in 1: length(pops)){
    cat(paste0("gzip /home/ktist/ph/data/vcf/maf00/",pops[i] ,"_maf00.ped \n"))
}
sink(NULL)


## MDS for minimum genotype freq 70% data 58,432 SNPs

sink("Output/Slurmscripts/plink_MDS70.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=plink_MDS70 \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e plinkMDA_byChromosome.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load plink")     
cat("\n\n")



cat(paste0("plink --vcf /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70.vcf --make-bed --out /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70.vcf \n")) 
cat(paste0("plink --bfile /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70 --recode --tab --out /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70 \n"))

cat(paste0("#Run MDS and PCA \n"))
cat(paste0("plink --file /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70 --mds-plot 10 --cluster --out /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70_mds_10 \n"))
cat(paste0("plink --file /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70 --pca --out /home/ktist/ph/data/vcf/maf05ID/ph_maf05_id_ov70_pca \n"))
cat("\n")
sink(NULL)


# Extract a region from bam file

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws07<-pop_info[pop_info$Population.Year=="PWS07",]


sink("Output/Slurmscripts/extract_chr7_pws07.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_chr7 \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_chr7.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load samtools")     
cat("\n\n")

for (i in 1: nrow(pws07)){
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/',pws07$Sample[i] , '.bam "chr7:18000000-28000000" >  /home/ktist/ph/data/bam/', pws07$Sample[i], '_chr7_subregion.bam \n'))
}
sink(NULL)


### extract_inversion_chr12

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
ca<-pop_info[pop_info$Population.Year=="CA17",]


sink("Output/Slurmscripts/extract_inversion.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_inversion \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_inversion.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load samtools")     
cat("\n\n")

for (i in 1: nrow(ca)){
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/', ca$Sample[i],'.bam "chr12:17748000-17754000" >  /home/ktist/ph/data/bam/', ca$Sample[i], '_chr12_breakpoint1.bam \n')) 
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/', ca$Sample[i],'.bam "chr12:25760000-25780000" >  /home/ktist/ph/data/bam/', ca$Sample[i], '_chr12_breakpoint2.bam \n')) 
}
sink(NULL)

