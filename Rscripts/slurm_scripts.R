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

#remove the head line and cat beagle files
sink("Output/Slurmscripts/cat_beagle_files")
for (i in 2:26){
    cat(paste0("sed -e '1, 1d' < PWS_SS_maf05_c",i,".BEAGLE.PL > PWS_SS_maf05_c",i,".BEAGLE.PL \n"))
}
cat("cat ") 
for (i in 1:26){
    cat(paste0("PWS_SS_maf05_c",i,".BEAGLE.PL "))
}
cat(paste0(" > PWS_SS_maf05_BEAGLE.PL \n"))
sink(NULL)


sink("Output/Slurmscripts/create_beagle_files")
for (i in 1:26){
    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/split_vcf/ph_noTB_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz  --chr chr",i))
    cat(paste0(" --out /home/ktist/ph/data/split_vcf/beagle/noTB/noTB_maf05_c",i," --BEAGLE-PL \n"))
}
sink(NULL)



cat("cp -R /home/ktist/ph/data/split_vcf/beagle/noTB/ /home/ktist/ph/data/split_vcf/beagle/noTB/noTB_backup")

sink("Output/Slurmscripts/runPCAngsd.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=pca \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e pca.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load pcaangsd")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("gzip /home/ktist/ph/data/split_vcf/beagle/noTB/noTB_maf05_c",i,".BEAGLE.PL \n"))
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/split_vcf/beagle/noTB/noTB_maf05_c",i,".BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/PCAngsd/noTB_maf05_chr",i," -threads 16 \n")) 
}    
    
sink(NULL)
    

#subset by chromosomes

sink("Output/Slurmscripts/subset_byChrom3.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=subset_byChrom3 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e subset_byChrom3.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load bcftools")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("bcftools view /home/ktist/ph/data/new_vcf/MD7000/PH_minDP600_maxDP7000_minQ20_minMQ30_NS0.5.vcf.gz -r chr",i," -Oz --threads 12 > /home/ktist/ph/data/new_vcf/MD7000/PH_MD7000_chr",i,".vcf.gz \n"))

}    

sink(NULL)

# Extract a region from chr6 bam file

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws07<-pop_info[pop_info$Population.Year=="PWS07",]


sink("Output/Slurmscripts/extract_chr6_pws07.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_chr6 \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_chr6.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load samtools")     
cat("\n\n")

for (i in 1: nrow(pws07)){
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/',pws07$Sample[i] , '.bam "chr6:22000000-32000000" >  /home/ktist/ph/data/bam/', pws07$Sample[i], '_chr6_subregion.bam \n'))
}
sink(NULL)


pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)
for (i in 1: length(pops)){
    sink(paste0("Data/vcfs/Population/",pops[i],".txt"))
    samples<-pop_info$Sample[pop_info$Population.Year==pops[i]]
    snames<-sapply(samples, paste, collapse = '\n')
    cat(snames, sep="\n")
    sink(NULL)
}    

sink("Output/Slurmscripts/subset_vcf.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=subset_vcf \n"))
cat(paste0("#SBATCH --mem=12G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e esubset_vcf.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load bcftools")     
cat("\n\n")

for (i in 1: length(pops)){
    cat(paste0('bcftools view -Oz -S /home/ktist/ph/data/new_vcf/MD7000/population/',pops[i],'.txt --threads 12  /home/ktist/ph/data/new_vcf/MD7000/PH_minDP600_maxDP7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz > /home/ktist/ph/data/new_vcf/MD7000/population/',pops[i],'_maf05.vcf.gz \n'))
}

sink(NULL)

## extract allele count data for each population            
sink("Output/Slurmscripts/extract_AC.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_AC \n"))
cat(paste0("#SBATCH --mem=12G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_AC.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load bcftools")     
cat("\n\n")

for (i in 1: length(pops)){
    cat(paste0("bcftools query -f '%INFO/AC  %INFO/AN\n' /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"_maf05.vcf.gz > /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"_AC.txt \n"))
}
sink(NULL)

#Extract allele count dta from vcf files
bcftools query -f '%INFO/AC  %INFO/AN\n ' mini.vcf > 
    

## extract read depth data for each population            
sink("Data/Slurmscripts/extract_coverage.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_coverage \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e extract_coverage.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load bcftools")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("bcftools query -f '"))
    cat(paste0("%CHROM  %POS  %INFO/DP\\n' "))
    cat(paste0("/home/jamcgirr/ph/data/combine_gvcfs/raw_variants_chr",i,"_1.vcf > /home/ktist/ph/data/new_vcf/raw_DP/chr",i,".1.info \n"))
    cat(paste0("bcftools query -f '"))
    cat(paste0("%CHROM  %POS  %INFO/DP\\n' "))
    cat(paste0("/home/jamcgirr/ph/data/combine_gvcfs/raw_variants_chr",i,"_2.vcf > /home/ktist/ph/data/new_vcf/raw_DP/chr",i,".2.info \n"))
    
}
sink(NULL)


sink("Data/Slurmscripts/zip.sh")
cat("zip /home/ktist/ph/data/new_vcf/raw_DP/dp2.zip ")
for (i in 11:20){
    cat(paste0("/home/ktist/ph/data/new_vcf/raw_DP/chr",i,".1.info /home/ktist/ph/data/new_vcf/raw_DP/chr",i,".2.info "))
}
cat("\n")
cat("zip /home/ktist/ph/data/new_vcf/raw_DP/dp3.zip ")

for (i in 21:26){
    cat(paste0("/home/ktist/ph/data/new_vcf/raw_DP/chr",i,".1.info /home/ktist/ph/data/new_vcf/raw_DP/chr",i,".2.info "))
}

cat("\n")    
sink(NULL)    
    


### Convert pruned vcf to beagle format

sink("Data/Slurmscripts/convert_beagle.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=convert_beagle \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e convert_beagle.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load vcftools")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/split_vcf/population/pruned2017.vcf.gz  --chr chr",i))
    cat(paste0(" --out /home/ktist/ph/data/split_vcf/population/beagle/pruned2017_c",i," --BEAGLE-PL \n"))
}

cat("\n")
#remove the head line and cat beagle files
for (i in 2:26){
    cat(paste0("sed -e '1, 1d' < /home/ktist/ph/data/split_vcf/population/beagle/pruned2017_c",i,".BEAGLE.PL > /home/ktist/ph/data/split_vcf/population/beagle/pruned2017_c",i,".2.BEAGLE.PL \n"))
}
cat("\n")

cat("cat ") 
for (i in 1:26){
    cat(paste0("/home/ktist/ph/data/split_vcf/population/beagle/pruned2017_c",i,".2.BEAGLE.PL "))
}
cat(paste0(" > /home/ktist/ph/data/split_vcf/population/beagle/pruned2017_BEAGLE.PL \n"))
sink(NULL)

## run selection scan using pcangsd

pops<-c("pws","ss","tb")

sink("Data/Slurmscripts/selection.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=selection \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e selection.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load bcftools")     
cat("\n\n")

for (i in 1:3){
    cat(paste0("bcftools index /home/ktist/ph/data/new_vcf/MD7000/population/ph_",pops[i], "_MD7000_maf05.vcf.gz \n" ))
    cat(paste0("bcftools view -R /home/ktist/ph/data/new_vcf/MD7000/population/in75MD7000.txt /home/ktist/ph/data/new_vcf/MD7000/population/ph_",pops[i], "_MD7000_maf05.vcf.gz > /home/ktist/ph/data/new_vcf/MD7000/population/pruned",pops[i], ".vcf \n"))
    cat(paste0("bgzip -c /home/ktist/ph/data/new_vcf/MD7000/population/pruned",pops[i], ".vcf  > /home/ktist/ph/data/new_vcf/MD7000/population/pruned",pops[i], ".vcf.gz \n"))

}
    
cat("\n\n")
cat("module load vcftools \n\n")

for (j in 1:3){
    for (i in 1:26){
        cat(paste0("vcftools --gzvcf /home/ktist/ph/data/new_vcf/MD7000/population/pruned",pops[j],".vcf.gz  --chr chr",i))
        cat(paste0(" --out /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_c",i," --BEAGLE-PL \n"))
    }
    
    
}

sink(NULL)


sink("Data/Slurmscripts/selection2.sh")
cat("#!/bin/bash -l \n")
cat(paste0("#SBATCH --job-name=selection2 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e selection2.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load bcftools \n\n")     

for (j in 1:3){
    for (i in 2:26){
        cat(paste0("sed -e '1, 1d' < /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_c",i,".2.BEAGLE.PL \n"))
    }}
cat("\n\n") 
for (j in 1:3){
    cat("cat ")
    cat(paste0("/home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_c1.BEAGLE.PL "))
    
    for (i in 2:26){
        cat(paste0("/home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_c",i,".2.BEAGLE.PL "))
    }
    cat(paste0(" > /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_BEAGLE.PL \n"))
}
cat("\n\n")

for (j in 1:3){
    cat(paste0("gzip /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_BEAGLE.PL \n"))
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned",pops[j],"_BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/selection/pruned.",pops[j],"_selection -selection -sites_save #-n $N \n\n"))
}
sink(NULL)


sink("Data/Slurmscripts/selection2.sh")
cat("#!/bin/bash -l \n")
cat(paste0("#SBATCH --job-name=selection2 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e selection2.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load bcftools \n\n")     

for (j in 1:892){
    f}
cat("\n\n") 



####### Create beagle geno files for runnign ngsLD (4.11.22)
sink("Data/Slurmscripts/convert_beagle2.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=convert_beagle2 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e convert_beagle2.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load vcftools")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/new_vcf/MD7000/PH_minDP600_maxDP7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz  --chr chr",i))
    cat(paste0(" --out /home/ktist/ph/data/new_vcf/MD7000/beagle/MD7000_maf05_c",i," --BEAGLE-PL \n"))
}
sink(NULL)

sink("Data/Slurmscripts/remove_headers.sh")

#remove the head line and cat beagle files
for (i in 1:26){
    cat(paste0("sed -e '1d' < /home/ktist/ph/data/new_vcf/MD7000/beagle/MD7000_maf05_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".BEAGLE.PL \n"))
}

#remove the first 3 columns
for (i in 1:26){
    cat(paste0("cut -f 4- /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".PL \n"))
    cat(paste0("bgzip /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".PL \n"))
}

sink(NULL)


## subdivide PL file into each pop.year
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)

sink("Data/Slurmscripts/prep_PLfiles.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=prep_PLfiles \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e prep_PLfiles.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")

for (j in 1:14){
    pop<-pops[j]
    for (i in 1:26){
        cat(paste0("awk -f /home/ktist/ph/tst.awk ph/data/new_vcf/MD7000/population/", pop,".txt /home/ktist/ph/data/new_vcf/MD7000/beagle/MD7000_maf05_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/MD7000_maf05_c",i,".",pop," \n"))
            }
}


sink(NULL)


sink("Data/Slurmscripts/zip_PLfiles.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=zip_PLfiles \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e zip_PLfiles.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")



sink("Data/Slurmscripts/remove_headers.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=remove_headers \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e remove_headers.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")


#remove the first 2 columns and zip beagle files
for (j in 1:14){
    pop<-pops[j]
    for (i in 1:26){
        cat(paste0("cut -f 3- /home/ktist/ph/data/new_vcf/MD7000/beagle/PopYr/MD7000_maf05_c",i,".", pop," > /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_c",i,".",pop," \n"))
        cat(paste0("sed -e '1d' < /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_c",i,".",pop," > /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".",pop,"  \n"))
        cat(paste0("bgzip /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".",pop," \n"))
    }
}

    
sink(NULL)



sink("Data/Slurmscripts/zip_pl.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=zip_pl \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e zip_pl.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
for (j in 1:14){
    pop<-pops[j]
    for (i in 1:26){
                cat(paste0("bgzip /home/ktist/ph/data/new_vcf/MD7000/beagle/noHeader/MD7000_maf05_noHeader_c",i,".",pop," \n"))
    }
}

sink(NULL)



### Run ngsLD

pos<-read.csv("Data/new_vcf/new_positions", sep="\t", header = F)
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)


sink("Data/Slurmscripts/run_ngsLD.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=run_ngsLD \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e run_ngsLD.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")

for (j in 1:14){
    pop<-pops[j]
    s<-nrow(pops_info[pops_info$Population.Year==pop,])
    for (i in 1:26){
        n<-nrow(pos[pos$V1==paste0("chr",i),])
        cat(paste0("./programs/ngsLD/ngsLD --geno /home/ktist/ph/data/new_vcf/MD7000/beagle/Pop.Yr/MD7000_maf05_c",i,"_", pop,".gPL.gz --max_kb_dist 0 --n_ind ",s," --n_sites ", n, " --prob --pos /home/ktist/ph/data/new_vcf/MD7000/beagle/positions/pos_chr",i,".txt --n_threads 12 --out /home/ktist/ph/data/ngsLD/ph_chr", i,".",pop, ".ld --extend_out \n"))
    }
    
}

sink(NULL)

sink("Data/Slurmscripts/pcangsd_selection.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=pcangsd_selection \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e pcangsd_selection.err  \n"))
cat(paste0("#SBATCH --time=48:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
for (i in 1:6){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/population/beagle/pruned_",pop1,"_",pop2,".BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/selection/pruned_",pop1,"_",pop2,"_selection -selection -sites_save \n"))
}
sink(NULL)

####### 4.18.22
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

sink("Data/Slurmscripts/calculateAF_eachPop2.sh")
for (i in 1: length(pops)){
    cat(paste0("vcftools --gzvcf Data/new_vcf/population/",pops[i], "_maf05.vcf.gz "))
    cat(paste0("--freq --out Data/new_vcf/population/",pops[i],"_freq \n"))
}
sink(NULL)
#






#####
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws96<-pops_info$Sample[pops_info$Population.Year=="PWS96"]
pws91<-pops_info$Sample[pops_info$Population.Year=="PWS91"]
pws07<-pops_info$Sample[pops_info$Population.Year=="PWS07"]
ca<-pops_info$Sample[pops_info$Population.Year=="CA17"]

sink("Data/Slurmscripts/Extract_Depth1.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=Extract_Depth1 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e Extract_Depth1.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load samtools \n") 

for (i in 1:10){
    cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/chr6_1cyp.bed /home/eoziolor/phpopg/data/align/", pws96[i],".bam > /home/ktist/ph/data/bam_depth/",pws96[i],"_chr6-cyp.txt \n"))
    cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws96[i],"_chr6-cyp.txt \n"))
    cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/chr6_1cyp.bed /home/eoziolor/phpopg/data/align/", pws91[i],".bam > /home/ktist/ph/data/bam_depth/",pws91[i],"_chr6-cyp.txt \n"))
    cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws91[i],"_chr6-cyp.txt \n"))
    cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/chr6_1cyp.bed /home/eoziolor/phpopg/data/align/", pws07[i],".bam > /home/ktist/ph/data/bam_depth/",pws07[i],"_chr6-cyp.txt \n"))
    cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws07[i],"_chr6-cyp.txt \n"))
    
}
sink(NULL)

sink("Data/Slurmscripts/Extract_Depth1_ca.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=Extract_Depth1 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e Extract_Depth1.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load samtools \n") 

for (i in 1:10){
    cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/chr6_1cyp.bed /home/eoziolor/phpopg/data/align/", ca[i],".bam > /home/ktist/ph/data/bam_depth/",ca[i],"_chr6-cyp.txt \n"))
    cat(paste0("gzip /home/ktist/ph/data/bam_depth/",ca[i],"_chr6-cyp.txt \n"))
}
sink(NULL)

#####  for differente genes in the ahR pathway
beds<-list.files("Output/CNV/ahr/bed")
for (j in 1: length(beds)) {
    gene<-gsub(".bed","",beds[j])
    sink(paste0("Data/Slurmscripts/extract_depth_",gene,".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=Extract_depth2 \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --ntasks=1 \n")) 
    cat(paste0("#SBATCH -e Extract_depth2.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    #cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    #cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load samtools \n\n") 

    for (i in 1:10){
        cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/", beds[j]," /home/eoziolor/phpopg/data/align/", pws96[i],".bam > /home/ktist/ph/data/bam_depth/",pws96[i],"_", gene,".txt \n"))
        cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws96[i],"_", gene,".txt  \n"))
        cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/", beds[j]," /home/eoziolor/phpopg/data/align/", pws91[i],".bam > /home/ktist/ph/data/bam_depth/",pws91[i],"_", gene,".txt \n"))
        cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws91[i],"_", gene,".txt  \n"))
        cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/", beds[j]," /home/eoziolor/phpopg/data/align/", pws07[i],".bam > /home/ktist/ph/data/bam_depth/",pws07[i],"_", gene,".txt \n"))
        cat(paste0("gzip /home/ktist/ph/data/bam_depth/",pws07[i],"_", gene,".txt  \n"))
        
        cat(paste0("samtools depth -b /home/ktist/ph/data/bam_depth/bed/", beds[j]," /home/eoziolor/phpopg/data/align/", ca[i],".bam > /home/ktist/ph/data/bam_depth/",ca[i],"_", gene,".txt \n"))
        cat(paste0("gzip /home/ktist/ph/data/bam_depth/",ca[i],"_", gene,".txt  \n"))
    }
    
    sink(NULL)
    
}

##CAlculate allele frequency for each pop from MD7000 vcf file.
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

sink("Data/Slurmscripts/calculateAF_eachPop.sh")
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=calculateAF
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error calculateAF.err
#SBATCH --time=48:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load vcftools\n\n")
for (i in 1: length(pops)){
    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i], "_maf05.vcf.gz "))
    cat(paste0("--freq --out /home/ktist/ph/data/AF/MD7000/",pops[i],"_md7000_maf05_freq \n"))
}
sink(NULL)



##CAlculate allele frequency using ANGSD for each pop from MD7000 vcf file.
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

sink("Data/Slurmscripts/calculateAF_angsd.sh")
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=calculateAF
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error calculateAF.err
#SBATCH --time=48:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load angsd\n\n")
for (i in 1: length(pops)){
    cat(paste0("angsd -out /home/ktist/ph/data/new_vcf/MD7000/AF/",pops[i]," -fai /home/jamcgirr/ph/data/c_harengus/c.harengus.fa.fai -doGlf 2 -doMaf 3 -doMajorMinor 4 -doPost 1 -doGeno 2 -vcf-pl /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i],"_maf05.vcf.gz -ref /home/jamcgirr/ph/data/c_harengus/c.harengus.fa \n"))
}
sink(NULL)


### calculate Fst using new datasets 

#first create folded SFS, and then calculate Fst in ANGSD
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)


pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","TB17")
sss<-c("SS96","SS06","SS17")
y17<-pops[grep("17",pops)]

#With VCF maf00
for (i in 1:length(pops)){
    
    sink(paste0("Data/Slurmscripts/",pops[i],"_angsd_SFS.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_SFS \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_SFS.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n") 

    cat(paste0("angsd -doSaf 1 -vcf-pl /home/ktist/ph/data/new_vcf/MD7000/",pops[i],"_filtered_snps.vcf.gz -out /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00 -anc /home/jamcgirr/ph/data/c_harengus/c.harengus.fa\n")) 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.saf.idx -P 16 > /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.unfolded.sfs \n")) 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS saf2theta /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.saf.idx -sfs /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.unfolded.sfs -outname /home/ktist/ph/data/angsd/SFS/fromVCF/", pops[i],"_maf00\n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromVCF/", pops[i],"_maf00.thetas.idx -win 50000 -step 10000 -outnames /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.thetas100kWindow.gz \n\n"))
    #folded SFS
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.saf.idx -P 16 -fold 1 > /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.folded.sfs \n")) 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS saf2theta /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.saf.idx -sfs /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.folded.sfs -outname /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pops[i],"_maf00\n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pops[i],"_maf00.thetas.idx -win 50000 -step 10000 -outnames /home/ktist/ph/data/angsd/SFS/fromVCF/folded_",pops[i],"_maf00.thetas50kWindow.gz \n"))
    
    sink(NULL)
}

#with VCF maf05
for (i in 1:length(pops)){
    
    sink(paste0("Data/Slurmscripts/",pops[i],"_angsd_SFS_maf05.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_SFS_maf05 \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_SFS_maf05.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n") 
    
    cat(paste0("angsd -doSaf 1 -vcf-pl /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i],"_maf05.vcf.gz -out /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05 -anc /home/jamcgirr/ph/data/c_harengus/c.harengus.fa\n")) 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05.saf.idx -P 16 > /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05.unfolded.sfs \n")) 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS saf2theta /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05.saf.idx -sfs /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05.unfolded.sfs -outname /home/ktist/ph/data/angsd/SFS/fromVCF/", pops[i],"_maf05\n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromVCF/", pops[i],"_maf05.thetas.idx -win 100000 -step 20000 -outnames /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf05.thetas100kWindow.gz \n")) 
    sink(NULL)
}




## Calculate unfolded SFS
sink("Data/Slurmscripts/unfolded_SFS.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=unfolded_SFS \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e unfolded_SFS.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n") 

for (i in 1:14){
    #create 
    cat(paste0("realSFS /home/ktist/ph/data/angsd/SFS/",pops[i],".saf.idx -P 16 > /home/ktist/ph/data/angsd/SFS/",pops[i],"_unfolded.sfs \n")) 
 }
sink(NULL)

#Make 2D unfolded SFS for population pairs
pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","TB17")
sss<-c("SS96","SS06","SS17")
y17<-pops[grep("17",pops)]

comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-t(comb4)

sink("Data/Slurmscripts/2DSFS_2.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=2DSFS \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e 2DSFS.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n\n") 


for (i in 1:nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/", pop1, ".saf.idx /home/ktist/ph/data/angsd/SFS/", pop2, ".saf.idx -P 16 > /home/ktist/ph/data/angsd/SFS/unfolded_", pop1,"_",pop2,".sfs \n")) 
}    
cat("\n")
for (i in 1:nrow(comb2)){
    
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    sink(paste0("Data/Slurmscripts/",pop1,pop2,"2DSFS.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pop1,pop2,"2DSFS \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pop1,pop2,"2DSFS.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n\n") 
    
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2, "_maf00.saf.idx -fold 1  -P 16 > /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs \n")) 
} 

#SS
for (i in 1:nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    sink(paste0("Data/Slurmscripts/",pop1,pop2,"2DSFS.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pop1,pop2,"2DSFS \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pop1,pop2,"2DSFS.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n\n") 
    
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2, "_maf00.saf.idx -fold 1  -P 16 > /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs \n")) 
    sink(NULL)
} 
   

for (i in 1:nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    sink(paste0("Data/Slurmscripts/",pop1,pop2,"2DSFS.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pop1,pop2,"2DSFS \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pop1,pop2,"2DSFS.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n\n") 
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2, "_maf00.saf.idx -fold 1  -P 16 > /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs \n")) 
    sink(NULL)
} 



## Fst calculation from VCF maf00
for (i in 1:nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    
    sink(paste0("Data/Slurmscripts/FstPWS_maf00_",i,".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=FstPWS",i," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e FstPWF",i,".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n") 
    
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2,"_maf00.saf.idx "))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs -fstout /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00 \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_50kWindow_maf00 \n"))
    sink(NULL)
}
    

for (i in 1:nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    sink(paste0("Data/Slurmscripts/FstTB_maf00_",i,".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=FstTB",i," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e FstTB",i,".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2,"_maf00.saf.idx "))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs -fstout /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00 \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_50kWindow_maf00 \n"))
    sink(NULL)
} 
for (i in 1:nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    sink(paste0("Data/Slurmscripts/FstSS_maf00_",i,".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=FstSS",i," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e FstSS",i,".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2,"_maf00.saf.idx "))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs -fstout /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00 \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_50kWindow_maf00 \n"))
    sink(NULL)} 
cat("\n")
for (i in 1:nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    sink(paste0("Data/Slurmscripts/Fst17_",i,".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=Fst17",i," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e Fst17",i,".err  \n"))
    cat(paste0("#SBATCH --time=144:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2,"_maf00.saf.idx "))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs -fstout /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00 \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_persite_maf00.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fromVCF/fst_",pop1, "_",pop2,"_50kWindow_maf00 \n"))
    sink(NULL)
    } 
cat("\n")

sink(NULL)


# 3D SFS & Fst fir 2017 files
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)
y17<-pops[grep("17",pops)]
comb<-combn(y17, 3)
comb<-t(comb)
#comb<-comb[!duplicated(apply(comb, 1, sort), MARGIN = 1), ]


sink("Data/Slurmscripts/Fst_3D.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=Fst_3D \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e Fst_3D.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n") 


for (i in 1:nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    pop3<-comb[i,3]
    
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/", pop1, ".saf.idx /home/ktist/ph/data/angsd/SFS/", pop2,".saf.idx  /home/ktist/ph/data/angsd/SFS/", pop3,".saf.idx " ))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/unfolded_", pop1,"_",pop2,".sfs -sfs /home/ktist/ph/data/angsd/SFS/unfolded_", pop1,"_",pop3,".sfs -sfs /home/ktist/ph/data/angsd/SFS/unfolded_", pop2,"_",pop3,".sfs -fstout /home/ktist/ph/data/angsd/SFS/fst_",pop1, "_",pop2,"_",pop3,"_persite \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fst_",pop1, "_",pop2,"_",pop3,"_persite.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fst_",pop1, "_",pop2,"_",pop3,"_50kWindow \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst print /home/ktist/ph/data/angsd/SFS/fst_unfolded_",pop1, "_",pop2,"_",pop3,".fst.idx > /home/ktist/ph/data/angsd/SFS/fst_unfolded_",pop1, "_",pop2,"_",pop3,".fst.txt \n\n"))
    
}

sink(NULL)


pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

for (i in 1: length(pops)){
    df<-pops_info[pops_info$Population.Year==pops[i],]
    df$bam<-paste0("/home/eoziolor/phpopg/data/align/",df$Sample,".bam")
    write.table(df[,"bam"], paste0("Data/bam/samplenames/", pops[i],".txt"), quote=F, row.names = F, col.names = F)
    
}



sink("Data/Slurmscripts/PWS07_filenames.sh")
for (i in 1:26){
    cat(paste0("/home/ktist/ph/data/angsd/SFS/fromBam/PWS07_unfolded_chr",i,".sfs.idx "))
}
cat(paste0("-outnames /home/ktist/ph/data/angsd/SFS/fromBam/PWS07_unfolded_combined "))
cat("\n\n")
sink(NULL)

sink("Data/Slurmscripts/Combine_sfs.sh")
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=Combine_sfs \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e Combine_sfs.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n") 

cat("realSFS cat ")
for (i in 1:26){
    cat(paste0("/home/ktist/ph/data/angsd/SFS/fromBam/PWS07_unfolded_chr",i,".sfs "))
}
cat(paste0("-outnames /home/ktist/ph/data/angsd/SFS/fromBam/PWS07_unfolded_combined "))
cat("\n\n")
sink(NULL)


## run sfs from bam files. Create .sfs file by dividing the genome into each chromosome
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)

for (i in 1:length(pops)){
    sink(paste0("Data/Slurmscripts/",pops[i],"_sfs_step1.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_sfs_step1\n"))
    cat(paste0("#SBATCH --mem=16G \n" )) 
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_sfs_step2.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n")
    
    cat(paste0("angsd -bam /home/ktist/ph/data/angsd/samples/",pops[i],".txt -doSaf 1 -doMajorMinor 1 -GL 2 -doMaf 3 -doCounts 1 -doGlf 3  -anc /home/jamcgirr/ph/data/c_harengus/c.harengus.fa -ref /home/jamcgirr/ph/data/c_harengus/c.harengus.fa -minMapQ 30 -minQ 20 -P 8 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minInd 10 -setMinDepth 10 -setMaxDepth 500 -out /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],"\n\n"))
    sink(NULL)
    
    
    sink(paste0("Data/Slurmscripts/",pops[i],"_sfs_step2.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_sfs_step2\n"))
    cat(paste0("#SBATCH --mem=16G \n" )) 
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_sfs_step2.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n")
    
    for (j in 1:26){
        cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromBam/",pops[i],".saf.idx -P 24 -r chr",j," > /home/ktist/ph/data/angsd/SFS/fromBam/",pops[i],"_unfolded_chr",j,".sfs \n")) 
    }
    sink(NULL)
    
}


#create folded sfs as well 
for (i in 1:length(pops)){
    sink(paste0("Data/Slurmscripts/",pops[i],"_sfs_step2_folded.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_sfs_step2_folded\n"))
    cat(paste0("#SBATCH --mem=16G \n" )) 
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_sfs_step2_folded.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n")
    
    for (j in 1:26){
        cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromBam/",pops[i],".saf.idx -P 24 -r chr",j," -fold 1 > /home/ktist/ph/data/angsd/SFS/fromBam/folded/",pops[i],"_folded_chr",j,".sfs \n")) 
    }
    sink(NULL)
    
}
for (i in 1:length(pops)){
    sink(paste0("Data/Slurmscripts/",pops[i],"_sfs_theta.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_sfs_theta\n"))
    cat(paste0("#SBATCH --mem=16G \n" )) 
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_sfs_theta.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n")
    
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/realSFS saf2theta /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],".saf.idx -sfs /home/ktist/ph/data/angsd/SFS/fromBam/",pops[i],"_unfolded.sfs -outname /home/ktist/ph/data/angsd/SFS/fromBam/",pops[i],"\n"))
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],".thetas.idx \n"))
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],".thetas.idx -win 50000 -step 10000 -outnames /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],"_50kwin_10kstep \n\n"))

    sink(NULL)
}

for (i in 1:length(pops)){
    sink(paste0("Data/Slurmscripts/",pops[i],"_sfs_step3_folded.sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"_sfs_step3_folded\n"))
    cat(paste0("#SBATCH --mem=16G \n" )) 
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"_sfs_step3_folded.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n")
    
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/realSFS saf2theta -fold 1 /home/ktist/ph/data/angsd/SFS/fromBam/", pops[i],".saf.idx -sfs /home/ktist/ph/data/angsd/SFS/fromBam/folded/",pops[i],"_folded.sfs -outname /home/ktist/ph/data/angsd/SFS/fromBam/folded/",pops[i],"\n"))
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromBam/folded/", pops[i],".thetas.idx \n"))
    cat(paste0("/home/jamcgirr/apps/angsd_sep_20/angsd/misc/thetaStat do_stat /home/ktist/ph/data/angsd/SFS/fromBam/folded/", pops[i],".thetas.idx -win 50000 -step 10000 -outnames /home/ktist/ph/data/angsd/SFS/fromBam/folded/", pops[i],"_50kwin_10kstep \n\n"))
    
    sink(NULL)
}



#Find minimum genotyping rate of 0.5

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops.info$Population.Year)

for (i in 1:length(pops)){
    no<-nrow(pops.info[pops.info$Population.Year==pops[i],])/2
    
    sink(paste0("Data/Slurmscripts/Subset_byGenotype_rate_",pops[i],".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=",pops[i],"__subset_genotype_rate\n"))
    cat(paste0("#SBATCH --nodes=4 \n" )) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH -e ",pops[i],"__subset_genotype_rate.err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load bcftools  \n")
    
    cat(paste0("bcftools view -S /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i],".txt /home/ktist/ph/data/new_vcf/MD7000/merged_filtered_snps.bcf | bcftools +fill-tags -- -t all,'DP=sum(DP)' | bcftools filter -Oz -i 'NS>",no, "' > /home/ktist/ph/data/new_vcf/MD7000/",pops[i],"_filtered_snps.vcf.gz \n"))
    cat(paste0("bcftools index /home/ktist/ph/data/new_vcf/MD7000/",pops[i],"_filtered_snps.vcf.gz \n"))
    sink(NULL)
}    


### Run PCAngsd for pruned new vcf file (NS0.5, MD7000, maf05)

sink("Data/Slurmscripts/runPCAngsd.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=pca \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --nodes=4 \n" )) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
cat(paste0("#SBATCH -e pca.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load pcaangsd")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("gzip /home/ktist/ph/data/new_vcf/MD7000/prune/beagle/MD7000_maf05_c",i,".BEAGLE.PL \n"))
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/prune/beagle/MD7000_maf05_c",i,".BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/PCAngsd/MD7000_maf05_chr",i," -threads 16 \n")) 
}    

sink(NULL)

sink("Data/Slurmscripts/runPCAngsd2.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=pca2 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --nodes=4 \n" )) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
cat(paste0("#SBATCH -e pca2.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load pcaangsd")     
cat("\n\n")

for (i in 1:26){
    cat(paste0("gzip -c /home/ktist/ph/data/new_vcf/MD7000/prune/beagle/subpop/npTB_pruned_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/prune/beagle/subpop/noTB_pruned_c",i,".BEAGLE.PL.gz \n"))
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/prune/beagle/subpop/noTB_pruned_c",i,".BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/PCAngsd/noTB_pruned_chr",i," -threads 16 \n")) 
}    

sink(NULL)


### Pairwise Fst calculations for all pairs
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)
y17<-pops[grep("17",pops)]
comb<-combn(y17, 2)
comb<-t(comb)


for (i in 1:nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    sink(paste0("Data/Slurmscripts/Fst",gsub("17",'',pop1),gsub("17",'',pop2),".sh"))
    cat("#!/bin/bash -l")
    cat("\n")
    cat(paste0("#SBATCH --job-name=Fst",i," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --ntasks=8 \n")) 
    cat(paste0("#SBATCH --nodes=4 \n"))
    cat(paste0("#SBATCH -e Fst",i,".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n\n")
    cat("module load angsd \n") 
    
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst index /home/ktist/ph/data/angsd/SFS/fromVCF/", pop1, "_maf00.saf.idx /home/ktist/ph/data/angsd/SFS/fromVCF/", pop2,"_maf00.saf.idx  " ))
    cat(paste0("-sfs /home/ktist/ph/data/angsd/SFS/fromVCF/folded_", pop1,"_",pop2,"_maf00.sfs -fstout /home/ktist/ph/data/angsd/SFS/fromVCF/2D/fst_folded_",pop1, "_",pop2,"_persite_maf00 \n"))
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS fst stats2 /home/ktist/ph/data/angsd/SFS/fromVCF/2D/fst_folded_",pop1, "_",pop2,"_persite_maf00.fst.idx -win 50000 -step 10000 > /home/ktist/ph/data/angsd/SFS/fromVCF/2D/fst_folded_",pop1, "_",pop2,"_50kWindow_maf00 \n"))
    sink(NULL)
    
}

## Create vcf file for each chromosome for PWS population:
sink(paste0("Data/Slurmscripts/subsetVCF_chrom.sh"))
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=vcfch \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
#cat(paste0("#SBATCH --nodes=4 \n"))
cat(paste0("#SBATCH -e vcfch.err  \n"))
cat(paste0("#SBATCH --time=144:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load bcftools \n") 

for (i in 1:26){
    cat(paste0("bcftools view -Oz -S /home/ktist/ph/data/new_vcf/MD7000/population/pws.txt --threads 16 /home/ktist/ph/data/new_vcf/MD7000/ph_chr",i,".vcf.gz > /home/ktist/ph/data/new_vcf/MD7000/PWS_chr",i,"_maf05.vcf.gz \n"))
               
}
sink(NULL)


#Estimate allele frequency from GL data in ANGSD

#!/bin/bash -l
sink(paste0("Data/Slurmscripts/maf01_pws.sh"))
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=maf01 \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
cat(paste0("#SBATCH --nodes=4 \n"))
cat(paste0("#SBATCH -e maf01.err  \n"))
cat(paste0("#SBATCH --time=144:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n") 

for (i in 1:length(pops)){
    cat(paste0("angsd -doGlf 2 -doMaf 1 -doMajorMinor 4 -doPost 1 -doGeno 2  -vcf-pl /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"_maf01.vcf.gz "))
    cat(paste0("-out /home/ktist/ph/data/angsd/MAF/",pops[i],"_maf01 -ref /home/jamcgirr/ph/data/c_harengus/c.harengus.fa \n"))
}
sink(NULL)

#####  Sbash sscripts to rename  files ####
files<-list.files("Data/new_vcf/PWS/", pattern = ".vcf")

sink("Data/new_vcf/PWS//renameFiles.sh")
for (i in 1: length(files)){
    nam<-files[i]
    nam<-gsub(".vcf", "",nam)
    nam2<-unlist(strsplit(nam, "_"))
    newname<-paste0(nam2[1],"_",nam2[3],"_",nam2[2],".vcf")
    cat(paste0('mv ',files[i]," ",newname, " \n"))
}
sink(NULL)

### Run PCAngsd

sink(paste0("Data/Slurmscripts/pcanagsd_pws.sh"))
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=pcanagsd_pws \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
cat(paste0("#SBATCH --nodes=4 \n"))
cat(paste0("#SBATCH -e pcanagsd_pws.err  \n"))
cat(paste0("#SBATCH --time=144:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load angsd \n") 

for (i in 1:26){
    cat(paste0("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/beagle/PWSonly_pruned_c",i,".BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/PCAngsd/PWSonly_maf05_chr",i," -threads 16 \n")) 
}
sink(NULL)


### heterozygosity stats info 
pops<-c("PWS91","PWS96","PWS07","PWS17")
windws<-read.table("Data/vcfs/windows.txt", sep="\t")
windws<-windws[windws$V1=="chr8",]
for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/het_stats8_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=het",pops[i],"\n"))
cat("#SBATCH --mem=8G
#SBATCH --ntasks=8
#SBATCH --nodes=4 
#SBATCH --cpus-per-task=12 \n")
    cat(paste0("#SBATCH --error het",pops[i],".err\n"))
cat("#SBATCH --error het.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load bcftools\n\n")
    cat(paste0("mkdir /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"/ \n"))
    for (j in 1: nrow(windws)){
        reg<-paste0(windws$V1[j],":",windws$V2[j],"-",windws$V3[j])
        cat(paste0("bcftools stats -r ", reg," -s - /home/ktist/ph/data/new_vcf/MD7000/population/PWSonly",gsub("PWS","",pops[i]), "_maf05.vcf.gz "))
        cat(paste0("| grep '^PSC' > /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"/",pops[i], "_stats_",j,"\n"))        
    }
    sink(NULL)
}
#


#add file name and concatinate all files

for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/catFiles_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat("#SBATCH --job-name=cat
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error cat.err
#SBATCH --time=24:00:00
#SBATCH -p high \n")
    cat("\n")
    cat(paste0("cd /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i],"/\n"))
    cat(paste0("for f in *; do sed -i "))
    cat('"s/$/\\t$f/" $f; done \n')
    cat(paste0("cat $(ls -t) > PWSonly_",pops[i], "_chr8_statsFile\n"))
    cat("cd ~ \n")
    sink(NULL)
}

sink(paste0("Data/Slurmscripts/indexVCF.sh"))
cat("#!/bin/bash -l")
cat("\n")
cat(paste0("#SBATCH --job-name=indexVCF \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n")) 
cat(paste0("#SBATCH --nodes=4 \n"))
cat(paste0("#SBATCH -e indexVCF.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
cat(paste0("#SBATCH --mail-type=ALL \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n\n")
cat("module load bcftools \n") 
for(i in 1:length(pops)){
    cat(paste0("bcftools index /home/ktist/ph/data/new_vcf/MD7000/population/PWSonly",gsub("PWS","",pops[i]), "_maf05.vcf.gz \n"))
    
}
sink(NULL)


### Create invariant VCF for estimating pi with Pixy 
pops<-c("PWS96","PWS07","PWS17")
for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/InvariantVCF_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=invar",pops[i],"\n"))
    cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 \n")
    cat(paste0("#SBATCH --error invar",pops[i],".err\n"))
    cat("#SBATCH --time=86:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load bcftools\n\n")
    for (j in 1:26){
        cat(paste0("bcftools mpileup -f /home/jamcgirr/ph/data/c_harengus/c.harengus.fa -b /home/ktist/ph/data/angsd/samples/",pops[i],".txt -r chr",j," | bcftools call -m -Oz -f GQ -o /home/ktist/ph/data/pixy/invariant_VCF/",pops[i], "_ch",j,".vcf.gz \n"))
    }
    sink(NULL)
}

pops<-c("TB91","TB96","TB06","TB17")
for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/InvariantVCF_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=invar",pops[i],"\n"))
    cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 \n")
    cat(paste0("#SBATCH --error invar",pops[i],".err\n"))
    cat("#SBATCH --time=86:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load bcftools\n\n")
    for (j in 1:26){
        cat(paste0("bcftools mpileup -f /home/jamcgirr/ph/data/c_harengus/c.harengus.fa -b /home/ktist/ph/data/angsd/samples/",pops[i],".txt -r chr",j," | bcftools call -m -Oz -f GQ -o /home/ktist/ph/data/pixy/invariant_VCF/",pops[i], "_ch",j,".vcf.gz \n"))
    }
    sink(NULL)
}
#
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)
sss<-c("SS96","SS06","SS17")
pops<-c(sss,pops[grep("17",pops)])
for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/InvariantVCF_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=invar",pops[i],"\n"))
    cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 \n")
    cat(paste0("#SBATCH --error invar",pops[i],".err\n"))
    cat("#SBATCH --time=86:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load bcftools\n\n")
    for (j in 1:26){
        cat(paste0("bcftools mpileup -f /home/jamcgirr/ph/data/c_harengus/c.harengus.fa -b /home/ktist/ph/data/angsd/samples/",pops[i],".txt -r chr",j," | bcftools call -m -Oz -f GQ -o /home/ktist/ph/data/pixy/invariant_VCF/",pops[i], "_ch",j,".vcf.gz \n"))
    }
    sink(NULL)
}



### Heterozygosity stats for all pops
plist<-c("group1","group2","group3")
windws<-read.table("Data/vcfs/windows.txt", sep="\t")
windws<-windws[windws$V1=="chr20",]

for (i in 1: length(plist)){
    sink(paste0("Data/Slurmscripts/het_stats20",plist[i],".sh"))
    cat("#!/bin/bash -l \n")
    cat(paste0("#SBATCH --job-name=het",plist[i]," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --ntasks=1 \n")) 
    cat(paste0("#SBATCH -e het",plist[i],".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc \n"))
    cat(paste0("#SBATCH --mail-type=ALL \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n")
    cat("module load bcftools\n\n")
    cat(paste0("mkdir /home/ktist/ph/data/new_vcf/MD7000/population/ch20", plist[i],"/ \n"))
    for (j in 1: nrow(windws)){
            reg<-paste0(windws$V1[j],":",windws$V2[j],"-",windws$V3[j])
            cat(paste0("bcftools stats -r ", reg," -s - /home/ktist/ph/data/new_vcf/MD7000/PWSonly_chr20_",plist[i], "_maf05.vcf.gz  "))
            cat(paste0("| grep '^PSC' > /home/ktist/ph/data/new_vcf/MD7000/population/ch20",plist[i],"/pws_ch20_",plist[i],"_maf05_stats_",j,"\n"))        
    }
    sink(NULL)
    }
 
#### Heterozygosity for one vcf file
### heterozygosity stats info 
windws<-read.table("Data/vcfs/windows.txt", sep="\t")
windws<-windws[windws$V1=="chr4",]

sink(paste0("Data/Slurmscripts/het_stats4_PWS.sh"))
cat("#!/bin/bash -l\n\n")
cat(paste0("#SBATCH --job-name=hetPWS \n"))
cat("#SBATCH --mem=8G
#SBATCH --ntasks=8
#SBATCH --nodes=4 
#SBATCH --cpus-per-task=12 \n")
cat(paste0("#SBATCH --error het4.err\n"))
cat("#SBATCH --error het.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load bcftools\n\n")
cat(paste0("mkdir /home/ktist/ph/data/new_vcf/MD7000/population/PWS_ch4/ \n"))
for (j in 1: nrow(windws)){
    reg<-paste0(windws$V1[j],":",windws$V2[j],"-",windws$V3[j])
    cat(paste0("bcftools stats -r ", reg," -s - /home/ktist/ph/data/new_vcf/MD7000/PWSonly_ch4_maf05.vcf.gz "))
    cat(paste0("| grep '^PSC' > /home/ktist/ph/data/new_vcf/MD7000/population/PWS_ch4/PWS_ch4_stats_",j,"\n"))        
}
sink(NULL)
#


pops<-read.csv("Data/Sample_metadata_892pops.csv")

sink(paste0("Data/Slurmscripts/NoReads.sh"))
cat("#!/bin/bash -l\n\n")
cat(paste0("#SBATCH --job-name=NoReads \n"))
cat("#SBATCH --mem=8G
#SBATCH --ntasks=8
#SBATCH --nodes=4 
#SBATCH --cpus-per-task=12 \n")
cat(paste0("#SBATCH --error NoReads.err\n"))
cat("#SBATCH --error het.err
#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load samtools\n\n")
for (i in 1: nrow(pops)){
    id<-pops$Sample[i]
    cat(paste0("samtools flagstat /home/eoziolor/phpopg/data/align/",id,".bam > /home/ktist/ph/data/",id,".stats.txt \n"))}
sink(NULL)
#


sink(paste0("Data/Slurmscripts/RunNGSadmix.sh"))
cat("#!/bin/bash -l\n\n")
cat(paste0("#SBATCH --job-name=NGSadmix \n"))
cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 ")
cat(paste0("#SBATCH --error NGSadmix.err\n"))
cat("#SBATCH --time=24:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
cat("module load angsd\n\n")


for (i in 5:10) {
    sink(paste0("Data/Slurmscripts/RunNGSadmix2.",i,".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=admix2.",i,"\n"))
    cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 \n")
    cat(paste0("#SBATCH --error admix2.",i,".err\n"))
    cat("#SBATCH --time=48:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
    cat("\n")
    cat("module load angsd\n\n")
    
    cat(paste0("NGSadmix -likes /home/ktist/ph/data/new_vcf/MD7000/beagle/PH_maf05_pruned_BEAGLE.PL.gz -K 2 -o /home/ktist/ph/data/NGSadmix/Ph_pruned_maf05_k2_run", i, "\n\n"))
    sink(NULL)
}

sink(NULL)

#Extract arntl2_a region from bam files

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pws07<-pop_info[pop_info$Population.Year=="PWS07",]

sink("Data/Slurmscripts/extract_arntl2_a_pws07.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=extract_arntl2_a \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n"))
cat(paste0("#SBATCH -e extract_arntl2_a.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load samtools")     
cat("\n\n")

for (i in 1: nrow(pws07)){
    cat(paste0('samtools view -b -h /home/eoziolor/phpopg/data/align/',pws07$Sample[i] , '.bam "chr3:3700000-3900000" >  /home/ktist/ph/data/bam/', pws07$Sample[i], '_arntl2_a_region.bam \n'))
}
sink(NULL)


# Calculate heterozygosity from sfs using angsd
pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)

sink("Data/Slurmscripts/hetero.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=hetero \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n"))
cat(paste0("#SBATCH -e hetero.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load angsd")     
cat("\n\n")

for (i in 1:length(pops)){
    cat(paste0("/home/jamcgirr/apps/angsd/misc/realSFS /home/ktist/ph/data/angsd/SFS/fromVCF/",pops[i],"_maf00.saf.idx > /home/ktist/ph/data/angsd/SFS/fromVCF/estml/",pops[i],"_maf00_est.ml \n"))
}
sink(NULL)


### heterozygosity stats info 
pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)
windws<-read.table("Data/vcfs/windows.txt", sep="\t")

for (i in 1: length(pops)){
    sink(paste0("Data/Slurmscripts/het_",pops[i],".sh"))
    cat("#!/bin/bash -l\n\n")
    cat(paste0("#SBATCH --job-name=het",pops[i],"\n"))
    cat("#!/bin/bash")
    cat("\n")
    cat(paste0("#SBATCH --job-name=het",pops[i]," \n"))
    cat(paste0("#SBATCH --mem=16G \n")) 
    cat(paste0("#SBATCH --ntasks=8 \n"))
    cat(paste0("#SBATCH -e het",pops[i],".err  \n"))
    cat(paste0("#SBATCH --time=72:00:00  \n"))
    cat(paste0("#SBATCH -p high  \n"))
    cat("\n")
    cat("module load bcftools\n\n")
    cat(paste0("mkdir /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"/ \n"))
    for (j in 1: nrow(windws)){
        reg<-paste0(windws$V1[j],":",windws$V2[j]+1,"-",windws$V3[j])
        cat(paste0("bcftools stats -r ", reg," -s - /home/ktist/ph/data/new_vcf/MD7000/population/",pops[i], "_maf00.vcf.gz "))
        cat(paste0("| grep '^PSC' > /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"/",pops[i], "_hetstats_",j,"\n"))        
    }
    sink(NULL)
}

#concatenate the files per pop
sink(paste0("Data/Slurmscripts/modify_statsFiles.sh"))
cat("#!/bin/bash -l\n\n")
cat(paste0("#SBATCH --job-name=add \n"))
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=add \n"))
cat(paste0("#SBATCH --mem=16G \n")) 
cat(paste0("#SBATCH --ntasks=8 \n"))
cat(paste0("#SBATCH -e add.err  \n"))
cat(paste0("#SBATCH --time=72:00:00  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
for (i in 1: length(pops)){
    cat(paste0("cd /home/ktist/ph/data/new_vcf/MD7000/population/", pops[i],"/ \n"))
    cat('for f in *; do sed -i "s/$/\t$f/" $f; done \n')
    cat(paste0("cat $(ls -t) > output_",pops[i],"\n"))
}
sink(NULL)




sink(paste0("Data/Slurmscripts/BeagleConvert_PCAangsd.sh"))
cat("#!/bin/bash -l\n\n")
cat(paste0("#SBATCH --job-name=BeagleConvert \n"))
cat("#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=4 \n")
cat(paste0("#SBATCH --error BeagleConvert.err\n"))
cat("#SBATCH --time=48:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")

cat("module load vcftools")     
cat("\n\n")
    
for (i in 1:26){
    cat(paste0("vcftools --gzvcf /home/ktist/ph/data/new_vcf/MD7000/population/PWS.SS.WA2017_maf05.vcf.gz  --chr chr",i))
    cat(paste0(" --out /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_c",i," --BEAGLE-PL \n"))
}

cat("\n")

#remove the head line and cat beagle files
for (i in 2:26){
    cat(paste0("sed -e '1, 1d' < /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_c",i,".BEAGLE.PL > /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_c",i,".2.BEAGLE.PL \n"))
}
cat("\n")

cat("cat /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_c1.BEAGLE.PL ") 
for (i in 2:26){
    cat(paste0("/home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_c",i,".2.BEAGLE.PL "))
}
cat(paste0(" > /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_BEAGLE.PL \n"))

cat("gzip /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_BEAGLE.PL \n\n")

cat("python /home/jamcgirr/apps/pcangsd/pcangsd.py -beagle /home/ktist/ph/data/new_vcf/MD7000/beagle/PWS.SS.WA2017_BEAGLE.PL.gz -o /home/ktist/ph/data/angsd/PCAngsd/PWS.SS.WA2017 -threads 16 \n")
sink(NULL)
    

