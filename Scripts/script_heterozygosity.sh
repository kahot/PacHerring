#!/bin/bash

#SBATCH --job-name=heterozygosity
#SBATCH --mem=12G 
#SBATCH --ntasks=1 
#SBATCH -e heterozygosity_%A_%a.err 
#SBATCH --time=24:00:00  
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high 

module load angsd 
#realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/BC17.saf.idx >/home/ktist/ph/data/angsd/BC17_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/CA17.saf.idx >/home/ktist/ph/data/angsd/CA17_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/PWS07.saf.idx >/home/ktist/ph/data/angsd/PWS07_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/PWS17.saf.idx >/home/ktist/ph/data/angsd/PWS17_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/PWS91.saf.idx >/home/ktist/ph/data/angsd/PWS91_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/PWS96.saf.idx >/home/ktist/ph/data/angsd/PWS96_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/SS06.saf.idx >/home/ktist/ph/data/angsd/SS06_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/SS17.saf.idx >/home/ktist/ph/data/angsd/SS17_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/SS96.saf.idx >/home/ktist/ph/data/angsd/SS96_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/TB06.saf.idx >/home/ktist/ph/data/angsd/TB06_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/TB17.saf.idx >/home/ktist/ph/data/angsd/TB17_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/TB91.saf.idx >/home/ktist/ph/data/angsd/TB91_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/TB96.saf.idx >/home/ktist/ph/data/angsd/TB96_maf05_est.ml 

realSFS /home/jamcgirr/ph/data/angsd/SFS/vcf/maf05/WA17.saf.idx >/home/ktist/ph/data/angsd/WA17_maf05_est.ml 



#run: sbatch script_heterozygosity.sh