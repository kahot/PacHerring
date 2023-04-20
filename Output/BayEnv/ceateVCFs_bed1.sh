#!/bin/bash 

vcftools --gzvcf Data/new_vcf/PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz --bed Output/BayEnv/Run1_outlier_bed_100kbuffer.bed --out Output/BayEnv/Run1_outlier_100kbuffer --recode --keep-INFO-all 
vcftools --gzvcf Data/new_vcf/PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz --bed Output/BayEnv/Run1_outlier_bed_200kbuffer.bed --out Output/BayEnv/Run1_outlier_200kbuffer --recode --keep-INFO-all 
