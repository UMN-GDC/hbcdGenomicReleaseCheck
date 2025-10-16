#!/bin/bash

# GOAL: remove remainder of controls and pilot data from the file for release
# AUTHOR: Christian Coffman
# Date Created: 25-10-09
# Last Updated: 25-10-09

# generate list of subjects to exclude
# NOTE: all deidentified IDs are purely numeric with single letter indicator on end
# so the only occurences of multiple letters in the file occur in the samples we want to exclucde


# Change: in ../HBCD.fam 
# line 414: NTC NTC  -> NTC_x NTC_x
# This was in order to make unique IDs for plink to use

#awk '/NTC/ {print $2} /IQ/ {print $2} /Pos/ {print $2} ' ../HBCD.fam > ../exclude_BR20.0.txt
echo "#FID IID PAT SEX PHENO" > ../remove_BR20.0.txt
grep -e NTC -e IQ -e Pos ../HBCD.fam >> ../remove_BR20.0.txt

# need to manually remove duplicate IDs in ../remove_BR20.0.txt

# exclude these individuals from the release data
#plink2 --bfile ../HBCD --make-bed --out HBCD_BR20.0 --exclude ../exclude_BR20.0.txt
plink2 --bfile ../HBCD --remove ../remove_BR20.0.txt --make-bed --out HBCD_BR20.0 

# checked this by searching for any occurence of two letters in a row with this pattern in vim
# /[a-zA-Z][a-zA-Z]


# generate the batch.info file
cp /projects/standard/basu_hbcd/shared/batch.info batch.info

Rscript generateBatchInfo.R
