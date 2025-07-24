#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lorincn@ccf.org
#SBATCH --job-name=task
#SBATCH --mem=50gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=bigmem

# Popcorn analysis (h2 + rg) files
cd /home/lorincn/isilon/Cheng-Noah/software/Popcorn
BASE="/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/t2d_popcorn_scores/"
RESULT="/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/t2d_popcorn_results/"
DATA="/home/lorincn/beegfs/lorincn/data/T2D/"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# EUR and others
## AFR
echo "EUR and AFR"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eur_afr_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_AFA_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eur_afr_correlation.txt

## EAS
echo "EUR and EAS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eur_eas_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_EAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eur_eas_correlation.txt

## SAS
echo "EUR and SAS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eur_sas_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_SAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eur_sas_correlation.txt

## HIS
echo "EUR and HIS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eur_his_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_HIS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eur_his_correlation.txt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# AFR and others
## EAS
echo "AFR and EAS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}afr_eas_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_AFA_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_EAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}afr_eas_correlation.txt

## SAS
echo "AFR and SAS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}afr_sas_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_AFA_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_SAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}afr_sas_correlation.txt

## HIS
echo "AFR and HIS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}afr_his_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_AFA_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_HIS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}afr_his_correlation.txt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# EAS and others
## SAS
echo "EAS and SAS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eas_sas_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_SAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eas_sas_correlation.txt

## HIS
echo "EAS and HIS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}eas_his_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_EAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_HIS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}eas_his_correlation.txt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SAS and others
## HIS
echo "SAS and HIS"
popcorn fit \
    -v 1 \
    --cfile ${BASE}sas_his_t2d_popcorn_scores.txt \
    --sfile1 ${DATA}1kgSNPs_Suzuki_SAS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    --sfile2 ${DATA}1kgSNPs_Suzuki_HIS_Metal_LDSC-CORR_Neff.v2.tsv.gz \
    ${RESULT}sas_his_correlation.txt


