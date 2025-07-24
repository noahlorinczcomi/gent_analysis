#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lorincn@ccf.org
#SBATCH --job-name=task
#SBATCH --mem=50gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=bigmem

# Popcorn 'scores' files
cd /home/lorincn/isilon/Cheng-Noah/software/Popcorn
module load python/3.8.6
REF="/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/"
BASE="/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/t2d_popcorn_scores/"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# EUR and others
## AFR
if [ ! -f "${BASE}eur_afr_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EUR \
        --bfile2 ${REF}AFR \
        --SNPs_to_store 20000 \
        ${BASE}eur_afr_t2d_popcorn_scores.txt
fi

## EAS
if [ ! -f "${BASE}eur_eas_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EUR \
        --bfile2 ${REF}EAS \
        --SNPs_to_store 20000 \
        ${BASE}eur_eas_t2d_popcorn_scores.txt
fi

## SAS
if [ ! -f "${BASE}eur_sas_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EUR \
        --bfile2 ${REF}SAS \
        --SNPs_to_store 20000 \
        ${BASE}eur_sas_t2d_popcorn_scores.txt
fi

## HIS
if [ ! -f "${BASE}eur_his_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EUR \
        --bfile2 ${REF}AMR \
        --SNPs_to_store 20000 \
        ${BASE}eur_his_t2d_popcorn_scores.txt
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# AFR and others
## EAS
if [ ! -f "${BASE}afr_eas_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}AFR \
        --bfile2 ${REF}EAS \
        --SNPs_to_store 20000 \
        ${BASE}afr_eas_t2d_popcorn_scores.txt
fi

## SAS
if [ ! -f "${BASE}afr_sas_t2d_popcorn_scores.txt" ]; then    
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}AFR \
        --bfile2 ${REF}SAS \
        --SNPs_to_store 20000 \
        ${BASE}afr_sas_t2d_popcorn_scores.txt
fi

## HIS
if [ ! -f "${BASE}afr_his_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}AFR \
        --bfile2 ${REF}AMR \
        --SNPs_to_store 20000 \
        ${BASE}afr_his_t2d_popcorn_scores.txt
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# EAS and others
## SAS
if [ ! -f "${BASE}eas_sas_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EAS \
        --bfile2 ${REF}SAS \
        --SNPs_to_store 20000 \
        ${BASE}eas_sas_t2d_popcorn_scores.txt
fi

## HIS
if [ ! -f "${BASE}eas_his_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}EAS \
        --bfile2 ${REF}AMR \
        --SNPs_to_store 20000 \
        ${BASE}eas_his_t2d_popcorn_scores.txt
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SAS and others
## HIS
if [ ! -f "${BASE}sas_his_t2d_popcorn_scores.txt" ]; then
    popcorn compute \
        -v 1 \
        --bfile1 ${REF}SAS \
        --bfile2 ${REF}AMR \
        --SNPs_to_store 20000 \
        ${BASE}sas_his_t2d_popcorn_scores.txt
fi

