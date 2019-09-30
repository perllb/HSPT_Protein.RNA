#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J 260.u.noS.STAR
#SBATCH -o 260.u.noS.STAR.out
#SBATCH -e 260.u.noS.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/common/genome/lunarc/indicies/star/human/hg38  --readFilesIn /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/260.R1.fastq.gz /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/260.R2.fastq.gz  --readFilesCommand gunzip -c  --outFilterMismatchNoverLmax 0.03  --outFilterMultimapNmax 1  --runThreadN 8 --outSAMattributes All  --outSAMtype BAM Unsorted  --outFileNamePrefix /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_hg38_STAR_unique.NOsjdb/hg38.unique.NOsjdb.260
