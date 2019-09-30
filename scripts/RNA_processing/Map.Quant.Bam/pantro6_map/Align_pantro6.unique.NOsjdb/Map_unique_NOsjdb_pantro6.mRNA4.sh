#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J mRNA4.u.noS.STAR
#SBATCH -o mRNA4.u.noS.STAR.out
#SBATCH -e mRNA4.u.noS.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/  --readFilesIn /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/mRNA4.R1.fastq.gz /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/mRNA4.R2.fastq.gz  --readFilesCommand gunzip -c  --outFilterMismatchNoverLmax 0.03  --outFilterMultimapNmax 1  --runThreadN 8 --outSAMattributes All  --outSAMtype BAM Unsorted  --outFileNamePrefix /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_pantro6_STAR_unique.NOsjdb/pantro6.unique.NOsjdb.mRNA4
