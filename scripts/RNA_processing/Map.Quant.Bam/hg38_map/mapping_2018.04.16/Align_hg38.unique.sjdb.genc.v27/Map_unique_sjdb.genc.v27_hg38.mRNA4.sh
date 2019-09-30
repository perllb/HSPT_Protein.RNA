#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J mRNA4.u.sjdb.STAR
#SBATCH -o mRNA4.u.sjdb.STAR.out
#SBATCH -e mRNA4.u.sjdb.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/common/genome/lunarc/indicies/star/human/hg38  --readFilesIn /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/mRNA4.R1.fastq.gz /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/mRNA4.R2.fastq.gz  --readFilesCommand gunzip -c  --outFilterMismatchNoverLmax 0.03  --outFilterMultimapNmax 1  --runThreadN 8 --outSAMattributes All  --sjdbGTFfile /projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.annotation.gtf  --outSAMtype BAM Unsorted  --outFileNamePrefix /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_hg38_STAR_unique.sjdb.genc.v27/hg38.unique.sjdb.genc.v27.mRNA4
