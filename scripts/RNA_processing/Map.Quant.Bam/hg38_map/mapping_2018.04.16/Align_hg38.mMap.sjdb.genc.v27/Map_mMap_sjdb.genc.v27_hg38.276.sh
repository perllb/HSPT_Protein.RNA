#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J 276.m.sjdb.STAR
#SBATCH -o 276.m.sjdb.STAR.out
#SBATCH -e 276.m.sjdb.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/common/genome/lunarc/indicies/star/human/hg38  --readFilesIn /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/276.R1.fastq.gz /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/276.R2.fastq.gz  --readFilesCommand gunzip -c  --outFilterMismatchNoverLmax 0.03  --runThreadN 8 --outSAMattributes All  --outSAMtype BAM Unsorted  --sjdbGTFfile /projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.annotation.gtf  --outFileNamePrefix /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_hg38_STAR_mMap.sjdb.genc.v27/hg38.mMap.genc.v27.276
