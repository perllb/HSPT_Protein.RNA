#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH -J 273.m.Sjdb.STAR
#SBATCH -o 273.m.Sjdb.STAR.out
#SBATCH -e 273.m.Sjdb.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/indices/  --readFilesIn /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/273.R1.fastq.gz /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq/273.R2.fastq.gz  --readFilesCommand gunzip -c  --outFilterMismatchNoverLmax 0.03  --outFilterMultimapNmax 10  --runThreadN 8 --outSAMattributes All  --outSAMtype BAM Unsorted  --sjdbGTFfile /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/annotation/gencode/gencode.v27/gencode.v27.annotation.exon.transcript_id_hg38.liftedTo.panTro6.gtf  --outFileNamePrefix /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_pantro6_STAR_mMap.sjdb/pantro6.mMap.sjdb.genc27.lifted.273
