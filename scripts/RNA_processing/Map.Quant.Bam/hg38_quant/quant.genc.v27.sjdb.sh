#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 19:00:00
#SBATCH -J mMap.sjdb.fastq
#SBATCH -o %j.mMap.sjdb.out
#SBATCH -e %j.mMap.sjdb.err


ml GCC/5.4.0-2.26 OpenMPI/1.10.3 Subread/1.5.0-p1

files=$(ls /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_hg38_STAR_mMap.sjdb.genc.v27/*out.bam)

# Gencode exon
featureCounts -T 8 -p -B \
  --primary -s 2 -F GTF -t exon -g gene_name \
    -a /projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.exon.gtf \
    -o /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Quant/hg38.s2.multi.Gencode27.Exon.sjdb.txt $files
