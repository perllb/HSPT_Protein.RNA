#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 19:00:00
#SBATCH -J pantro.uMap.sjdb.fastq
#SBATCH -o %j.pantro.uMap.sjdb.out
#SBATCH -e %j.pantro.uMap.sjdb.err


ml GCC/5.4.0-2.26 OpenMPI/1.10.3 Subread/1.5.0-p1

files=$(ls /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_pantro6_STAR_unique.sjdb/*out.bam)

# Gencode exon
featureCounts -T 8 -p -B \
  --primary -s 2 -F GTF -t exon -g gene_name \
    -a /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/annotation/gencode/gencode.v27/gencode.v27.exon.gene_name_hg38.liftedTo.panTro6.gtf \
    -o /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Quant/pantro6.s2.unique.Gencode27.Exon.sjdb.txt $files
