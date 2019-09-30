#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 19:00:00
#SBATCH -J uMap.NOsjdb.fastq
#SBATCH -o %j.uMap.NOsjdb.out
#SBATCH -e %j.uMap.NOsjdb.err

ml GCC/5.4.0-2.26 OpenMPI/1.10.3 Subread/1.5.0-p1

files=$(ls /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Aligned_hg38_STAR_unique.NOsjdb/*out.bam)

# Gencode exon
featureCounts -T 8 -p -B \
  --primary -s 2 -F GTF -t Retroelement -g RepName \
    -a /projects/fs1/medpvb/no_backup/genomicData/hg38/RepeatMasker/not_in_exon.gencode/Retro.hg38.uniqID.NOT.IN.EXON.gencode.gtf \
    -o /projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/Quant/hg38.s2.unique.RepeatMasker.NOsjdb.txt $files
