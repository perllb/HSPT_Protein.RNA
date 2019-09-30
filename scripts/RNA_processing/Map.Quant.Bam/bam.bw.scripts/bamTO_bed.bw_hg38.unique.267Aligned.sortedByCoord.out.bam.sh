#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J bamTOu_hg38.unique.267Aligned.sortedByCoord.out.bam
#SBATCH -o bamTOu_hg38.unique.267Aligned.sortedByCoord.out.bam.out
#SBATCH -e bamTOu_hg38.unique.267Aligned.sortedByCoord.out.bam.err


# Generate bed, bg and bw from sorted.bam 

module load foss/2016b
module load ucsc-tools/3.4.3
module load BEDTools/2.26.0

cd /projects/fs1/medpvb/backup/projects/ChimpHuman/SeqRound1/Aligned_hg38_STAR_unique

genomeCoverageBed -bg -split -ibam hg38.unique.267Aligned.sortedByCoord.out.bam -g /projects/fs1/common/genome/lunarc/indicies/star/human/hg38  > hg38.unique.267.bed

sortBed -i hg38.unique.267.bed > tmp.unq.hg38.unique.267.txt
mv tmp.unq.hg38.unique.267.txt hg38.unique.267.bed

bedGraphToBigWig hg38.unique.267.bed /projects/fs1/common/genome/lunarc/genomes/human/hg38/hg38.chrom.sizes.txt hg38.unique.267.bw 

