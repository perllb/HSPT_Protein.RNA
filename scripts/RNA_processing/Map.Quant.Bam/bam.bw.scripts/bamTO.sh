## generate scripts to convert RNAseq bam-files to bed and bw

cd /projects/fs1/medpvb/backup/projects/ChimpHuman/SeqRound1/Aligned_hg38_STAR_unique 


for bamfile in *bam

do echo $bamfile

newname=$(echo $bamfile | sed 's/Aligned.sortedByCoord.out.bam//g')

echo "#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J bamTOu_$bamfile
#SBATCH -o bamTOu_$bamfile.out
#SBATCH -e bamTOu_$bamfile.err


# Generate bed, bg and bw from sorted.bam 

module load foss/2016b
module load ucsc-tools/3.4.3
module load BEDTools/2.26.0

cd /projects/fs1/medpvb/backup/projects/ChimpHuman/SeqRound1/Aligned_hg38_STAR_unique

genomeCoverageBed -bg -split -ibam $bamfile -g /projects/fs1/common/genome/lunarc/indicies/star/human/hg38  > $newname.bed

sortBed -i $newname.bed > tmp.unq.$newname.txt
mv tmp.unq.$newname.txt $newname.bed

bedGraphToBigWig $newname.bed /projects/fs1/common/genome/lunarc/genomes/human/hg38/hg38.chrom.sizes.txt $newname.bw 
" > /projects/fs1/medpvb/backup/projects/ChimpHuman/SeqRound1/scripts/bam.bw.scripts/bamTO_bed.bw_$bamfile.sh

done

