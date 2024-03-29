#########################################
# Generate scripts for STAR alignments  #
# Map to hg38                           #
# unique mapping                        #
# Allow 4.5 mismatches                  #
#########################################


path='/projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq'
fastq='/projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq'

# star output dir
starout=$path/Aligned_hg38_STAR_unique.sjdb.genc.v27
mkdir $starout
# script output dir
scrout=$path/scripts/mapping_2018.04.16/Align_hg38.unique.sjdb.genc.v27
mkdir $scrout

cd $fastq

files1=$(ls *R1.fastq.gz)
samples=$(echo $files1 | sed 's/.R1.fastq.gz//g')

echo $samples

for sample in $samples;

do echo "SAMPLE: $sample"

echo "#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J $sample.u.sjdb.STAR
#SBATCH -o $sample.u.sjdb.STAR.out
#SBATCH -e $sample.u.sjdb.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/common/genome/lunarc/indicies/star/human/hg38 \
 --readFilesIn $fastq/${sample}.R1.fastq.gz $fastq/${sample}.R2.fastq.gz \
 --readFilesCommand gunzip -c \
 --outFilterMismatchNoverLmax 0.03 \
 --outFilterMultimapNmax 1 \
 --runThreadN 8 --outSAMattributes All \
 --sjdbGTFfile /projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.annotation.gtf \
 --outSAMtype BAM Unsorted \
 --outFileNamePrefix $starout/hg38.unique.sjdb.genc.v27.$sample" > $scrout/Map_unique_sjdb.genc.v27_hg38.$sample.sh


done
