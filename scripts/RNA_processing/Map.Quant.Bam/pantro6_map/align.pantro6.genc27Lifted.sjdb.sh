#########################################
# Generate scripts for STAR alignments  #
# Map to panTro6 (Clint)                #
# unique mapping                        #
# Allow 4.5 mismatches                  #
#########################################


path='/projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq'
fastq='/projects/fs1/medpvb/backup/projects/ChimpHuman/RNAseq/fastq'

# star output dir
starout=$path/Aligned_pantro6_STAR_unique.sjdb
mkdir $starout
# script output dir
scrout=$path/scripts/pantro6_map/Align_pantro6.unique.sjdb
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
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu=11000
#SBATCH -C mem256GB
#SBATCH -J $sample.u.Sjdb.STAR
#SBATCH -o $sample.u.Sjdb.STAR.out
#SBATCH -e $sample.u.Sjdb.STAR.err

ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/indices/ \
 --readFilesIn $fastq/${sample}.R1.fastq.gz $fastq/${sample}.R2.fastq.gz \
 --readFilesCommand gunzip -c \
 --outFilterMismatchNoverLmax 0.03 \
 --outFilterMultimapNmax 1 \
 --runThreadN 8 --outSAMattributes All \
 --outSAMtype BAM Unsorted \
 --sjdbGTFfile /projects/fs1/medpvb/no_backup/genomicData/pantro/panTro6/annotation/gencode/gencode.v27/gencode.v27.annotation.exon.transcript_id_hg38.liftedTo.panTro6.gtf \
 --outFileNamePrefix $starout/pantro6.unique.sjdb.genc27.lifted.$sample" > $scrout/Map_unique_sjdb.gencv27lifted_pantro6.$sample.sh


done
