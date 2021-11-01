#!/bin/bash

# Request Resources
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job Name
#SBATCH -J Bowtie2_CUTRUN

#Specify Output and Input Files
#SBATCH -o outputbowtie2.out
#SBATCH -e errorbowtie2.err

#Load bowtie2
module load bowtie2
module load samtools

for read1 in *R1_001_val_1.fq.gz; do
	read2=$(echo $read1 | sed 's/R1_001_val_1.fq.gz/R2_001_val_2.fq.gz/g');
	INPUT1=$(basename ${read1} | cut -f 1 -d ".")
	INPUT2=$(basename ${read2} | cut -f 1 -d ".")
	OUTPUT=$(basename ${read1} | cut -f 1 -d "_")
	echo "input1=$INPUT1"
	echo "input2=$INPUT2"
	echo "output=$OUTPUT"
	gunzip ${INPUT1}.fq.gz
	gunzip ${INPUT2}.fq.gz
	new1=${INPUT1}.fq
	new2=${INPUT2}.fq
	echo "new=$new1"
	echo "new=$new2"
	bowtie2 -x ~/data/genomes/Mouse/GRCm39/genome -1 ${new1} -2 ${new2} -S ${OUTPUT}_unsorted.sam
	# Gzip copy of fastq
	gzip ${new1}
	gzip ${new2}
	# Convert to bam
	unsorted=${OUTPUT}_unsorted.sam
	echo "unsorted=$unsorted"
	samtools view -h -b ${unsorted} > ${OUTPUT}_unsorted.bam
	bam=${OUTPUT}_unsorted.bam
	echo "bam=$bam"
	# Sort for name
	samtools sort -n ${bam} -o ${OUTPUT}_sorted.bam
	sorted=${OUTPUT}_sorted.bam
	echo "sorted=$sorted"
	# Mark duplicates and remove unmapped:
	samtools fixmate -r -m ${sorted} ${OUTPUT}_fixmate.bam
	fixed=${OUTPUT}_fixmate.bam
	echo "fixed=$fixed"
	# Sort for position:
	samtools sort ${fixed} -o ${OUTPUT}_fixsort.bam
	fixsort=${OUTPUT}_fixsort.bam
	echo "fixsort=$fixsort"
	# Remove marked duplicates:
	samtools markdup -l 40 -r -s ${fixsort} ${OUTPUT}_dedup.bam
	dedup=${OUTPUT}_dedup.bam
	echo "dedup=$dedup"
	# Remove low quality maps
	samtools view -h -q 10 ${dedup} > ${OUTPUT}_filtered.bam
	filter=${OUTPUT}_filtered.bam
	echo "filter=$filter"
	# Sort for position again
	samtools sort ${filter} -o ${OUTPUT}_filtsort.bam
	# Indexing position sorted BAM file
	samtools index ${OUTPUT}_filtsort.bam
	# Remove sam file
	rm ${OUTPUT}_unsorted.sam
done
