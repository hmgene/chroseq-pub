#!/bin/bash 
# 02/21/2024 Modified by Hyunmin Kim (hxk728@case.edu)
BASEDIR=$( cd "$(dirname "$0")" ; pwd -P )
export PATH=$PATH:$BASEDIR/bin


############################################
# Declaring variables
HPC='/mnt/vstor/SOM_GENE_BEG33'
ENVIRONMENT=${HPC}/RNA_seq/hg38
SLURM_CPUS_PER_TASK=10

MYHOME=${ENVIRONMENT}/DATA
MYHOME=`realpath ${2:-$MYHOME}`; 
MANAGE=${ENVIRONMENT}/manage_samples
SCRIPT=${ENVIRONMENT}/scripts
IGVTOOLS=${HPC}/software/IGVTools_2.4.19_rs38

REF=${ENVIRONMENT}/ref/RSEM_hg38/star_hg38
RSEM=${ENVIRONMENT}/ref/RSEM_hg38/refseq_hg38/refseq_hg38
CHR=${ENVIRONMENT}/ref/RSEM_hg38/transcripts_table.txt
CHRGENELEVEL=${ENVIRONMENT}/ref/RSEM_hg38/genes_table.txt
GENOME=${ENVIRONMENT}/ref/RSEM_hg38/GCF_000001405.39_GRCh38.p13_genomic.primary_assembly.fna

SAMPLE=$MYHOME/$1
BAM="$1.bam"
TBAM="$1.rs.bam"

export PATH="/mnt/vstor/SOM_GENE_BEG33/mamba/miniforge3/bin:$PATH"
source activate proseq

#############################################
echo "Sample: $SAMPLE/$1_R1.fastq.gz"
cd $SAMPLE
#########################

# map BAMS if not done previously
if [ -s "$SAMPLE/$1.bam" -a -s "$SAMPLE/$1.rs.bam" ]
then 
	echo "mapping completed previously"
else
	if [ -f "$SAMPLE/${1}_R1.fastq.gz" -a -f "$SAMPLE/${1}_R2.fastq.gz" ]
	then
		echo "running STAR in paired-end mode"
		echo "with $SAMPLE/$1_R1.fastq.gz and $SAMPLE/$1_R2.fastq.gz"
		STAR --genomeDir $REF --readFilesIn $SAMPLE/$1_R1.fastq.gz $SAMPLE/$1_R2.fastq.gz --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --chimSegmentMin 12  --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10  --alignMatesGapMax 100000  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2 --outSAMunmapped Within --quantMode TranscriptomeSAM --limitSjdbInsertNsj 10000000 --outBAMsortingThreadN 8  --limitBAMsortRAM 8000000000 --outBAMsortingBinsN 100

	else
		echo "running STAR in single-end mode"
		STAR --genomeDir $REF --readFilesIn  $MYHOME/$1/$1_R1.fastq.gz  --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimSegmentMin 12  --chimJunctionOverhangMin 12  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2  --outSAMunmapped Within --quantMode TranscriptomeSAM  --limitSjdbInsertNsj 10000000 
	fi

	if [ -f "$1Log.final.out" ]
	then
		echo "mapping completed"
	else
		echo "mapping incomplete, please make sure the input files are in the $DATA directory "
		exit;
	fi
	samtools sort $1Aligned.out.bam > $1Aligned.sortedByCoord.out.bam

	mv $1Aligned.toTranscriptome.out.bam $1.rs.bam
	mv $1Aligned.sortedByCoord.out.bam $1.bam
fi

# count transcripts and genes if not done previously
if [ -f "$SAMPLE/$1.genes.results" ]
then 
	echo "RSEM counting completed previously"
else
	conda deactivate 
	if [ -f "$SAMPLE/$1_R1.fastq.gz" -a -f "$SAMPLE/$1_R2.fastq.gz" ]

	then
		rsem-calculate-expression --no-bam-output --paired-end -p $SLURM_CPUS_PER_TASK --estimate-rspd --bam $TBAM $RSEM $1
	else 
		rsem-calculate-expression --no-bam-output -p $SLURM_CPUS_PER_TASK --estimate-rspd --bam $TBAM $RSEM $1
	fi
 
#	check if complete
	if [ -f "$SAMPLE/$1.genes.results" ]
	then
		echo "Counts generated"
	else
		echo "Please check the input bam file"
	fi
fi

# gene level to EDEN format
join  -t $'\t' <(sort -k1,1  $CHRGENELEVEL) <(sort -k1,1 $1.genes.results) | tac > $1_geneswithcoordinates.txt
	echo "sort and joined chromosomal locations to gene level results"
awk 'BEGIN {FS="\t"; OFS="\t"} {print $2,$3,$4,$1,$9}' $1_geneswithcoordinates.txt  > $1_genetemp.txt
	echo "rearranged column order to temp gene file"
{ printf 'Chr\tStart\tStop\tGeneID\tTPM\n';cat $1_genetemp.txt ; } > $1.gene.TPM.txt

rm $1_genetemp.txt $1_geneswithcoordinates.txt

# isoforms to COLTRON format
awk 'NR == 1; NR > 1 {print $0 |"sort -k1"}' $1.isoforms.results > $1sorted.isoform.results
	echo "sorted isoform level results"
join  -t $'\t' <(sort -k1,1 $CHR) <(sort -k1,1 $1sorted.isoform.results) | tac > $1_isoformswithcoordinates.txt
	echo "joined chromosomal locations to isoform level results"
awk 'BEGIN {FS="\t"; OFS="\t"} {print $3,$4,$5,$2,$1,$10}' $1_isoformswithcoordinates.txt  > $1_temp.txt 
	echo "rearranged column order to temp isoform file"
{ printf 'Chr\tStart\tStop\tGeneID\tTranscriptID\tTPM\n';cat $1_temp.txt ; } > $1.transcript.TPM.txt

rm $1_temp.txt $1_isoformswithcoordinates.txt $1sorted.isoform.results

# index BAM and make TDF if not done previously
if [ -f "$MYHOME/$1/$1.tdf" ]
then 
	echo "TDF made previously"
else
	samtools index $BAM
	java -jar $IGVTOOLS/lib/igvtools.jar count $BAM  $1.tdf $GENOME
	echo "tdf file generated"
fi

# make RPM tdf
if [ -f "$MYHOME/$1/$1.RPM.tdf" ]
then 
	echo "RPM scaled TDF made previously"
else

	total_mapped_reads=`samtools view -c -F 260 $BAM`
	echo "total mapped reads: $total_mapped_reads"
	scale_factor=`echo "scale=5; 1000000/$total_mapped_reads" | bc `
	echo "scale TDF factor: $scale_factor"
  	
	perl $BASEDIR/scaleTDF.pl -i $1.tdf -o $1.RPM.tdf -f $scale_factor  
	echo "RPM scaled TDF complete"
fi

###################################	
#mv $1.gene.TPM.txt $SAMPLE
#mv $1.transcript.TPM.txt $SAMPLE
#mv $1.genes.results $SAMPLE
#mv $1.isoforms.results $SAMPLE
#mv $1.tdf $SAMPLE
#mv $1.RPM.tdf $SAMPLE
#mv $1.bam.bai $SAMPLE
#mv $1.bam $SAMPLE
###################################

chmod -R 775 $SAMPLE
chgrp beg33 -R $SAMPLE
echo "permissions changed"
#mv $MANAGE/slurm-$SLURM_JOBID.out $MANAGE/slurm.out/
echo "pipeline completed!"
