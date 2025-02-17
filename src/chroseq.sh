#!/bin/bash
usage="$BASH_SOURCE <sample_name> [output]"
if [ $# -lt 1 ];then echo "$usage";exit ;fi

BASEDIR=$( cd "$(dirname "$0")" ; pwd -P )
export PATH=$PATH:$BASEDIR


#get_ChROseq_data.sh
#Creates folder and links fastq
#run_ChROseq.sh
#1. Set up names and variables
#2. Load modules for dependencies
#3. Checks for BAM.  If it doesn’t exist, then check for 
#fastqs R1 and R2. If both exist, runs proseq2.1.bsh in 
#Paired End mode. Otherwise, run proseq2.0.bsh in 
#Single End mode.
#4. Runs RNAseq pipeline to make a new bam file and 
#also to quantify transcripts per million
#1. Makes stranded BAMs
#2. Makes stranded TDF files
#!/bin/bash



#rna_home=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/hg38/DATA
chromInfo=`realpath $BASEDIR/chromInfo.txt`


bwaIndex=/mnt/vstor/SOM_GENE_BEG33/mamba/db/hg38/hg38.fa
input=$1;
if [ ! -f $input ];then 
	input=/mnt/vstor/SOM_GENE_BEG33/data/*/${1%_R1.fastq.gz}_R1.fastq.gz; ## find in the structure 
	if [ ! -f $input ];then echo "$input does not exists";exit;fi
fi
input=`realpath $input`
sample=${1##*/};sample=${sample%_R1.*}
output=${2:-/mnt/vstor/SOM_GENE_BEG33/ChRO_seq/hg38/DATA/$sample}; mkdir -p $output; output=`realpath $output`;
rna_home=$output/rna 
source activate proseq
o=$output/temp/tmp/passQC/${sample}_dedup_QC_end_1.fastq.gz
if [ ! -f $o -a  -f $input -a -f "${input/_R1/_R2}" ];then
	#bash $BASEDIR/proseq2.0.bsh \
	echo "paired"
	proseq2.0 \
	--thread=20 -i $bwaIndex -c $chromInfo -PE -T $output/temp -O $output --UMI1=6 -I ${input%_R1.*} \
	--RNA3=R1_5prime --UMI1=6 --ADAPT2=TGGAATTCTCGGGTGCCAAGG --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC -3 
elif [ ! -f $o -a -f $input -a ! -f "${input/_R1/_R2}" ];then
	#bash $BASEDIR/proseq2.0.bsh \
	echo "single"
	proseq2.0 \
	--thread=20 -i $bwaIndex -c $chromInfo -SE -G -T $output/temp -O $output --UMI1=6 -I ${input%\.f*}
else
	echo "skip => $o"
fi

i=$output/temp/tmp/passQC/${sample}_dedup_QC_end_1.fastq.gz 
o=$rna_home/${sample}_trim/${sample}_trim_R1.fastq.gz
i2=${i/1.fastq/2.fastq}; o2=${o/1.fastq/2.fastq}
if [ ! -h $o ];then
	mkdir -p ${o%/*}; ln -s $i $o; 
	if [ -f $i2 ];then ln -s $i2 $o2; fi
fi

#bash $BASEDIR/bin/rnaseq_hg38.sh ${sample}_trim $rna_home 
o=$rna_home/${sample}_trim/${sample}_trim;
if [ ! -f $o.bam ];then
	$BASEDIR/rnaseq_hg38 ${sample}_trim $rna_home 
fi
if [ -f "$o.tdf" ]
then	
	echo "RNA-seq pipeline completed for ChROseq $o, will split by strand."

		# Set variables
		BAM=$o.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.bam"
		BAMF1=$o.F1.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.F1.bam"
		BAMF2=$o.F2.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.F2.bam"
		BAMF=$o.F.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.F.bam"
		BAMR1=$o.R1.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.R1.bam"
		BAMR2=$o.R2.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.R2.bam"
		BAMR=$o.R.bam #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.R.bam"
		TDFF=$o.F.tdf #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.F.tdf"
		TDFR=$o.R.tdf #"/data/khanlab/projects/ChIP_seq/RNA_DATA/$TRIMSAMPLE/$TRIMSAMPLE.R.tdf"
		
		if [ -f "$BAMR" ]
		then	
			echo "$BAM split by strand previously."
		else
			# Forward strand. 0x10 - SEQ being reverse complemented
				# 1. alignments of the second in pair if they map to the forward strand
				# 2. alignments of the first in pair if they map to the reverse  strand
				samtools view -b -f 128 -F 16 $BAM > $BAMF1
				samtools index $BAMF1

				samtools view -b -f 80 $BAM > $BAMF2
				samtools index $BAMF2

				# Combine alignments that originate on the forward strand.
				samtools merge -f $BAMF $BAMF1 $BAMF2
				samtools index $BAMF
				
			# Reverse strand
				# 1. alignments of the second in pair if they map to the reverse strand
				# 2. alignments of the first in pair if they map to the forward strand
				samtools view -b -f 144 $BAM > $BAMR1
				samtools index $BAMR1

				samtools view -b -f 64 -F 16 $BAM > $BAMR2
				samtools index $BAMR2

				# Combine alignments that originate on the reverse strand.
				samtools merge -f $BAMR $BAMR1 $BAMR2
				samtools index $BAMR
			
			#clean up
			echo "cleaning up: rm $BAMR1 $BAMR2 $BAMF1 $BAMF2 ${BAMR1}.bai ${BAMR2}.bai ${BAMF1}.bai ${BAMF2}.bai"
			rm $BAMR1 $BAMR2 $BAMF1 $BAMF2 ${BAMR1}.bai ${BAMR2}.bai ${BAMF1}.bai ${BAMF2}.bai
			echo "F and R BAMs completed and indexed."
		fi
		
		# making IGV viewable TDFs
		echo "Making F and R TDFs now."
		#x=`realpath $BASEDIR/ncbi_chrom.txt`
		x=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/hg38/ref/RSEM_hg38/GCF_000001405.39_GRCh38.p13_genomic.primary_assembly.fna

		igvtools count $BAMF $TDFF $x		
		igvtools count $BAMR $TDFR $x		
		#java -jar /usr/local/apps/igvtools/2.3.98/igvtools.jar count $BAMF  $TDFF $GENOME
		#java -jar /usr/local/apps/igvtools/2.3.98/igvtools.jar count $BAMR  $TDFR $GENOME
		
		chgrp -R beg33 $o* 
		chmod -R g+w $o*
		echo "F and R permissions opened"

else
	echo "RNA-seq pipeline not completed yet for $o.tdf" 
	echo "Rerun ChRO-seq pipeline later to split trim BAMS by strand" 
fi
		for f in $o*.{tdf,txt};do
			n2=${f##*/};
			echo " cp $f $output/${n2/_trim/} "
			cp $f $output/${n2/_trim/}
		done


#pipeline end, open permissions!
mamba deactivate
chmod -R 775 ${input%/*}
echo "permissions changed"
echo "ChROseq pipeline completed!"


