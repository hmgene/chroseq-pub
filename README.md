# chroseq-pub
- source code reference: https://github.com/Danko-Lab/proseq2.0

- make a virtual env
```
mamba create -n proseq -f environment.yml
mamba activate proseq
```
- install proseq2 and dependent tools 
https://github.com/Danko-Lab/proseq2.0

The pipelines depend on several common bioinformatics tools:
```
 cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)
 fastx_trimmer (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
 seqtk (https://github.com/lh3/seqtk)
 prinseq-lite.pl (https://sourceforge.net/projects/prinseq/files/standalone/)
 bwa (https://sourceforge.net/projects/bio-bwa/files/)
 samtools version: 1.9 (http://www.htslib.org/download/)
 bedtools v2.28.0 (http://bedtools.readthedocs.org/en/latest/)
 bedops (https://bedops.readthedocs.io/en/latest/)
 bedGraphToBigWig (from the Kent source utilities http://hgdownload.cse.ucsc.edu/admin/exe/)
```
- make PATH to the above tools 
- modify hard-coded environments in chroseq.sh
```
output=/mnt/vstor/SOM_GENE_BEG33/ChRO_seq/hg38/DATA/$sample
rna_home=$output/rna
# RSEM index
x=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/hg38/ref/RSEM_hg38/GCF_000001405.39_GRCh38.p13_genomic.primary_assembly.fna
# BWA index
bwaIndex=/mnt/vstor/SOM_GENE_BEG33/mamba/db/hg38/hg38.fa

```

- run
```
bash chroseq.sh <sample> [<output_dir>]
```
