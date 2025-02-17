# chroseq-pub
- source code reference: https://github.com/Danko-Lab/proseq2.0

- install prerequisite
```
mamba env create -f proseq_env.yml
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
Copy them to the {mamba_base}/proseq/bin folder
{mamba_base} := mamba base after activating the environmen

- run
```
mamba activate proseq
bash chroseq.sh 
```
