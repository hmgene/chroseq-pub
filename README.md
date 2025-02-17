# chroseq-pub
- source code reference: https://github.com/Danko-Lab/proseq2.0

- install prerequisite

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

- run
```
bash chroseq.sh 
```
