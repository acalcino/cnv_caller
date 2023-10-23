# cnv_caller

cnv_caller is a bunch of wrappers for Numbat (https://github.com/kharchenkolab/numbat) that simplifies some of the options. It also includes a bit of pre-processing and does a few calculations on the output. The idea is to use the python scripts in this order:

1. pileup_and_phaser.py
2. cnv_caller.py
3. bam_subsetter.py
4. finalser.py

The 'resources' folder contains the pileup_and_phase.R script from Numbat as well as the hg38 genetic map from Eagle (https://alkesgroup.broadinstitute.org/Eagle/). If you don't already have the 1000 Genomes reference panel located in the 'resources' folder, this will be downloaded when you run pileup_and_phaser.py for the first time. Same thing for the genome1k snp vcf file (genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf). 
