#!/bin/bash

/home/hui/software/bwa-0.7.10/bwa mem -T 20 -t 16 sbgenome/Sbicolor_255_v2.0 sc03_sc.mlc.cns.sgl.utg.scga7.importdb.fa > sc03_sc.mlc.cns.sgl.utg.scga7.importdb.sam 
/home/hui/software/samtools-1.1/samtools view -bhS sc03_sc.mlc.cns.sgl.utg.scga7.importdb.sam > sc03_sc.mlc.cns.sgl.utg.scga7.importdb.bam 
/home/hui/software/samtools-1.1/samtools sort -m 12000000000 sc03_sc.mlc.cns.sgl.utg.scga7.importdb.bam sc03_sc.mlc.cns.sgl.utg.scga7.importdb.sorted
/home/hui/software/samtools-1.1/samtools mpileup -g -C 50 -f sbgenome/Sbicolor_255_v2.0 sc03_sc.mlc.cns.sgl.utg.scga7.importdb.sorted.bam -o sc03_sc.mlc.cns.sgl.utg.scga7.importdb.bcf 
/home/hui/software/bcftools-1.1/bcftools call -m -v -O v -o sc03_sc.mlc.cns.sgl.utg.scga7.importdb.vcf sc03_sc.mlc.cns.sgl.utg.scga7.importdb.bcf
/home/hui/software/bcftools-1.1/bcftools call -m -O v -o sc03_sc.mlc.cns.sgl.utg.scga7.importdb.allsites.vcf sc03_sc.mlc.cns.sgl.utg.scga7.importdb.bcf 
/home/hui/software/samtools-1.1/samtools mpileup -C 50 -f sbgenome/Sbicolor_255_v2.0 sc03_sc.mlc.cns.sgl.utg.scga7.importdb.sorted.bam > sc03_sc.mlc.cns.sgl.utg.scga7.importdb.pileup
