#!/bin/bash

module load minimap2
module load bedtools
module load samtools
module load vcflib
FREEBAYES="/src/freebayes-v1.3.1"
module load datamash

mkdir -p cascade/variants/frameshifts
for fasta in $(find cascade/data -name "*.fa"); do
    asm=$(basename ${fasta} | sed -e 's/.fa$//')
    minimap2 -t 32 -ax splice ${fasta} hopCDS.fasta > cascade/variants/tmp.sam
    samtools view -b -T ${fasta} cascade/variants/tmp.sam | samtools sort > cascade/variants/${asm}.CDS.bam
    
    samtools view cascade/variants/${asm}.CDS.bam | awk '$3 !~ /^[0-9][0-9][0-9][0-9][0-9][0-9]F_[0-9][0-9][0-9]$/' >cascade/variants/${asm}.tmp.CDS.sam 
    samtools view -b -T ${fasta} cascade/variants/${asm}.tmp.CDS.sam >cascade/variants/${asm}.primary.CDS.bam

    samtools index cascade/variants/${asm}.primary.CDS.bam 
    samtools faidx ${fasta}

    $FREEBAYES -m 60 -C 1 --fasta-reference ${fasta} cascade/variants/${asm}.primary.CDS.bam | 
    vcfallelicprimitives > cascade/variants/${asm}.CDS.vcf

    bedtools intersect -nonamecheck -wa -wb \ 
    -a <(bedtools bamtobed -i cascade/variants/${asm}.primary.CDS.bam | sort -k1,1 -k2,2g) \
	-b <(awk '($1 !~ /^#/) { L=length($5)-length($4); print $1 "\t" $2-1 "\t" $2-1+length($4) "\t" L; }' cascade/variants/${asm}.CDS.vcf | sort -k1,1 -k2,2g) | 
    awk '{ print $4 "\t" $10; }' | sort -k1,1 | datamash -g1 sum 2 |
    awk '{ L=($2<0 ? 0-$2 : $2); print $1 "\t" L "\t" (L%3); }' | 
    awk '($3 != 0) { print $1; }' | cut -d'-' -f1 | sort -u > cascade/${asm}.frameshift-genes.txt
done
