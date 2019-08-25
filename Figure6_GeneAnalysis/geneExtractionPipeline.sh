1. Align haplotigs to primary contigs with lastz. Template command below -->
lastz primaryContig.fa haplotig.fa --gfextend --hspthresh=20000 --chain --gapped --output=output.txt --format=maf+ --inner=10000 --identity=80 --strand=plus

2. Extract sequences from alignment blocks that correspond to gene coordinates
python /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/scripts/getAlignmentSequenceFromMAF_v6_updated.py /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/primaryFirstMAF/000000F_vs_000000F_001_lastzPlusStrand.txt 000000F_vs_000000F_001_lastzPlusStrand.fasta /nfs0/BB/Hendrix_Lab/Hops/version5/primary/GeneModels/augustusOutputSeqs_withMatthews/000000F.gff /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/cdsFromAlignment/indexFiles/000000F_vs_000000F_001_lastzPlusStrand.index

3. Split primary contig-specific gene sequence FASTA files into individual gene sequence FASTA files
python /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/cdsFromAlignment/macseAnalysis/scripts/splitFasta.py ../seqFileList.txt

4. Remove frameshifts with MACSE sub-program, exportAlignment
java -jar /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/cdsFromAlignment/MACSE/macse_v2.03.jar -prog exportAlignment -align /nfs0/BB/Hendrix_Lab/Hops/version5/primary/GeneAnalysis/withMatthews/singleGeneFasta/singleGenePair.fasta -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT out_NT.fasta -out_AA out_AA.fasta 

5. Further refine alignment sequences for dN/dS
python /nfs0/BB/Hendrix_Lab/Hops/GenomeAlignments/cdsFromAlignment/macseAnalysis/scripts/processExportAlignmentOutFile_v3.py ../exportAlignment/out_NT.fasta

# Proceed to dNdSPipeline.sh
