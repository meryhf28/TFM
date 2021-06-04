# Trimmomatic

D=/Desktop/ArchivosRNAseq
trimmomatic=/Desktop/ArchivosRNAseq/software/Trimmomatic-0.39
list='F1M_2-1_1 F1M_2-1_2 F2M_1-1_1 F2M_1-1_2 M1M_1-1_1 M1M_1-1_2 M2M_1-1_1 M2M_1-1_2'
		
for sample in ${list}
do
	
  java -jar ${trimmomatic}/trimmomatic-0.39.jar PE -phred33 ${D}/reads/${sample}_1.fastq ${D}/reads/${sample}_2.fastq ${sample}_fp.fastq ${sample}_fu.fastq ${sample}_rp.fastq ${sample}_ru.fastq ILLUMINACLIP: ${trimmomatic}/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 
done
