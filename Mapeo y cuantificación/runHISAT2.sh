# HISAT2

D=/Desktop/ArchivosRNAseq
hisat2=/Desktop/ArchivosRNAseq/software/hisat2-2.1.0
REF=/Desktop/ArchivosRNAseq/reference/solea_v4.1_representative.fasta
list='F1M_2-1_1 F1M_2-1_2 F2M_1-1_1 F2M_1-1_2 M1M_1-1_1 M1M_1-1_2 M2M_1-1_1 M2M_1-1_2'
	

${hisat2}/hisat2-build ${REF} ${REF}
	
for sample in ${list}
do

	${hisat2}/hisat2 -p 3 --dta -x ${REF} -1 ${D}/reads/${sample}_1_fp.fastq -2 ${D}/reads/${sample}_2_fp.fastq -S ${sample}.sam
	
	samtools sort -@ 4 -o ${sample}.bam ${sample}.sam
	
done
	
rm *sam
