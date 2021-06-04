# HTSeq-count
		
D=/Desktop/ArchivosRNAseq
REFgff=/Desktop/ArchivosRNAseq/annotation/solea_v4.1.mod.gtf
HTSeq=/Desktop/ArchivosRNAseq/software/htseq-count-0.13.5
list='F1M_2-1_1 F1M_2-1_2 F2M_1-1_1 F2M_1-1_2 M1M_1-1_1 M1M_1-1_2 M2M_1-1_1 M2M_1-1_2'
		
for sample in ${list}
do
	
  ${HTSeq}/htseq-count -f bam ${D}/mapping/${sample}.bam --stranded=no ${REFgff} > ${sample}.gtf 
	  
done
