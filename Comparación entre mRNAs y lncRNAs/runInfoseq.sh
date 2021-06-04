# Infoseq

D=/Desktop/ArchivosRNAseq
emboss=/Desktop/ArchivosRNAseq/software/emboss-6.6.0.0
list='lncRNA_longest mRNA_longest'
		
for sample in ${list}
do

  ${emboss}/infoseq ${D}/fasta/${sample}.fasta -auto -only -name -length -pgc > ${sample}.characteristics.txt
  
done
