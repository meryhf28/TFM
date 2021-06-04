# Stringtie
	
D=/Desktop/ArchivosRNAseq
Stringtie=/Desktop/ArchivosRNAseq/software/stringtie-2.1.1
list='F1M_2-1_1 F1M_2-1_2 F2M_1-1_1 F2M_1-1_2 M1M_1-1_1 M1M_1-1_2 M2M_1-1_1 M2M_1-1_2'
	
# Obtener un gtf con las posiciones de los transcritos para cada muestra

for sample in ${list}
do

	${Stringtie}/stringtie ${D}/mapping/${sample}.bam -p 3 -o ${sample}.gtf 
  
done

# Combinar los archivos gtf individuales en un único gtf que contenga las transcripciones de todas las muestras

${Stringtie}/stringtie --merge ${sample}.gtf -p 3 -o ./all.gtf 
