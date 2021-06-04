# Número de transcritos

D=/Desktop/ArchivosRNAseq
list='lncRNA mRNA'
			
for sample in ${list}
do
	
  grep '"1"' ${D}/comparacion/candidates_{$sample}_todos.gtf|cut -f9|cut -d';' -f1|sort|uniq -c|grep -v '1' > ${D}/comparacion/isoformas_${sample}
  awk '{print $1}' ${D}/comparacion/isoformas_${sample} > ${D}/comparacion/numero_isoformas_${sample}

done

for i in 1 2 3 4 5 6 7 8
do

  grep '$i' ${D}/comparacion/numero_isoformas_${sample} | wc -l

done

# Número de exones por transcrito

for sample in ${list}
do

  cut -f9 ${D}/comparacion/candidates_${sample}_todos.gtf|cut -d';' -f2|sort|uniq -c > ${D}/comparacion/exones_${sample}_transcritos
  awk '{print $1}' ${D}/comparacion/exones_${sample}_transcritos > ${D}/comparacion/ exones_${sample}

done

for i in 1 2 3 4 5 6 7 
do

  grep '$i' ${D}/comparacion/exones_${sample} | wc -l

done
