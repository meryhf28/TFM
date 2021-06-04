# FASTA
			
D=/Desktop/ArchivosRNAseq
REF=/Desktop/ArchivosRNAseq/reference/solea_v4.1_unigenes.fasta
bedtools=/Desktop/ArchivosRNAseq/software/bedtools-2.27.1
FEELnc=/Desktop/ArchivosRNAseq/feelnc_codpot_out
	
# Se crea el fichero intervalos.bed para extraer la información de las coordenadas del fichero all.candidates.gtf.lncRNA.gtf:
awk -F "\t" -v OFS="\t" '{print $1, $4 - 1, $5}' ${FEELnc}/all_candidates.gtf.lncRNA.gtf > ${D}/fasta/intervals.bed 

# Se extraen las secuencias:
${bedtools}/bedtools getfasta -fi ${REF} -bed {D}/fasta/intervals.bed -fo {D}/fasta/tmp.fasta 

# Se crea el fichero unigene_pos_transcriptId.txt que combina los identificadores de secuencias "unigene" y gen de lncRNA:
awk -F "\t" -v OFS="\t" '{printf(">%s:%d-%d\n", $1, $4 - 1, $5)}' ${FEELnc}/all_candidates.gtf.lncRNA.gtf > {D}/fasta/unigene_pos.tmp
awk -F "\t" -v OFS="\t" '{print $NF}' ${FEELnc}/all_candidates.gtf.lncRNA.gtf | awk '{print $4}' | tr -d '";' > {D}/fasta/transcriptId.tmp 
paste {D}/fasta/unigene_pos.tmp {D}/fasta/transcriptId.tmp > {D}/fasta/unigene_pos_transcriptId.txt
rm {D}/fasta/unigene_pos.tmp {D}/fasta/transcriptId.tmp

# Se obtiene el fichero description.tmp con el transcript.id:
cut -f 1 {D}/fasta/unigene_pos_transcriptId.txt | cut -d ":" -f 1 > {D}/fasta/left.tmp
cut -f 2 {D}/fasta/unigene_pos_transcriptId.txt | cut -d "." -f 2 > {D}/fasta/rigth.tmp
paste -d "_" {D}/fasta/left.tmp {D}/fasta/rigth.tmp > {D}/fasta/unigene_names.txt
cut -f 2 {D}/fasta/unigene_pos_transcriptId.txt | cut -d "." -f 3 > {D}/fasta/mstrg_id.txt
paste {D}/fasta/unigene_names.txt {D}/fasta/mstrg_id.txt | sed 's/\t/./' > {D}/fasta/description.tmp 

# Se obtiene el fichero sequence.tmp, se pegan la descripción y la secuencia y se obtiene el fichero lncRNA.fasta:
sed -n '2~2p' {D}/fasta/tmp.fasta > {D}/fasta/sequence.tmp
paste {D}/fasta/description.tmp {D}/fasta/sequence.tmp > {D}/fasta/desc_seq.tmp
awk -F "\t" -v OFS="\t" '{if (array[$1]) array[$1] = array[$1]$2; else array[$1] = $2}; END {for (i in array) printf("%s\n%s\n", i, array[i])}' {D}/fasta/desc_seq.tmp > {D}/fasta/tmp.renamed.fasta
mv {D}/fasta/tmp.renamed.fasta {D}/fasta/lncRNA.fasta
rm {D}/fasta/desc_seq.tmp {D}/fasta/sequence.tmp {D}/fasta/description.tmp {D}/fasta/tmp.fasta
rm {D}/fasta/mstrg_id.txt {D}/fasta/unigene_names.txt 
rm {D}/fasta/rigth.tmp {D}/fasta/left.tmp

# Se obtiene el fichero lncRNA_longest.fasta:
bioawk -c fastx '{print ">" $name "\t" $seq "\t" length($seq)}' {D}/fasta/lncRNA.fasta | sort -k1b,1 > {D}/fasta/seqs.sorted 
cut -f 1 -d " " {D}/fasta/seqs.sorted | cut -d "." -f 1,2 > {D}/fasta/gene_id.txt 
paste {D}/fasta/seqs.sorted {D}/fasta/gene_id.txt | sort -k4b,4 -k3nr,3 | uniq -f 3 | sort -k4V,4 | awk '{print $4 "\n" $2}' | cut -d '_' -f 1,2,3 > {D}/fasta/lncRNA_longest.fasta 
rm {D}/fasta/seqs.sorted {D}/fasta/gene_id.txt
		
