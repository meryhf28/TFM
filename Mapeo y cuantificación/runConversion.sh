# Conversión de txt a gff3

D=/Desktop/ArchivosRNAseq
REF1=/Desktop/ArchivosRNAseq/reference/solea_v4.1_rep_annot_file.txt.gz
REF2=/Desktop/ArchivosRNAseq/reference/solea_v4.1.gtf

zcat ${REF1} \
        | awk -F "\t" -v OFS="\t" 'NR>1 {print "##sequence-region", $1, 1, length($2)}' \
                > seqRegion.tmp
zcat ${REF1}  \
        | awk -F "\t" -v OFS="\t" 'NR>1 {if ($4 == "") $4 = "."; print $1, "AUTOFACT", "mRNA", 1, length($2), $4, "?", ".", "ID="$1}' > body.tmp
echo "##gff-version 3" \
        | cat - seqRegion.tmp body.tmp \
                > tmp.tmp && mv tmp.tmp solea_v4.1.gff3
rm body.tmp seqRegion.tmp

# Conversión de gtf a gtf modificado

grep -v "^#" ${REF2} | cut -f 1-8 | sed 's/mRNA/transcript/' > to_gtf.transcript.txt
grep -v "^#" ${REF2} | cut -f 1-8 | sed 's/mRNA/exon/' > to_gtf.exon.txt

cut -f 9 ${REF2} | awk '$0 !~ "^#" {print "gene_id " $2 "; transcript_id " gensub(/solea_v4.1_unigene[[:digit:]]*/, "&.1", 1, $2) ";"}' > right_fields.txt
awk '{print $0 " exon_number " "\"1\"" ";"}' right_fields.txt > exons.txt

paste to_gtf.transcript.txt right_fields.txt > to_gtf.transcript_full.txt
paste to_gtf.exon.txt exons.txt > to_gtf.exon_full.txt
paste -d "\n" to_gtf.transcript_full.txt to_gtf.exon_full.txt > solea_v4.1.mod.gtf

