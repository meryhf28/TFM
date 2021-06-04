# Conversión de txt a gff3
zcat solea_v4.1_rep_annot_file.txt.gz \
        | awk -F "\t" -v OFS="\t" 'NR>1 {print "##sequence-region", $1, 1, length($2)}' \
                > seqRegion.tmp
zcat solea_v4.1_rep_annot_file.txt.gz \
        | awk -F "\t" -v OFS="\t" 'NR>1 {if ($4 == "") $4 = "."; print $1, "AUTOFACT", "mRNA", 1, length($2), $4, "?", ".", "ID="$1}' > body.tmp
echo "##gff-version 3" \
        | cat - seqRegion.tmp body.tmp \
                > tmp.tmp && mv tmp.tmp solea_v4.1.gff3
rm body.tmp seqRegion.tmp
