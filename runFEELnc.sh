#FEELnc
			
D=/Desktop/ArchivosRNAseq
REFgff=/Desktop/ArchivosRNAseq/annotation/solea_v4.1.mod.gtf
REF=/Desktop/ArchivosRNAseq/reference/solea_v4.1_representative.fasta
FEELnc=/Desktop/ArchivosRNAseq/software/feelnc_install_dir
	
# FEELnc_filter.pl		

${FEELnc}/FEELnc_filter.pl -i ${D}/assembly/all.gtf -a ${REFgff} > all_candidates.gtf 
		  
# Modificación del transcriptoma de referencia: Se cambias las ‘x’ como nucleótido por ‘n’.
awk '$0 ~ "^>" {print $0}; $0 !~ "^>" {print gensub(/x/,"N","G")}' ${REF} > ${D}/reference/res.fasta

# FEELnc_codpot.pl
${FEElnc}/FEElnc_codpot.pl -i all_candidates.gtf -a ${REFgff} -g ${D}/reference/res.fasta –outdir=".feelnc_codpot_out" --mode=shuffle

# FEELnc_classifier.pl
${FEElnc}/FEElnc_classifier.pl -i ./feelnc_codpot_out/all_candidates.gtf.lncRNA.gtf -a ${REFgff} > all_candidates.gtf.lncRNA_classes.txt
		

