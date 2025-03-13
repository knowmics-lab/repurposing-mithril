
#EX_DRUG_FILE='zuclopenthixol_24h.perturbations.txt'
#f3=gene id f4=gene_name
GENE_NAMES=$(more  'zuclopenthixol_24h.perturbations.txt'  | cut -f4  | sort -u |  head -n10)

#sort sorts, -u only keeps uniques. Can't do the same for gene_IDS because the order would change by sorting them
mkdir ../bash_output



# Initialize gene-wise metanalysis files with header:
for GENE_NAME in $GENE_NAMES; do
	GENE_ID=$(more  'zuclopenthixol_24h.perturbations.txt'  | grep $GENE_NAME | cut -f3  | head -n1)
	METANALYSIS_FILENAME=$GENE_NAME'_'$GENE_ID'_metanalysis.txt'
	# Create header for gene wise metanalysis file with mith3 data:
	echo 'gene'$'\t''drug'$'\t''DE_log2_FC_6h'$'\t''std.error_6h' $'\t''t.value_6h'$'\t'\
	         'p.value_6h'$'\t' 'adj.p.value_6h' $'\t' 'DE_log2_FC_24h'$'\t' 'std.error_24h'$'\t' 't.value_24h'$'\t'\
	         'p.value_24h'$'\t' 'adj.p.value_24h'$'\t''DE_log2_FC_6h_24h'$'\t''std.error_6h_24h'$'\t''t.value_6h_24h'$'\t''p.value_6h_24h'\
	          > ../bash_output/$METANALYSIS_FILENAME

done


for drug_file in *6h_24h.perturbations.txt;do

  # extract drug string from drug file name:
  arrIN=(${drug_file//_/ });
  drug=$(echo ${arrIN[0]}); 
  echo $drug;

  for GENE_NAME in $GENE_NAMES; do
	  # Take relevant values from line with corresponding gene. Every gene data is repeated for multiple pathways, so only take the first appearance
	  pert_24h=$(grep $GENE_NAME $drug'_24h.perturbations.txt' | cut -f5 | head -n1);
		pvalue_24h=$(grep $GENE_NAME $drug'_24h.perturbations.txt'| cut -f7 | head -n1);
		adj_pvalue_24h=$(grep $GENE_NAME $drug'_24h.perturbations.txt' | cut -f8 | head -n1);

		pert_6h=$(grep $GENE_NAME $drug'_6h.perturbations.txt'  | cut -f5 | head -n1);
		pvalue_6h=$(grep $GENE_NAME $drug'_6h.perturbations.txt' | cut -f7 | head -n1);
		adj_pvalue_6h=$(grep $GENE_NAME $drug'_6h.perturbations.txt' | cut -f8 | head -n1);

		pert_6h_24h=$(grep $GENE_NAME $drug'_6h_24h.perturbations.txt'| cut -f5 | head -n1);
		pvalue_6h_24h=$(grep $GENE_NAME $drug'_6h_24h.perturbations.txt' | cut -f7 | head -n1);

		# Append line to gene-wise file

		#probably not very efficient to create these variables again
		GENE_ID=$(more  'zuclopenthixol_24h.perturbations.txt'  | grep $GENE_NAME | cut -f3  | head -n1)
		METANALYSIS_FILENAME=$GENE_NAME'_'$GENE_ID'_metanalysis.txt'


		echo $GENE_ID$'\t'$drug$'\t'$pert_6h$'\t''NA'$'\t''NA'$'\t'$pvalue_6h$'\t'$adj_pvalue_6h$'\t'\
		$pert_24h$'\t''NA'$'\t''NA'$'\t'$pvalue_24h$'\t'$adj_pvalue_24h$'\t'$pert_6h_24h\
		$'\t''NA'$'\t''NA'$'\t'$pvalue_6h_24h  >> ../bash_output/$METANALYSIS_FILENAME
	done
done