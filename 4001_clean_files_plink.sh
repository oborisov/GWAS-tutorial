%%bash
# cleaning temporary plink files
bfile=""
bfile_short=$(echo $bfile | sed 's/_checkedsex.*//g')
rm ${bfile_short}*.bed ${bfile_short}*.bim \
${bfile_short}*.fam ${bfile_short}*.log ${bfile_short}*.hh
