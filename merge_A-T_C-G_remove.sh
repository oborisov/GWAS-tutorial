%%bash
# merge 2 data sets
bfile1=""
bfile2=""
bfile2_name=""
plink --bfile ${bfile1} \
--bmerge ${bfile2} \
--out ${bfile1}_${bfile2_name}
plink --bfile ${bfile1} \
--flip ${bfile1}_${bfile2_name}.missnp \
--make-bed --out ${bfile1}_filt
plink --bfile ${bfile1}_filt \
--bmerge ${bfile2} \
--out ${bfile1}_${bfile2_name}
# Remove ambigous variants after merging
bfile=""
plink --bfile ${bfile} \
--exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile}.bim) \
--allow-no-sex --make-bed --out ${bfile}_noambig
