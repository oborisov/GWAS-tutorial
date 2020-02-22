%%bash
# merge 2 data sets
bfile1=""
bfile2=""
plink --bfile ${bfile1} \
--bmerge ${bfile2} \
--out ${bfile1}_bfile2
plink --bfile ${bfile1} \
--flip ${bfile1}_bfile2.missnp \
--make-bed --out ${bfile1}_filt
plink --bfile ${bfile1}_filt \
--bmerge ${bfile2}_Bonn \
--out {bfile1}_bfile2
# Remove ambigous variants after merging
bfile=""
plink --bfile ${bfile} \
--exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile}.bim) \
--allow-no-sex --make-bed --out ${bfile}_noambig
