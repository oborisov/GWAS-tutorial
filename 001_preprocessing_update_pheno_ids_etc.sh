### --update-ids cannot be used in the same run as --update-parents or --update-sex.
%%bash
bfile="/home/borisov/nsCLP/LKG_2018_Sharknado/LKG_2018_Sharknado"

echo "" > ${bfile}_updpheno_temp
echo "" > ${bfile}_updsex_temp


plink --file  \
--make-bed --out ${bfile}

plink --bfile ${bfile} \
--pheno <(awk '{print $1,$2,2}' ${bfile}.fam) \
--update-ids ${bfile}_updpheno_temp \
--make-bed --out ${bfile}_updids

plink --bfile ${bfile}_updids \
--update-sex ${bfile}_updsex_temp \
--make-bed --out ${bfile}_updids_updsex

rm ${bfile}_updpheno_temp ${bfile}_updsex_temp ${bfile}_updids.*
