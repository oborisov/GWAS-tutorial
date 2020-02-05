%%bash
# checking sex of samples using X-chromosome
bfile=""
plink --bfile ${bfile} --indep-pairwise 1000 50 0.2 --out ${bfile}_pruning
plink --bfile ${bfile} --extract <(cat ${bfile}_pruning.prune.in) \
--make-bed --out ${bfile}_pruning
plink --bfile ${bfile}_pruning --check-sex --out ${bfile}_sexcheck
rm ${bfile}_pruning*

%%R
# listing samples with incorrect sex assignment
bfile=""
bfile=fread(paste0(bfile, "_sexcheck.sexcheck"))
bfile[STATUS != "OK"]

%%bash
# Updating sex using SNPSEX (determined by X-chromosome)
bfile=""
plink --bfile ${bfile} \
--update-sex <(awk '{if ($5 == "PROBLEM") print $1,$2,$4}' ${bfile}_sexcheck.sexcheck) \
--make-bed --out ${bfile}_checkedsex

