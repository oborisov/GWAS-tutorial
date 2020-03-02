%%bash
# check if X is split into 23 and 25 (pseudoautosomal region), split if needed
bfile=""
plink --bfile ${bfile} --split-x hg19 --make-bed --out ${bfile}_splitX

%%bash
# checking sex of samples using X-chromosome
bfile="_splitX"
plink --bfile ${bfile} --chr 23 --indep-pairwise 1000 50 0.2 --out ${bfile}_pruning
plink --bfile ${bfile} --extract <(cat ${bfile}_pruning.prune.in) \
--make-bed --out ${bfile}_pruning
plink --bfile ${bfile}_pruning --check-sex --out ${bfile}_sexcheck
rm ${bfile}_pruning*
Rscript -e 'fread(paste0(commandArgs(T), "_sexcheck.sexcheck"))[STATUS != "OK"]' ${bfile}

%%bash
bfile=""
# Updating sex using SNPSEX (determined by X-chromosome)
# Removing SNPSEX=0
plink --bfile ${bfile} \
--update-sex <(awk '{if ($5 == "PROBLEM") print $1,$2,$4}' ${bfile}_sexcheck.sexcheck) \
--remove <(awk '{if ($4 == 0) print $1,$2}' ${bfile}_sexcheck.sexcheck) \
--make-bed --out ${bfile}_checkedsex

