exclude_ambigous_AT_GC () {
  bfile=$1
  plink --bfile ${bfile} \
  --exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile}.bim) \
  --allow-no-sex \
  --make-bed --out ${bfile}_noambig
}
