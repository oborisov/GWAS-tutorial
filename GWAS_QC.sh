# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
# individual missingness (mind) excluding
# geno - SNP missingness
# maf - minor allele frequency (MAF)
# hwe - deviations from Hardy–Weinberg equilibrium (HWE)
# user can set maf, geno, mind by passing 2nd 3rd and 4th argument to maf_geno_mind_hwe
maf_geno_mind_hwe () {
  bfile=$1; maf_set=$2; geno_set=$3; mind_set=$4
  maf_set="${maf_set:-0.01}"
  geno_set="${geno_set:-0.02}"
  mind_set="${mind_set:-0.02}"
  echo $bfile $maf_set $geno_set $mind_set
  plink --bfile ${bfile} \
  --maf ${maf_set} \
  --geno 0.2 \
  --hwe 1e-6 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02
  plink --bfile ${bfile}_geno02 \
  --hwe 1e-10 include-nonctrl \
  --mind 0.2 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02
  plink --bfile ${bfile}_geno02_mind02 \
  --geno ${geno_set} \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002
  plink --bfile ${bfile}_geno02_mind02_geno002 \
  --mind ${mind_set} \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002_mind002
}




# relatedness
# ethnic outliers (see population stratification).
# inconsistencies in assigned and genetic sex of subjects (see sex discrepancy)

king_relatedness () {
  bfile=$1
  awk '{print $2,$2,$3,$4,$5,$6}' ${bfile}.fam > \
  ${bfile}1.fam
  cp ${bfile}.fam ${bfile}.famBK
  cp ${bfile}1.fam ${bfile}.fam
  king -b ${bfile}.bed \
  --related --degree 2 --prefix ${bfile} --cpus 22
  if [ -f ${bfile}.kin0 ]; then
    first_option_cases=$(grep -wFf <(awk '{if ($6 == 2) print $1,$2}' ${bfile}.fam | tr " " "_") \
    <(awk '{if ($14 != "UN") print $1,$2}' ${bfile}.kin0 | tr " " "_") | wc -l)
    second_option_cases=$(grep -wFf <(awk '{if ($6 == 2) print $1,$2}' ${bfile}.fam | tr " " "_") \
    <(awk '{if ($14 != "UN") print $3,$4}' ${bfile}.kin0 | tr " " "_") | wc -l)
    if [ ${second_option_cases} -ge ${first_option_cases} ]; then
      plink --bfile ${bfile} \
      --remove <(awk '{if ($14 != "UN") print $1,$2}' ${bfile}.kin0) \
      --allow-no-sex \
      --make-bed --out ${bfile}_norelated
    else
      plink --bfile ${bfile} \
      --remove <(awk '{if ($14 != "UN") print $3,$4}' ${bfile}.kin0) \
      --allow-no-sex \
      --make-bed --out ${bfile}_norelated
    fi
  else
    for suff in bed bim fam; do mv ${bfile}.${suff} ${bfile}_norelated.${suff}; done
  fi
  awk '{print 0,$2,$3,$4,$5,$6}' ${bfile}_norelated.fam > \
  ${bfile}_norelated1.fam
  mv ${bfile}_norelated1.fam ${bfile}_norelated.fam
}
plink_sex_check () {
  bfile=$1
  plink --bfile ${bfile} \
  --allow-no-sex \
  --check-sex --out ${bfile}
}
exclude_ambigous_AT_GC () {
  bfile=$1
  plink --bfile ${bfile} \
  --exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile}.bim) \
  --allow-no-sex \
  --make-bed --out ${bfile}_noambig
}
plink_pca () {
  bfile=$1
  plink2 --bfile ${bfile} \
  --allow-no-sex \
  --pca --out ${bfile}
}


# heterozygosity rate - to be added

