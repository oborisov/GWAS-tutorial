%%bash
sex_check () {
  bfile=$1
  plink --bfile ${bfile} \
  --check-sex --out ${bfile}_sexcheck
}
sex_check >/dev/null 2>&1
