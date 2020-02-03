%%bash
bfile=""
sex_check () {
  bfile=$1
  plink --bfile ${bfile} \
  --check-sex --out ${bfile}_sexcheck
}
sex_check ${bfile}

%%R
bfile=""
bfile=fread(paste0(bfile, "_sexcheck.sexcheck"))
bfile[STATUS != "OK"]

%%bash
bfile=""
plink --bfile ${bfile} \
--update-sex <(awk '{if ($5 == "PROBLEM") print $1,$2,$4}' ${bfile}_sexcheck.sexcheck) \
--make-bed --out ${bfile}_checkedsex

