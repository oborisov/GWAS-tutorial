# from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
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
