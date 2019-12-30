exclude_ambigous_AT_GC () {
  bfile=$1
  plink --bfile ${bfile} \
  --exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile}.bim) \
  --allow-no-sex \
  --make-bed --out ${bfile}_noambig
}

maf_geno_mind_hwe () {
  bfile=$1
  plink --bfile ${bfile} \
  --maf 0.01 --geno 0.2 --hwe 1e-6 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02
  plink --bfile ${bfile}_geno02 \
  --hwe 1e-10 include-nonctrl --mind 0.2 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02
  plink --bfile ${bfile}_geno02_mind02 \
  --geno 0.02 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002
  plink --bfile ${bfile}_geno02_mind02_geno002 \
  --mind 0.02 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002_mind002
}
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
plink_pca () {
  bfile=$1
  plink2 --bfile ${bfile} \
  --allow-no-sex \
  --pca --out ${bfile}
}
plink_regression () {
  bfile=$1
  plink --bfile ${bfile} \
  --logistic --ci 0.95 --allow-no-sex \
  --covar ${bfile}.eigenvec \
  --out ${bfile}
  head -1 ${bfile}.assoc.logistic > ${bfile}.assoc.logistic_pc
  awk '{if ($5 == "ADD") print $0}' ${bfile}.assoc.logistic >> ${bfile}.assoc.logistic_pc
  # qqman R script
  echo '
  options(bitmapType="cairo")
  library(qqman)
  library(QCEWAS)
  args=commandArgs(T)
  assoc_file=args[1]
  assoc=fread(assoc_file)
  assoc[CHR == "X", CHR := 23]
  assoc[, CHR := as.numeric(CHR)]
  assoc=assoc[!is.na(assoc$P)]
  assoc=assoc[P > 0]
  my_lambda=P_lambda(assoc$P)
  fwrite(data.table(filename=args, lambda=P_lambda(assoc$P)), paste0(args, "_lambda"), sep="\t")
  jpeg(paste0(assoc_file, "_manh.jpeg"), width = 12, height = 6, units = "in", res = 600)
  print(manhattan(assoc, chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T))
  title(main = assoc_file)
  dev.off()
  jpeg(paste0(assoc_file, "_qq.jpeg"), width = 6, height = 6, units = "in", res = 600)
  print(qq(assoc$P))
  title(main = assoc_file, sub = paste0("Lambda=", round(my_lambda,4)))
  dev.off()
  ' > ${bfile}.assoc.logistic.R
  Rscript ${bfile}.assoc.logistic.R ${bfile}.assoc.logistic_pc
}
qqman_R_script () {
  assoc_file_full_name=$1
  echo '
  options(bitmapType="cairo")
  library(qqman)
  library(QCEWAS)
  args=commandArgs(T)
  assoc_file=args[1]
  assoc=fread(assoc_file)
  assoc[CHR == "X", CHR := 23]
  assoc[, CHR := as.numeric(CHR)]
  assoc=assoc[!is.na(assoc$P)]
  assoc=assoc[P > 0]
  my_lambda=P_lambda(assoc$P)
  fwrite(data.table(filename=args, lambda=P_lambda(assoc$P)), paste0(args, "_lambda"), sep="\t")
  jpeg(paste0(assoc_file, "_manh.jpeg"), width = 12, height = 6, units = "in", res = 600)
  print(manhattan(assoc, chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T))
  title(main = assoc_file)
  dev.off()
  jpeg(paste0(assoc_file, "_qq.jpeg"), width = 6, height = 6, units = "in", res = 600)
  print(qq(assoc$P))
  title(main = assoc_file, sub = paste0("Lambda=", round(my_lambda,4)))
  dev.off()
  ' > ${assoc_file_full_name}_qqman.R
  Rscript ${assoc_file_full_name}_qqman.R ${assoc_file_full_name}
}
