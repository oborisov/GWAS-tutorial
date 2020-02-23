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
echo '#!/bin/bash
plink_regression () {
  bfile=$1
  plink2 --bfile ${bfile} \
  --pca \
  --out ${bfile}
  lambda_pc_iterate () {
  bfile=$1; covar_name=$2
  plink2 --bfile ${bfile} \
  --glm hide-covar cols=p,chrom,pos \
  --covar ${bfile}.eigenvec \
  --covar-col-nums ${covar_name} \
  --out ${bfile}_PC${covar_name}
  Rscript -e "median(qchisq(1 - (fread(commandArgs(T))[[5]]), 1))/qchisq(0.5,1)" \
  ${bfile}_PC${covar_name}.PHENO1.glm.logistic > \
  ${bfile}_PC${covar_name}.PHENO1.glm.logistic_lambda
  }
  for covar_name in 3 3-4 3-5 3-6 3-7 3-8 3-9 3-10 3-11 3-12; do
  lambda_pc_iterate ACE_HNR_full_QC ${covar_name} &
  done
  
  for covar_name in 3 3-4 3-5 3-6 3-7 3-8 3-9 3-10 3-11 3-12; do
  echo $covar_name
  cat ${bfile}_PC${covar_name}.PHENO1.glm.logistic_lambda
  done
}
plink_regression ${bfile}
' > ${bfile}_glm.sh
sbatch --wait --partition=medium --time=24:00:00 --cpus-per-task 16 --mem=16G ${bfile}_glm.sh
rm ${bfile}_glm.sh
