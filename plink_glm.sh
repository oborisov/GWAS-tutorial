## Regression

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

### QQ

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

### Iterate PC
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

### QQ
options(bitmapType="cairo")
library(qqman)
library(RColorBrewer)
sumstats=fread("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic")
colnames(sumstats)=c("CHR", "BP", "SNP", "A1", "P")
sumstats[CHR == "X", CHR := 23]
sumstats[, CHR := as.numeric(CHR)]
sumstats=sumstats[!is.na(CHR)]
sumstats=sumstats[!is.na(P)]
sumstats=sumstats[P > 0]
jpeg("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic_manh.jpeg", width = 12, height = 6, units = "in", res = 200)
print(manhattan(sumstats, chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T, col=brewer.pal(8, "Dark2")))
title(main = "ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic")
dev.off()
my_lambda=median(qchisq(1 - (sumstats[[5]]), 1))/qchisq(0.5,1)
jpeg("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic_qq.jpeg", width = 6, height = 6, units = "in", res = 200)
print(qq(sumstats[[5]]))
title(main = "ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic", sub = paste0("Lambda=", round(my_lambda,4)))
dev.off()

### Lambda
%%bash
### Computing lambda
bfile=""
# plink logistic with 4 pc
bfile="${bfile}_ref.phased.sample_concat"
plink --bfile ${bfile} \
--logistic \
--covar ${bfile}_eigen.eigenvec \
--covar-name PC1, PC2, PC3, PC4 \
--out ${bfile}_lambda
# select p values for phenotype only
awk '{if ($5 == "ADD") print $9}' ${bfile}_lambda.assoc.logistic > ${bfile}_lambda_P
# remove temp files
rm ${bfile}_lambda.assoc.logistic
# Computing lambda
Rscript -e 'bfile=commandArgs(T); dt=data.table(lambda_4pc=median(qchisq(1 - fread(paste0(bfile, "_lambda_P"), header=F)[[1]], 1), na.rm=T) / qchisq(0.5, 1)); system(paste0("rm ", paste0(bfile, "_lambda_P"))); fwrite(dt, paste0(bfile, "_lambda")); dt' ${bfile} > ${bfile}_lambda
#####################################################
