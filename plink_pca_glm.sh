%%bash
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
