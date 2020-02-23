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

