%%bash
bfile=""
plink --bfile ${bfile} \
--logistic \
--out ${bfile}_test_nopc
plink --bfile ${bfile} \
--covar ${bfile}_eigen.eigenvec \
--covar-name PC1, PC2, PC3, PC4 \
--logistic \
--out ${bfile}_test_4pc

%%R
bfile=""
pvec_nopc=fread(paste0(bfile, "_test_nopc", ".assoc.logistic"))
pvec_4pc=fread(paste0(bfile, "_test_4pc", ".assoc.logistic"))
pvec_4pc=pvec_4pc[TEST=="ADD"]
print(data.table(lambda_nopc=median(qchisq(1 - pvec_nopc[[9]], 1), na.rm=T) / qchisq(0.5, 1),
          lambda_4pc=median(qchisq(1 - pvec_4pc[[9]], 1), na.rm=T) / qchisq(0.5, 1)))
print(pvec_nopc[order(P)][1:10])
print(pvec_4pc[order(P)][1:10])
