%%bash
bfile="/home/borisov/LUTO/Dutch/dcg_checkedsex_geno02_mind02_geno002_mind002_norelated_pca_pca_pca_pca_pca"
plink --bfile ${bfile} \
--logistic \
--out ${bfile}_test_nopc
plink --bfile ${bfile} \
--covar ${bfile}_eigen.eigenvec \
--covar-name PC1, PC2, PC3, PC4 \
--logistic \
--out ${bfile}_test_4PC

%%R
bfile="/home/borisov/LUTO/Dutch/ioe24_hoe24_NBS_updatedsexpheno_onlywithpheno_SNPSEX_checkedsex_geno02_mind02_geno002_mind002_norelated_pca_pca_pca_pca_pca"
suffix="_test_nopc"
pvec=fread(paste0(bfile, suffix, ".assoc.logistic"))
median(qchisq(1 - pvec[[9]], 1), na.rm=T) / qchisq(0.5, 1)

%%R
pvec_nopc[order(P)][1:5]
%%R
pvec_4pc[order(P)][1:5]
