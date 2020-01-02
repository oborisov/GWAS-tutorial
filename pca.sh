# performing pca with plink2
# visualizing 2 first principal components
# identifying strong outliers
plink_pca () {
  bfile=$1
  plink2 --bfile ${bfile} \
  --allow-no-sex \
  --pca --out ${bfile}
}
