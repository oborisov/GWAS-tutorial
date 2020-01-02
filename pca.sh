# performing pca with plink2
# visualizing 2 first principal components
plink_pca () {
  bfile=$1
  plink2 --bfile ${bfile} \
  --pca \
  --out ${bfile}
}

# identifying strong outliers - calculating percentiles of principal component 1 for every sample
echo '
library(ggplot2)
myargs=commandArgs[1]
eigenvec_file=fread(myargs)
p=ggplot(eigenvec_file, aes(x=PC1, y=PC2, label=IID))+
geom_point()
'

outliers_percentiles () {
  eigenvec_file=$1
  
}
```R
```

