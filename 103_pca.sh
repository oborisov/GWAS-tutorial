# performing pca with plink2
# visualizing 2 first principal components

pruning_pca_fun () {
  bfile=$1
  plink --bfile ${bfile} \
  --indep-pairwise 1000 50 0.2 \
  --out ${bfile}_pruning >/dev/null 2>&1
  plink --bfile ${bfile} \
  --extract <(cat ${bfile}_pruning.prune.in) \
  --make-bed --out ${bfile}_pruned
  plink2 --bfile ${bfile}_pruned \
  --pca \
  --out ${bfile}_eigen
}

# plot PC
%%R
eigenvec=fread(".eigenvec")
eigenvec[, cc_status := "controls"]
eigenvec[grep("lkg",IID, ignore.case=T), cc_status := "cases"]
table(eigenvec$cc_status)
ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status))+
geom_point()



# identifying strong outliers - calculating percentiles of principal component 1 for every sample
echo '
library(ggplot2)
myargs=commandArgs[1]
eigenvec_file=fread(myargs)
p=ggplot(eigenvec_file, aes(x=PC1, y=PC2, label=IID))+
geom_point()
'
