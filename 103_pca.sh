# performing pca with plink2
%%bash
bfile=""
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
pruning_pca_fun ${bfile}
cat ${bfile}_eigen.eigenval

# visualizing 2 first principal components
%%R
eigenvec=fread(".eigenvec")
eigenvec[, cc_status := "controls"]
eigenvec[grep("lkg",IID, ignore.case=T), cc_status := "cases"]
table(eigenvec$cc_status)
ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status))+
geom_point()

# identifying strong outliers - calculating number of SD
%%R
for (x in 3:4) {
    ind=x+15
    mycol=paste0("sd_for_PC", x-2)
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - mean(eigenvec[[x]])) / sd(eigenvec[[x]]))]
}

# if there are pc outliers, copy files with "_pca" suffix
cp ${bfile}.bed ${bfile}_pca.bed
cp ${bfile}.bim ${bfile}_pca.bim
cp ${bfile}.fam ${bfile}_pca.fam
