# performing pca with plink2
%%bash
bfile=""
pruning_pca_fun () {
  bfile=$1
  plink --bfile ${bfile} \
  --indep-pairwise 1000 50 0.2 \
  --out ${bfile}_pruning
  plink --bfile ${bfile} \
  --extract <(cat ${bfile}_pruning.prune.in) \
  --make-bed --out ${bfile}_pruned
  plink2 --bfile ${bfile}_pruned \
  --pca \
  --out ${bfile}_eigen
}
pruning_pca_fun ${bfile}
cat ${bfile}_eigen.eigenval

%%R
library(ggrepel)
# visualizing 2 first principal components
bfile=""
eigenvec=fread(paste0(bfile, "_eigen.eigenvec"))
eigenvec[, cc_status := "controls"]
eigenvec[grep("lkg",IID, ignore.case=T), cc_status := "cases"]
for (x in 3:4) {
    ind=x+15
    mycol=paste0("sd_for_PC", x-2)
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - mean(eigenvec[[x]])) / sd(eigenvec[[x]]))]
}
six_sd_iids=eigenvec[sd_for_PC1 > 6 | sd_for_PC2 > 6]$IID
ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status, label = IID))+
geom_point(color = ifelse(eigenvec$IID %in% six_sd_iids, "red", "grey50")) +
geom_label_repel(data=eigenvec[IID %in% six_sd_iids])

%%R
# More than 6 SD outliers based on PC1 and PC2:
eigenvec[IID %in% six_sd_iids]



# if there are no pc outliers, copy files with "_pca" suffix
cp ${bfile}.bed ${bfile}_pca.bed
cp ${bfile}.bim ${bfile}_pca.bim
cp ${bfile}.fam ${bfile}_pca.fam
# or remove pc outliers (using _eigen.rm)
plink --bfile ${bfile} \
--remove ${bfile}_eigen.rm \
--make-bed --out ${bfile}_pca
