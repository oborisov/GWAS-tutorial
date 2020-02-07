%%bash
# performing pca with plink2
bfile=""
plink --bfile ${bfile} \
--indep-pairwise 1000 50 0.2 \
--out ${bfile}_pruned
plink --bfile ${bfile} \
--extract <(cat ${bfile}_pruned.prune.in) \
--make-bed --out ${bfile}_pruned
salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 \
srun plink2 --bfile ${bfile}_pruned --pca --out ${bfile}_eigen
rm ${bfile}_pruned*
cat ${bfile}_eigen.eigenval

%%R
# visualizing 2 first principal components
bfile=""
n_sd=6
library(ggrepel)
eigenvec=fread(paste0(bfile, "_eigen.eigenvec"))
fam=fread(paste0(bfile, ".fam"), header=F)
colnames(fam)[6]="cc_status"
eigenvec=merge(eigenvec, fam[,c(2,6)], by.x="IID", by.y="V2")
eigenvec[, cc_status := as.factor(cc_status)]
for (x in 3:4) {
    ind=x+15
    mycol=paste0("sd_for_PC", x-2)
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - median(eigenvec[[x]])) / IQR(eigenvec[[x]])+0.5)]
}
sd_iids=eigenvec[sd_for_PC1 > n_sd | sd_for_PC2 > n_sd]$IID
ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status, label = IID))+
geom_point() + # color = ifelse(eigenvec$IID %in% sd_iids, "red", "grey50")
geom_label_repel(data=eigenvec[IID %in% sd_iids])


%%R
# More than n_sd SD outliers based on PC1 and PC2:
fwrite(eigenvec[IID %in% sd_iids][,c(2,1)], paste0(bfile, "_eigen.rm"), col.names=F, sep=" ")
eigenvec[IID %in% sd_iids][,c(1,2,13,14,15)]

%%bash
# removing pc outliers (using _eigen.rm)
bfile=""
plink --bfile ${bfile} \
--remove <(grep -wFf ${bfile}_eigen.rm ${bfile}.fam) \
--make-bed --out ${bfile}_pca


%%bash
# if there are no pc outliers, copy files with "_pca" suffix
bfile=""
cp ${bfile}.bed ${bfile}_pca.bed
cp ${bfile}.bim ${bfile}_pca.bim
cp ${bfile}.fam ${bfile}_pca.fam
