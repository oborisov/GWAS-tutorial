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

%%R
# PC1&PC2 plot with centoid pc computed outliers
bfile=""
n_PC=4
n_sd=6
use_SD=T
if (use_SD) {center_fun <- mean; var_fun <- sd} else {center_fun <- median; var_fun <- IQR}
system(paste0("rm ", bfile, "_pruned*"))
library(ggrepel)
eigenvec=fread(paste0(bfile, "_eigen.eigenvec"))
eigenval=fread(paste0(bfile, "_eigen.eigenval"), header=F)
fam=fread(paste0(bfile, ".fam"), header=F)
colnames(fam)[6]="cc_status"
PC_means=data.table(t(sapply(eigenvec[,3:12], center_fun)))
eigenvec[, distance := apply(eigenvec[,3:(n_PC+2)], 1, function(vector_person_PCs) { apply(PC_means[,1:n_PC], 1, function(vector_PCs_means) { dist(rbind(vector_person_PCs, vector_PCs_means)) }) })]
eigenvec=merge(eigenvec, fam[,c(2,6)], by.x="IID", by.y="V2")
eigenvec[, cc_status := as.factor(cc_status)]
sd_iids=eigenvec[distance > n_sd*var_fun(distance)]$IID
print(ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status, label = IID))+
geom_point() + geom_label_repel(data=eigenvec[IID %in% sd_iids]) +
ggtitle(paste0("PC1=", eigenval[1], " PC2=", eigenval[2], " PC3=", eigenval[3], " PC4=", eigenval[4])))
fwrite(eigenvec[IID %in% sd_iids][,c(2,1)], paste0(bfile, "_eigen.rm"), col.names=F, sep=" ")
eigenvec[IID %in% sd_iids]

%%R
bfile=""
# PC1&PC2 plot with independent pc computed outliers
n_sd=6
use_SD=T
if (use_SD) {center_fun <- mean; var_fun <- sd} else {center_fun <- median; var_fun <- IQR}
system(paste0("rm ", bfile, "_pruned*"))
library(ggrepel)
eigenvec=fread(paste0(bfile, "_eigen.eigenvec"))
eigenval=fread(paste0(bfile, "_eigen.eigenval"), header=F)
fam=fread(paste0(bfile, ".fam"), header=F)
colnames(fam)[6]="cc_status"
eigenvec=merge(eigenvec, fam[,c(2,6)], by.x="IID", by.y="V2")
eigenvec[, cc_status := as.factor(cc_status)]
for (x in 3:4) {
    ind=x+15
    mycol=paste0("sd_for_PC", x-2)
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - center_fun(eigenvec[[x]])) / var_fun(eigenvec[[x]])+0.5)]
}
sd_iids=eigenvec[sd_for_PC1 > n_sd | sd_for_PC2 > n_sd]$IID
print(ggplot(eigenvec, aes(x=PC1, y=PC2, color=country, label = IID))+
geom_point() +
geom_label_repel(data=eigenvec[IID %in% sd_iids]) +
ggtitle(paste0("PC1=", eigenval[1], " PC2=", eigenval[2], " PC3=", eigenval[3], " PC4=", eigenval[4])))
# More than n_sd SD outliers based on PC1 and PC2:
fwrite(eigenvec[IID %in% sd_iids][,c(2,1)], paste0(bfile, "_eigen.rm"), col.names=F, sep=" ")
eigenvec[IID %in% sd_iids][,c(1,2,13,14,15)]


%%bash
# removing pc outliers (using _eigen.rm)
bfile=""
plink --bfile ${bfile} \
--remove <(grep -wFf <(awk '{print $2}' ${bfile}_eigen.rm) ${bfile}.fam) \
--make-bed --out ${bfile}_pca

%%bash
# if there are no pc outliers, copy files with "_pca" suffix
bfile=""
cp ${bfile}.bed ${bfile}_pca.bed
cp ${bfile}.bim ${bfile}_pca.bim
cp ${bfile}.fam ${bfile}_pca.fam
