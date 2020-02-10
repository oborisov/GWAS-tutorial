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
# visualizing 2 first principal components
bfile=""
n_sd=6
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
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - median(eigenvec[[x]])) / IQR(eigenvec[[x]])+0.5)]
}
sd_iids=eigenvec[sd_for_PC1 > n_sd | sd_for_PC2 > n_sd]$IID
ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status, label = IID))+
geom_point() +
geom_label_repel(data=eigenvec[IID %in% sd_iids]) +
ggtitle(paste0("PC1=", eigenval[1], " PC2=", eigenval[2], " PC3=", eigenval[3], " PC4=", eigenval[4]))

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
rm ${bfile}_eigen.rm

%%bash
# if there are no pc outliers, copy files with "_pca" suffix
bfile=""
cp ${bfile}.bed ${bfile}_pca.bed
cp ${bfile}.bim ${bfile}_pca.bim
cp ${bfile}.fam ${bfile}_pca.fam

%%bash
### anchoring
bfile="/home/borisov/nsCLP/LKG_2018_Sharknado/LKG_2018_Sharknado_updids_updsex_checkedsex_geno02_mind02_geno002_mind002_norelated"

# pruning
plink --bfile ${bfile} \
--indep-pairwise 1000 50 0.2 \
--out ${bfile}_pruned
plink --bfile ${bfile} \
--extract <(cat ${bfile}_pruned.prune.in) \
--make-bed --out ${bfile}_pruned

# extracting pruned SNPs from 1000 Genomes
echo ${bfile}_pruned > ${bfile}_1kg_temp
for chr in {1..22}; do
plink --bfile /home/borisov/software/1000GP_Phase3/vcf/chr${chr} \
--extract <(awk '{print $2}' ${bfile}_pruned.bim) \
--make-bed --out ${bfile}_1kg_temp_chr${chr}
echo ${bfile}_1kg_temp_chr${chr} >> ${bfile}_1kg_temp
done

# merging 1000 genomes and data
plink --merge-list ${bfile}_1kg_temp --out ${bfile}_1kg_temp
plink --bfile ${bfile}_1kg_temp --geno 0 --make-bed --out ${bfile}_1kg_temp_geno0

# extracting PC from the merged data
salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 \
srun plink2 --bfile ${bfile}_1kg_temp_geno0 --pca --out ${bfile}_1kg_temp_geno0_eigen

# removing temp files
rm ${bfile}_1kg_temp

