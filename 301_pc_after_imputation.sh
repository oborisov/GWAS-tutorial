%%bash
# Merging chromosomes
fam="_checkedsex_geno02_mind02_geno002_mind002_norelated_pca.fam"
sample="_geno02_mind02_geno002_mind002_norelated_pca_ref_chr1.phased.sample"
gen_prefix='_geno02_mind02_geno002_mind002_norelated_pca_ref_chr${chr}.phased.impute2'

# the following run automatically
sample_out=$(echo ${sample} | sed 's/_chr1//')
rm ${sample}_concat 2> /dev/null
for chr in {1..22}; do
eval gen=${gen_prefix}
plink2 --gen ${gen}.gz --sample ${sample} \
--oxford-single-chr ${chr} \
--update-sex <(paste -d ' ' <(awk '{print $1,$2}' $fam | tr " " "_" | awk '{print $1,$1}') <(awk '{print $5}' $fam)) \
--pheno <(paste -d ' ' <(awk '{print $1,$2}' $fam | tr " " "_" | awk '{print $1,$1}') <(awk '{print $6}' $fam)) \
--make-bed --out ${gen} &
echo ${gen} >> ${sample}_concat
done
wait

# merging all chromosomes into a single file
salloc --mem=32000M --time=5:00:00 \
srun plink --merge-list ${sample}_concat \
--out ${sample}_concat

# cleaning from temp plink per-chromosome files
for chr in {1..22}; do eval gen=${gen_prefix}
rm ${gen}.bed ${gen}.bim ${gen}.fam ${gen}.log
done; rm ${sample}_concat

bfile="${sample}_concat"
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
library(ggrepel)
n_sd=6
bfile=""
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

