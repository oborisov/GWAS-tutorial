%%bash
# Merging chromosomes
bfile=""

# the following run automatically
fam=${bfile}.fam
sample=$(echo ${bfile}_ref_chr1.phased.sample | sed 's/_chr1//')
gen_prefix=${bfile}'_ref_chr${chr}.phased.impute2'
rm ${sample}_concat 2> /dev/null
for chr in {1..22}; do
eval gen=${gen_prefix}
plink2 --gen ${gen}.gz --sample ${bfile}_ref_chr1.phased.sample \
--oxford-single-chr ${chr} \
--update-sex <(paste -d ' ' <(awk '{print $1,$2}' $fam | tr " " "_" | awk '{print $1,$1}') <(awk '{print $5}' $fam)) \
--pheno <(paste -d ' ' <(awk '{print $1,$2}' $fam | tr " " "_" | awk '{print $1,$1}') <(awk '{print $6}' $fam)) \
--extract <(zcat ${gen}_info.gz | awk '{if ($7 > 0.8) print $2}') \
--make-bed --out ${gen} > /dev/null 2>&1
echo ${gen} >> ${sample}_concat
done


salloc --mem=32000M --time=5:00:00 --cpus-per-task=1 \
srun plink --merge-list ${sample}_concat \
--out ${sample}_concat
#######################################################

%%bash
bfile=""

# the following run automatically
fam=${bfile}.fam
sample=$(echo ${bfile}_ref_chr1.phased.sample | sed 's/_chr1//')
gen_prefix=${bfile}'_ref_chr${chr}.phased.impute2'
# cleaning from temp plink per-chromosome files
for chr in {1..22}; do eval gen=${gen_prefix}
eval gen=${gen_prefix}
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
#######################################################

%%R
bfile=""
bfile=paste0(bfile, "_ref.phased.sample_concat")
# PC1&PC2 plot with independent pc computed outliers
n_PC=2
n_sd=6
use_SD=T
if (use_SD) {center_fun <- mean; var_fun <- sd} else {center_fun <- median; var_fun <- IQR}

# library
library(ggrepel)

# fam, eigenvec, eigenval
eigenvec=fread(paste0(bfile, "_eigen.eigenvec"))
eigenval=fread(paste0(bfile, "_eigen.eigenval"), header=F)
fam=fread(paste0(bfile, ".fam"), header=F)
colnames(fam)[6]="cc_status"
eigenvec=merge(eigenvec, fam[,c(2,6)], by.x="IID", by.y="V2")
eigenvec[, cc_status := as.factor(cc_status)]
for (x in 3:(2+n_PC)) { mycol=paste0("sd_for_PC", x-2)
    eigenvec[, (mycol) := round(abs(eigenvec[[x]] - center_fun(eigenvec[[x]])) / var_fun(eigenvec[[x]])+0.5)] }
column_number_vector=c(14:(13+n_PC))
eigenvec_melt = melt(eigenvec[,c(1,13:ncol(eigenvec)), with=F], id.vars = c("IID"), 
                     measure.vars = colnames(eigenvec[,c(14:ncol(eigenvec)), with=F]))
print(eigenvec_melt[value > n_sd])

print(ggplot(eigenvec, aes(x=PC1, y=PC2, color=cc_status, label = IID))+
geom_point() +
geom_label_repel(data=eigenvec[IID %in% eigenvec_melt[value > n_sd]$IID]) + theme_bw())
ggsave(paste0(bfile, "_eigen.eigenvec_pca.jpeg"))

# More than n_sd SD outliers based on PC1 and PC2:
fwrite(eigenvec[IID %in% eigenvec_melt[value > n_sd]$IID][,c(2,1)], paste0(bfile, "_eigen.rm"), col.names=F, sep=" ")

# Eigenvalues
colnames(eigenval)="Eigenvalues"
eigenval$PC=factor(paste0("PC", 1:10), levels=paste0("PC", 1:10))
eigenval$PC_temp=as.integer(1:10)
eigenval[, perc_var_explained := round(eigenval[[1]]*100/sum(eigenval[[1]]),1)]

# Variance explained by PCs
print(eigenval[,c(2,1,4)])

# Scree plot
print(ggplot(eigenval, aes(x=PC, y=Eigenvalues)) +
geom_point() + geom_line(aes(x=PC_temp, y=Eigenvalues)) + theme_bw() +
ggtitle(paste0("PC1=", eigenval[1,4], "%; PC2=", eigenval[2,4], "%; PC3=", eigenval[3,4], "%; PC4=", eigenval[4,4], "%")))
ggsave(paste0(bfile, "_eigen.eigenval_scree.jpeg"))

# remove temp files
system(paste0("rm ", bfile, "_pruned*"))
