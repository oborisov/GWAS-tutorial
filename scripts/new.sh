%%bash
bfile="/home/borisov/ARM/Dutch_cases"
#### checking genomic build
#http://ncbi.nlm.nih.gov/snp/rs1000002  
#GRCh37.p13 chr 3 NC_000003.11:g.183635768C>T  
#GRCh38.p12 chr 3 NC_000003.12:g.183917980C>T  
grep -wF rs1000002 ${bfile}.bim

%%bash
bfile="/home/borisov/ARM/Dutch_cases"
# checking sex of samples using X-chromosome
# http://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
bfile_temp=${bfile}
if plink --bfile ${bfile} --split-x hg19 --make-bed --out ${bfile}_splitX; then bfile_temp=${bfile}_splitX; else echo -e "\n$(tput bold)Dataset already contains an XY region$(tput sgr0)\n"; fi
plink --bfile ${bfile_temp} --chr 23 --maf 0.01 --extract <(grep rs ${bfile_temp}.bim) \
--indep-pairphase 20000 2000 0.5 --out ${bfile_temp}_pruned
plink --bfile ${bfile_temp} --extract ${bfile_temp}_pruned.prune.in \
--check-sex 0.2 0.8 --out ${bfile}_checkedsex_report
awk '{if ($5 != "OK") print $0}' ${bfile}_checkedsex_report.sexcheck
plink --bfile ${bfile_temp} \
--remove <(awk '{if ($4 == 0) print $1,$2}' ${bfile}_checkedsex_report.sexcheck) \
--update-sex <(awk '{if ($5 == "PROBLEM") print $1,$2,$4}' ${bfile}_checkedsex_report.sexcheck) \
--chr 1-23,25 --set-hh-missing --make-bed --out ${bfile}_checkedsex
if [ -f ${bfile}_splitX.bed ]; then rm ${bfile}_splitX*; fi
Rscript -e 'dat=fread(commandArgs(T)); print(summary(dat[F < 0.5]$F)); print(summary(dat[F > 0.5]$F))' ${bfile}_checkedsex_report.sexcheck

%%bash
## Standard QC for genotyping data: maf, missingness of samples and variants, hwe
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
bfile="/home/borisov/ARM/Dutch_cases_checkedsex"
maf_set=0.001
geno_set=0.02
mind_set=0.02
hwe_set=0.0001
echo $bfile $maf_set $geno_set $mind_set
plink --bfile ${bfile} \
--maf ${maf_set} \
--geno 0.2 \
--hwe 1e-6 \
--allow-no-sex \
--make-bed --out ${bfile}_geno02
plink --bfile ${bfile}_geno02 \
--hwe ${hwe_set} include-nonctrl \
--mind 0.2 \
--allow-no-sex \
--make-bed --out ${bfile}_geno02_mind02
plink --bfile ${bfile}_geno02_mind02 \
--geno ${geno_set} \
--allow-no-sex \
--make-bed --out ${bfile}_geno02_mind02_geno002
plink --bfile ${bfile}_geno02_mind02_geno002 \
--mind ${mind_set} \
--allow-no-sex \
--make-bed --out ${bfile}_QC
rm ${bfile}_geno02.b* ${bfile}_geno02_mind02.b* ${bfile}_geno02_mind02_geno002.b*
if [ -f ${bfile}_geno02_mind02.irem ]; then cat ${bfile}_geno02_mind02.irem; fi
if [ -f ${bfile}_QC.irem ]; then cat ${bfile}_QC.irem; fi

%%bash
bfile1="/home/borisov/ARM/Dutch_cases_checkedsex_QC"
bfile2="/home/borisov/ARM/Dutch_controls_checkedsex_QC"
bfile_merged="/home/borisov/ARM/Dutch_cases_controls"
plink --bfile ${bfile1} --bmerge ${bfile2} --out ${bfile_merged} 2> /dev/null
if [ -f ${bfile_merged}.missnp ]; then
plink --bfile ${bfile1} --flip ${bfile_merged}.missnp --make-bed --out ${bfile1}_filt
plink --bfile ${bfile1}_filt --bmerge ${bfile2} --out ${bfile_merged} 2> /dev/null
fi
plink --bfile ${bfile_merged} \
--exclude <(awk '{if ($5 == "A" && $6 == "T" || $5 == "T" && $6 == "A" || $5 == "C" && $6 == "G" || $5 == "G" && $6 == "C") print $2}' ${bfile_merged}.bim) \
--make-bed --out ${bfile_merged}_noambig

%%bash
## Standard QC for genotyping data: maf, missingness of samples and variants, hwe
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
bfile="/home/borisov/ARM/Dutch_cases_checkedsex"
maf_set=0.001
geno_set=0.02
mind_set=0.02
hwe_set=0.0001
echo $bfile $maf_set $geno_set $mind_set
plink --bfile ${bfile} \
--maf ${maf_set} \
--geno 0.2 \
--hwe 1e-6 \
--allow-no-sex \
--make-bed --out ${bfile}_geno02
plink --bfile ${bfile}_geno02 \
--hwe ${hwe_set} include-nonctrl \
--mind 0.2 \
--allow-no-sex \
--make-bed --out ${bfile}_geno02_mind02
plink --bfile ${bfile}_geno02_mind02 \
--geno ${geno_set} \
--allow-no-sex \
--make-bed --out ${bfile}_geno02_mind02_geno002
plink --bfile ${bfile}_geno02_mind02_geno002 \
--mind ${mind_set} \
--allow-no-sex \
--make-bed --out ${bfile}_QC
rm ${bfile}_geno02.b* ${bfile}_geno02_mind02.b* ${bfile}_geno02_mind02_geno002.b*
if [ -f ${bfile}_geno02_mind02.irem ]; then cat ${bfile}_geno02_mind02.irem; fi
if [ -f ${bfile}_QC.irem ]; then cat ${bfile}_QC.irem; fi


%%bash
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC"
degree=3
# Checking relatedness using king
# from http://people.virginia.edu/~wc9c/KING/manual.html
# updating fids
awk '{print $1,$2,$2,$2}' ${bfile}.fam > ${bfile}.fam_temp
plink --bfile ${bfile} --update-ids ${bfile}.fam_temp --make-bed --out ${bfile}_updids; rm ${bfile}.fam_temp
# king
#salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 srun \
king -b ${bfile}_updids.bed --kinship --degree ${degree} --prefix ${bfile}_updids # --cpus 20
# removing relatives
if [ $degree == 3 ]; then kinship=0.0442; fi
if [ $degree == 2 ]; then kinship=0.0884; fi
rm ${bfile}_updids.related_excluded 2> /dev/null
tail -n +2 ${bfile}_updids.kin0 | awk -v kinship=${kinship} '{if ($NF > kinship) print $0}' | \
while read l; do if [[ $(grep -wFf <(echo ${l} | awk '{print $1,$2}') ${bfile}_updids.fam | awk '{print $6}') -eq 1 ]]; then \
echo ${l} | awk '{print $1,$2}' >> ${bfile}_updids.related_excluded; \
else echo ${l} | awk '{print $3,$4}' >> ${bfile}_updids.related_excluded; fi; done
plink --bfile ${bfile}_updids --remove ${bfile}_updids.related_excluded --make-bed --out ${bfile}_norelated
echo "Removed sample(s):"
grep -wFf ${bfile}_updids.related_excluded ${bfile}_updids.fam
# temoving temp files
rm ${bfile}_updids*

%%bash
# anchoring to 1kG
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated"
# extracting SNPs from 1000 Genomes
echo ${bfile} > ${bfile}_1kg_temp
for chr in {1..22}; do
plink --bfile /home/borisov/software/1000GP_Phase3/vcf/chr${chr} --extract <(awk '{print $2}' ${bfile}.bim) --make-bed --out ${bfile}_1kg_temp_chr${chr} > /dev/null 2>&1
echo ${bfile}_1kg_temp_chr${chr} >> ${bfile}_1kg_temp
done
# merging 1000 genomes and data
plink --merge-list ${bfile}_1kg_temp --out ${bfile}_1kg_temp 2> /dev/null
plink --bfile ${bfile} --exclude ${bfile}_1kg_temp.missnp --make-bed --out ${bfile}_1kg_temp_filt
suffix=$(echo $bfile | sed 's/.*\///'); cat ${bfile}_1kg_temp | sed "s/${suffix}$/${suffix}_1kg_temp_filt/" > ${bfile}_1kg_temp2
plink --merge-list ${bfile}_1kg_temp2 --out ${bfile}_1kg_temp 2> /dev/null
# filter geno 0.001, prunning
plink2 --bfile ${bfile}_1kg_temp --geno 0.001 --indep-pairwise 1000 50 0.2 --out ${bfile}_1kg_temp_pruned
# pca in plink2 #salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 srun \
plink2 --bfile ${bfile}_1kg_temp --extract ${bfile}_1kg_temp_pruned.prune.in --pca --out ${bfile}_1kg_all_eigen
# removing temp files
rm ${bfile}*1kg_temp*

%%R
### Visualizing
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated"
kg_sample_file="/home/borisov/software/1000GP_Phase3/1000GP_Phase3.sample"
n_PC=4; bfile_population="EUR"
# library
library(ggrepel)
# sample, fam, eigenvec, eigenval
kg_sample=fread(kg_sample_file)
fam=fread(paste0(bfile, ".fam"))
eigenvec=fread(paste0(bfile, "_1kg_all_eigen.eigenvec"))
eigenval=fread(paste0(bfile, "_1kg_all_eigen.eigenval"))
# Processing
eigenvec=merge(eigenvec, kg_sample, by.x="IID", by.y="sample", all.x=T)
eigenvec[is.na(population), population := "samples"]
eigenvec[is.na(group), group := "samples"]
PC_means=eigenvec[group != "samples"][,list(PC1_mean=mean(PC1), PC2_mean=mean(PC2), PC3_mean=mean(PC3),PC4_mean=mean(PC4), PC5_mean=mean(PC5), PC6_mean=mean(PC6), PC7_mean=mean(PC7),PC8_mean=mean(PC8), PC9_mean=mean(PC9),PC10_mean=mean(PC10)),by=group]
distances_dt=data.table(IID=eigenvec[group == "samples"]$IID, t(apply(eigenvec[group == "samples"][,3:(n_PC+2)], 1, function(person) { apply(PC_means[,2:(n_PC+1)],1, function(x) { dist(rbind(person,x))  })})))
colnames(distances_dt) = c("IID", PC_means$group)
distances_dt$PC_group = PC_means$group[apply(distances_dt[,2:6], 1, which.min)]
# Print outliers IDs to stdout
print(distances_dt[PC_group != bfile_population][,c(1,ncol(distances_dt)), with=F])
fwrite(eigenvec[IID %in% distances_dt[PC_group != bfile_population]$IID, c(2,1)], paste0(bfile, "_1kG.rm"), sep=" ", col.names=F, quote=F, na=NA)
# PCA plot. Highlight outliers IDs
print(ggplot(eigenvec, aes(x=PC1, y=PC2, color=group, label=IID)) +
geom_point() + geom_label_repel(data=eigenvec[IID %in% distances_dt[PC_group != bfile_population]$IID]))
ggsave(paste0(bfile, "_1kg_all_eigen.eigenvec_pca.jpeg"))
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
ggsave(paste0(bfile, "_1kg_all_eigen.eigenval_scree.jpeg"))
# remove temp files
system(paste0("rm ", bfile, "_1kg_temp*"))

%%bash
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated"
# remove 1kg outliers
plink --bfile ${bfile} --remove ${bfile}_1kG.rm --make-bed --out ${bfile}_1kG
# pca
bfile="${bfile}_1kG"
plink2 --bfile ${bfile} --indep-pairwise 1000 50 0.2 --out ${bfile}_pruned > /dev/null 2>&1
#salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 srun \
plink2 --bfile ${bfile} --extract ${bfile}_pruned.prune.in --pca --out ${bfile}_eigen

%%R
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated_1kG"
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

%%bash
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated_1kG"
# remove 1kg outliers
plink --bfile ${bfile} --remove ${bfile}_eigen.rm --make-bed --out ${bfile}_pca
# pca
bfile="${bfile}_pca"
plink2 --bfile ${bfile} --indep-pairwise 1000 50 0.2 --out ${bfile}_pruned > /dev/null 2>&1
#salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 srun \
plink2 --bfile ${bfile} --extract ${bfile}_pruned.prune.in --pca --out ${bfile}_eigen
# glm
plink2 --bfile ${bfile} \
--covar ${bfile}_eigen.eigenvec --covar-name PC1,PC2,PC3,PC4 \
--glm hide-covar --out ${bfile}_test_4pc

%%R
bfile="/home/borisov/ARM/Dutch_cases_controls_noambig_QC_norelated_1kG_pca"
lapply(c("qqman", "RColorBrewer"), require, character.only = TRUE)
pvec_4pc=fread(paste0(bfile, "_test_4pc.PHENO1.glm.logistic"))
print(pvec_4pc[order(P)][1:10])
# adjusting CHR and p
pvec_4pc[`#CHROM` == "X" | `#CHROM` == "XY", `#CHROM` := 23]
pvec_4pc[, `#CHROM` := as.numeric(`#CHROM`)]
pvec_4pc=pvec_4pc[!is.na(pvec_4pc$`#CHROM`) & !is.na(pvec_4pc$P) & P > 0]
lambda=median(qchisq(1 - pvec_4pc$P, 1), na.rm=T) / qchisq(0.5, 1)
# producing Manhattan plot
print(manhattan(rbind(pvec_4pc[P<5e-2], pvec_4pc[P>5e-2][seq(1,nrow(pvec_4pc[P>5e-2]),10)]), chr="#CHROM", bp="POS", p="P", snp="SNP", annotatePval = 1, annotateTop = T, col=brewer.pal(8, "Dark2")))
# producing Q-Q plot
print(qq(pvec_4pc$P))
title(main = paste0("Lambda=", round(lambda,4)))


%%R
bfile=commandArgs(T)
# libraries
lapply(c("qqman", "RColorBrewer"), require, character.only = TRUE)
# fread summary statistics
sumstats=rbindlist(sapply(list.files(path=gsub("(.*)/.*", "\\1", bfile), pattern=paste0(gsub(".*/(.*)", "\\1", bfile), ".*.glm.logistic$"), full.names=T), fread, simplify=F))
colnames(sumstats)[1:3]=c("CHR", "BP", "SNP")
# adjusting CHR and p
sumstats[CHR == "X", CHR := "23"]
sumstats[CHR == "XY", CHR := "25"]
sumstats[, CHR := as.numeric(CHR)]
# adjusting SNP
sumstats[, SNP := gsub(".*(rs[0-9]*).*", "\\1", SNP)]
sumstats=sumstats[!is.na(sumstats$P) & P > 0]
# calculating lambda https://www.biostars.org/p/298847/
lambda <- median(qchisq(1 - sumstats$P, 1), na.rm=T) / qchisq(0.5, 1)
# producing Manhattan plot
jpeg(paste0(bfile, "_imputed_manh.jpeg"), width = 12, height = 6, units = "in", res = 600)
manhattan_obj <- manhattan(rbind(sumstats[P<5e-2], sumstats[P>5e-2][seq(1,nrow(sumstats[P>5e-2]),10)]), chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T, col=brewer.pal(8, "Dark2"))
title(main = paste0(bfile, "_imputed_manhattan.jpeg"))
dev.off()
saveRDS(manhattan_obj, file = paste0(bfile, "_imputed_manhattan.rds"))
# producing Q-Q plot
jpeg(paste0(bfile, "_imputed_qq.jpeg"), width = 12, height = 6, units = "in", res = 600)
qq_obj <- qq(sumstats$P)
print(qq_obj)
title(main = paste0(bfile, "_imputed_manhattan.jpeg"), sub = paste0("Lambda=", round(lambda,4)))
dev.off()
saveRDS(qq_obj, file = paste0(bfile, "_imputed_qq.rds"))
# writing suggestive associations into a file
fwrite(sumstats[P < 1e-5], paste0(bfile, "_imputed_suggestive.tsv"), sep="\t", na=NA, quote=F)
fwrite(sumstats, paste0(bfile, "_imputed_all.tsv"), sep="\t", na=NA, quote=F)
system(paste0("gzip -f ", bfile, "_imputed_suggestive.tsv ", bfile, "_imputed_all.tsv"))


