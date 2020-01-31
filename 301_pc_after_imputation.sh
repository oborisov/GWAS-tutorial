# Merging chromosomes
%%bash
fam="_checkedsex_geno02_mind02_geno002_mind002_norelated_pca.fam"
sample="_geno02_mind02_geno002_mind002_norelated_pca_ref_chr1.phased.sample"
gen_prefix='_geno02_mind02_geno002_mind002_norelated_pca_ref_chr${chr}.phased.impute2'

# the following run automatically
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
plink --merge-list ${sample}_concat \
--out ${sample}_concat

# cleaning from temp plink per-chromosome files
for chr in {1..22}; do eval gen=${gen_prefix}
rm ${gen}.bed ${gen}.bim ${gen}.fam ${gen}.log
done; rm ${sample}_concat
