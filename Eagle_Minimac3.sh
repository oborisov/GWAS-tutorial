%%bash
#convert data to reference strand to impute sequencing panel
# Check that FIDs or IIDs do not have underscores "_"
bfile="/home/borisov/BOADICEA/GSP_2020-041-ILL_DIA_N_1/GSP_2020-041-ILL_DIA_N_1_splitX_HNR_geno02_mind02_geno002_mind002"
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
# the following run automatically
rm ${bfile}_jobs.list >/dev/null 2>&1
echo '#!/bin/bash
bfile=$1
reference_fasta=$2
export BCFTOOLS_PLUGINS=${3}
chr=$4
echo "fixref chromosome ${chr}"
plink --bfile ${bfile} --chr ${chr} --recode vcf --out ${bfile}_chr${chr}
bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
' > ${bfile}_fixref.sh
for chr in {1..22}; do
sbatch --parsable --mem=16G --time=05:00:00 --cpus-per-task=8 \
    ${bfile}_fixref.sh \
    ${bfile} \
    ${reference_fasta} \
    ${bcftools_plugins} \
    ${chr} >> ${bfile}_jobs.list
done
squeue | grep -wFf ${bfile}_jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done
# removing input vcfs (temrorary files)
for chr in {1..22}; do rm ${bfile}_chr${chr}.vcf ${bfile}_chr${chr}.log; done
rm ${bfile}_fixref.sh
echo "fixref is done"

%%bash
# ChrX
bfile="/home/borisov/BOADICEA/GSP_2020-041-ILL_DIA_N_1/GSP_2020-041-ILL_DIA_N_1_splitX_HNR_geno02_mind02_geno002_mind002"
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
n_threads=12

## Preprocessing
# remove males haploid heterozygous calls (should be done when X is splitted into PAR and nonPAR)
plink --bfile ${bfile} --set-hh-missing --make-bed --out ${bfile}_hh_missing

# to vcf
plink --bfile ${bfile}_hh_missing --chr 25 --recode vcf --out ${bfile}_chr25
plink --bfile ${bfile}_hh_missing --chr 23 --filter-females --recode vcf --out ${bfile}_chr23_females
plink --bfile ${bfile}_hh_missing --chr 23 --filter-males --recode vcf --out ${bfile}_chr23_males

# fixref
bcftools annotate --rename-chrs <(echo "25 X") ${bfile}_chr25.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_PAR.vcf.gz
bcftools annotate --rename-chrs <(echo "23 X") ${bfile}_chr23_females.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz
bcftools annotate --rename-chrs <(echo "23 X") ${bfile}_chr23_males.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz
bcftools index ${bfile}_ref_chrX_PAR.vcf.gz; bcftools index ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz; bcftools index ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz

### Phasing
# Phasing PAR
sbatch --job-name eagle_BOADICEA_X \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chrX_PAR.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr23.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom 23 \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chrX_PAR_phased" > /dev/null
# Phasing Female.Non.PAR
sbatch --job-name eagle_BOADICEA_X \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr23.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom 23 \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chrX_Female.Non.PAR_phased" > /dev/null
# Male.Non.PAR do not need to be phased as they are hemizygous

while squeue | grep eagle_BO > /dev/null; do sleep 5; done

### Imputation 
# Imputing PAR
sbatch --job-name Minimac3_BOADICEA_X \
--cpus-per-task=1 --mem=32G --partition=long --time=720:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_PAR_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_PAR_phased_imputed" > /dev/null

# Female.Non.PAR
sbatch --job-name Minimac3_BOADICEA_X \
--cpus-per-task=1 --mem=32G --partition=long --time=720:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_Female.Non.PAR_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_Female.Non.PAR_phased_imputed" > /dev/null

# Male.Non.PAR
sbatch --job-name Minimac3_BOADICEA_X \
--cpus-per-task=1 --mem=32G --partition=long --time=720:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_Male.Non.PAR.gwas.data_imputed" > /dev/null

while squeue | grep Minimac3 > /dev/null; do sleep 5; done

rm ${bfile}*chr23* ${bfile}*chr25* ${bfile}*hh_missing*


%%bash
bfile="/home/borisov/BOADICEA/GSP_2020-041-ILL_DIA_N_1/GSP_2020-041-ILL_DIA_N_1_splitX_HNR_geno02_mind02_geno002_mind002"
n_threads=12
for chr in {1..22}; do
bgzip -f ${bfile}_ref_chr${chr}.vcf
tabix ${bfile}_ref_chr${chr}.vcf.gz
sbatch --job-name eagle_BOADICEA_chr${chr} \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chr${chr}.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom ${chr} \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chr${chr}_phased" > /dev/null
done
while squeue | grep eagle_BO > /dev/null; do sleep 5; done

%%bash
# imputation with minimac3
bfile="/home/borisov/BOADICEA/GSP_2020-041-ILL_DIA_N_1/GSP_2020-041-ILL_DIA_N_1_splitX_HNR_geno02_mind02_geno002_mind002"
n_cpus=1 # --cpus ${n_cpus}
for chr in {1..22}; do
sbatch --job-name Minimac3_BOADICEA_chr${chr} \
--cpus-per-task=${n_cpus} --mem=32G --partition=long --time=720:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chr${chr}_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chr${chr}_phased_imputed" > /dev/null
done
while squeue | grep Minimac3 > /dev/null; do sleep 5; done

%%bash
bfile="/home/borisov/BOADICEA/GSP_2020-041-ILL_DIA_N_1/GSP_2020-041-ILL_DIA_N_1_splitX_HNR_geno02_mind02_geno002_mind002"
rm ${bfile}_ref_phased_imputed.dose_subset.vcf 2> /dev/null
for chr in {1..22}; do

    echo ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf_temp/0002.vcf.gz >> \
    ${bfile}_ref_phased_imputed.dose_subset.vcf

    cat <(zgrep \# ${bfile}_ref_chr${chr}_phased_imputed.dose.vcf.gz | \
    sed 's/#CHROM/##FILTER=<ID=GENOTYPED,Description="Site was genotyped">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">\n#CHROM/') \
    <(paste -d "\t" <(zgrep -v \# ${bfile}_ref_chr${chr}_phased_imputed.dose.vcf.gz | awk '{print "chr",$1}' | sed 's/\ //g') \
    <(zgrep -v \# ${bfile}_ref_chr${chr}_phased_imputed.dose.vcf.gz | cut -f 2-)) > \
    ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf
    
    sbatch --job-name BOADICEA_postprocessing_chr${chr} --time=05:00:00 --wrap="
        bgzip -f ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf
        tabix ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf.gz
        bcftools isec -p ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf_temp \
        ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf.gz \
        /home/borisov/BOADICEA/project_flow/sample_BCAC_313_dosages.txt.gz -Oz
        "
done

while squeue | grep BOADICEA > /dev/null; do sleep 5; done


#concat and R2 filtering
bcftools concat -f ${bfile}_ref_phased_imputed.dose_subset.vcf -Ou |
bcftools view -Ov -i 'R2>0.4' ${bfile}_ref_chr${chr}_phased_imputed.dose_fixedheaderchr.vcf.gz > \
${bfile}_ref_phased_imputed.dose_BOADICEA.vcf


#concat
bcftools concat -f ${bfile}_ref_phased_imputed.dose_subset.vcf -Ov > \
${bfile}_ref_phased_imputed.dose_BOADICEA.vcf
# reformat vcf with plink
plink --vcf ${bfile}_ref_phased_imputed.dose_BOADICEA.vcf \
--keep-allele-order \
--recode vcf --out ${bfile}_ref_phased_imputed.dose_BOADICEA_format

rm -r ${bfile}_ref*_phased_imputed.dose_fixedheaderchr.vcf_temp

