#########################################################################################
### Fixing reference allele #############################################################
#########################################################################################
%%bash
# Convert data to reference strand to impute sequencing panel
bfile=""
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
for chr in {1..22}; do
    export BCFTOOLS_PLUGINS=${bcftools_plugins}
    sbatch --job-name fixref_chr${chr} \
    --cpus-per-task=1 --mem=16G --time=05:00:00 \
    --wrap="
    plink --bfile ${bfile} --chr ${chr} --recode vcf --out ${bfile}_chr${chr}
    bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chr${chr}.vcf.gz
    bcftools index ${bfile}_ref_chr${chr}.vcf.gz
    " > /dev/null
done
while squeue | grep fixref > /dev/null; do sleep 5; done
for chr in {1..22}; do rm ${bfile}_chr${chr}.vcf ${bfile}_chr${chr}.log; done # removing input vcfs (temrorary files)

#########################################################################################
### PHASING #############################################################################
#########################################################################################
%%bash
bfile=""
# Phasing
n_threads=12
for chr in {1..22}; do
sbatch --job-name eagle_chr${chr} \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chr${chr}.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom ${chr} \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chr${chr}_phased" > /dev/null
done
while squeue | grep eagle > /dev/null; do sleep 5; done

#########################################################################################
### IMPUTATION ##########################################################################
#########################################################################################
%%bash
# imputation with minimac3
bfile=""
n_cpus=1 # --cpus ${n_cpus}
for chr in {1..22}; do
sbatch --job-name Minimac3_chr${chr} \
--cpus-per-task=${n_cpus} --mem=48G --partition=medium --time=24:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chr${chr}_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chr${chr}_phased_imputed" > /dev/null
done
while squeue | grep Minimac3 > /dev/null; do sleep 5; done


#########################################################################################

%%bash
### Chromosome X
bfile=""
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
n_threads=12

### Preprocessing
# Excluding samples with undetermined sex
# Removing males haploid heterozygous calls (should be done when X is splitted into PAR and nonPAR)
plink --bfile ${bfile} \
--remove <(awk '{if ($5 == 0) print $0}' ${bfile}.fam) \
--set-hh-missing \
--make-bed --out ${bfile}_hh_missing

# Convering data to vcf
plink --bfile ${bfile}_hh_missing --chr 25 --recode vcf --out ${bfile}_chr25
plink --bfile ${bfile}_hh_missing --chr 23 --filter-females --recode vcf --out ${bfile}_chr23_females
plink --bfile ${bfile}_hh_missing --chr 23 --filter-males --recode vcf --out ${bfile}_chr23_males

# Fixing reference allele
bcftools annotate --rename-chrs <(echo "25 X") ${bfile}_chr25.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_PAR.vcf.gz
bcftools annotate --rename-chrs <(echo "23 X") ${bfile}_chr23_females.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz
bcftools annotate --rename-chrs <(echo "23 X") ${bfile}_chr23_males.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -Oz -o ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz
bcftools index ${bfile}_ref_chrX_PAR.vcf.gz; bcftools index ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz; bcftools index ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz

# Remove temporary files
rm ${bfile}*chr23* ${bfile}*chr25* ${bfile}*hh_missing*

### Phasing
# Phasing PAR
sbatch --job-name eagle_chrX \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chrX_PAR.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr23.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom 23 \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chrX_PAR_phased" > /dev/null
# Phasing Female.Non.PAR
sbatch --job-name eagle_chrX \
--cpus-per-task=${n_threads} --mem=16G --time=05:00:00 \
--wrap="eagle \
--geneticMapFile /home/borisov/software/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz \
--vcfTarget ${bfile}_ref_chrX_Female.Non.PAR.vcf.gz \
--vcfRef /home/borisov/software/1000GP_Phase3/vcf/ALL.chr23.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_Ob.bcf \
--chrom 23 \
--numThreads ${n_threads} \
--outPrefix ${bfile}_ref_chrX_Female.Non.PAR_phased" > /dev/null
# Male.Non.PAR do not need to be phased as they are hemizygous

while squeue | grep eagle > /dev/null; do sleep 5; done

### Imputation 
# Imputing PAR
sbatch --job-name Minimac3_chrX \
--cpus-per-task=1 --mem=48G --partition=medium --time=24:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_PAR_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_PAR_phased_imputed" > /dev/null
# Imputing Female.Non.PAR
sbatch --job-name Minimac3_chrX \
--cpus-per-task=1 --mem=48G --partition=medium --time=24:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_Female.Non.PAR_phased.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_Female.Non.PAR_phased_imputed" > /dev/null
# Imputing Male.Non.PAR
sbatch --job-name Minimac3_chrX \
--cpus-per-task=1 --mem=48G --partition=medium --time=24:00:00 \
--wrap="Minimac3 \
--refHaps /home/borisov/software/G1K_P3_M3VCF_FILES_WITH_ESTIMATES/X.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${bfile}_ref_chrX_Male.Non.PAR.vcf.gz \
--noPhoneHome \
--rsid \
--lowMemory \
--prefix ${bfile}_ref_chrX_Male.Non.PAR.gwas.data_imputed" > /dev/null
while squeue | grep Minimac3 > /dev/null; do sleep 5; done

### Copy statistics
