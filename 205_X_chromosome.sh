%%bash
bfile=""
my_superpopulation=EUR # EUR, AMR, SAS, EAS, AFR
## Reassign default variables if run in other account / on other machine
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins/"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
dir1000genomes=/home/borisov/software/1000GP_Phase3/
n_threads=16
chunk_size=5000000
echo ${my_superpopulation} > ${bfile}_group.list

# the following run automatically
# Convert data to reference strand to impute sequencing panel; # convert to vcf and change chr23 to chrX; covert alleles to reference strand using +fixref bcftools; sort files by coordinates; remov duplicated variants (bcftools)
plink --bfile ${bfile} --chr 23 --recode vcf --out ${bfile}_chr23
bcftools annotate --rename-chrs <(echo "23 X") ${bfile}_chr23.vcf | bcftools +fixref -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chrX.vcf
plink --vcf ${bfile}_ref_chrX.vcf --keep-allele-order --const-fid --update-sex <(awk '{print 0,$1,"_",$2,$5}' ${bfile}.fam | sed 's/\ _\ /_/') --pheno <(awk '{print 0,$1,"_",$2,$6}' ${bfile}.fam | sed 's/\ _\ /_/') --make-bed --out ${bfile}_ref_chrX_const_fid
awk 'OFS="\t"{print $1,$2,$2,$2}' ${bfile}_ref_chrX_const_fid.fam > ${bfile}_ref_chrX_const_fid.fam_temp
plink --bfile ${bfile}_ref_chrX_const_fid --update-ids temp --make-bed --out ${bfile}_ref_chrX
rm ${bfile}_ref_chrX_const_fid.fam_temp ${bfile}_ref_chrX_const_fid*

### Phasing
# 1. shapeit in "-check" mode to double check problematic variants (duplicates, strand missmatches) # 2. shapeit in phasing mode
echo '#!/bin/bash
bfile=$1; dir1000genomes=$2; chr=$3; n_threads=$4; shapeit -check --chrX -B ${bfile}_ref_chrX -M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt --output-log ${bfile}_ref_chr${chr}.alignment -R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz ${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz ${dir1000genomes}/1000GP_Phase3.sample --include-grp ${bfile}_group.list; if [ -f ${bfile}_ref_chr${chr}.alignment.snp.strand.exclude ]; then shapeit --chrX --thread ${n_threads} -B ${bfile}_ref_chrX -M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt -O ${bfile}_ref_chr${chr}.phased -R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz ${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz ${dir1000genomes}/1000GP_Phase3.sample --include-grp ${bfile}_group.list --exclude-snp ${bfile}_ref_chr${chr}.alignment.snp.strand.exclude; else shapeit --chrX --thread ${n_threads} -B ${bfile}_ref_chrX -M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt -O ${bfile}_ref_chr${chr}.phased -R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz ${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz ${dir1000genomes}/1000GP_Phase3.sample --include-grp ${bfile}_group.list; fi' > ${bfile}_phasing.sh
for chr in X_PAR1 X_PAR2 X_NONPAR; do sbatch --parsable --partition=medium --time=24:00:00 --cpus-per-task ${n_threads} --mem=16G ${bfile}_phasing.sh ${bfile} ${dir1000genomes} ${chr} ${n_threads} >> ${bfile}_jobs.list; done
squeue | grep -wFf ${bfile}_jobs.list > /dev/null; while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done

### Imputation
echo '#!/bin/bash
bfile=$1; dir1000genomes=$2; chr=$3 lower=$4; upper=$5; my_superpopulation=$6; impute2 -chrX -m ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt -h ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz -l ${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz -use_prephased_g -known_haps_g ${bfile}_ref_chr${chr}.phased.haps -sample_g ${bfile}_ref_chr${chr}.phased.sample -int ${lower} ${upper} -Ne 20000 -o ${bfile}_ref_chr${chr}.phased.impute2.${upper} -filt_rules_l "${my_superpopulation}<0.01" "TYPE==LOWCOV"' > ${bfile}_imputation.sh
for chr in X_PAR1 X_PAR2 X_NONPAR; do chromosome_min_coord=`tail -n +2 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | head -1 | cut -d " " -f 1`; chromosome_max_coord=`tail -n -1 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | cut -d " " -f 1`; for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do upper=`echo ${lower} + ${chunk_size} | bc`; sbatch --parsable --partition=medium --time=24:00:00 --cpus-per-task 1 --mem=16G ${bfile}_imputation.sh ${bfile} ${dir1000genomes} ${chr} ${lower} ${upper} ${my_superpopulation} >> ${bfile}_jobs.list; sleep 1; done; done
squeue | grep -wFf ${bfile}_jobs.list > /dev/null; while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done; echo "imputation is done"

### Concat
rm ${bfile}_ref_chrX.phased.impute2 2>/dev/null; rm ${bfile}_ref_chrX.phased.impute2_info 2>/dev/null
for chr in X_PAR1 X_NONPAR X_PAR2; do 
    chromosome_min_coord=`tail -n +2 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | head -1 | cut -d " " -f 1`; 
    chromosome_max_coord=`tail -n -1 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | cut -d " " -f 1`; 
    for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do 
    upper=`echo ${lower} + ${chunk_size} | bc`;
    cat ${bfile}_ref_chr${chr}.phased.impute2.${upper} >> ${bfile}_ref_chrX.phased.impute2;
    if [ ! -s ${bfile}_ref_chrX.phased.impute2_info ]; then
    head -1 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info > ${bfile}_ref_chrX.phased.impute2_info;
    fi;
    tail -n +2 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info >> ${bfile}_ref_chrX.phased.impute2_info;
    rm ${bfile}_ref_chr${chr}.phased.impute2.${upper} ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info ${bfile}_ref_chr${chr}.phased.impute2.${upper}_samples;
    mv ${bfile}_ref_chr${chr}.phased.impute2.${upper}_summary ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info_by_sample ${bfile}_ref_chr${chr}.phased.impute2.${upper}_warnings ${bfile}_imputation_quality
    done
done
gzip -f ${bfile}_ref_chrX.phased.impute2_info; gzip -f ${bfile}_ref_chrX.phased.impute2;

echo "X chromosome is done"
