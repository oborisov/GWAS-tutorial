
# Split data into single chromosomes and converting to vcf
split_chr () {
  bfile=$1
  plink --bfile ${study_dir}/${bfile} --chr ${i} --recode vcf --out ${study_dir}/${bfile}_chr${i}
}
for i in {1..26}; do
split_chr
done


# Processing autosomes
# coverting alleles to plus strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)

bcftools_fixref_sort_rm_dup () {
  bfile=$1
  path_to_BCFTOOLS_PLUGINS=$2
  reference_fasta=$3
  common_variants=$4
  export BCFTOOLS_PLUGINS=$2
  for chr in {1..22}; do
    bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | bcftools +fixref -- -d -f ${reference_fasta} -i ${common_variants} | bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf &
  done
}

bcftools_fixref_sort_rm_dup \
  /home/borisov/ACE/ACE_HNR_full_QC_noGSA_rs78130450 \
  /home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins/ \
  /home/borisov/software/human_g1k_v37.fasta \
  /home/borisov/software/common_all_20180423.vcf.gz

