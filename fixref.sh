
# Dividing data into single chromosomes and converting to vcf

for i in {1..26}; do
  plink --bfile ${study_dir}/${bfile} --chr ${i} --recode vcf --out ${study_dir}/${bfile}_chr${i}
done

# Processing autosomes
# coverting alleles to plus strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)
export BCFTOOLS_PLUGINS=/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins/

bcftools_fixref_sort_rm_dup () {
  bfile=$1
  export_BCFTOOLS_PLUGINS="$2"
  reference_fasta=$3
  common_variants=$4
  bcftools +fixref ${bfile}_chr${i}.vcf -Ov -o ${bfile}_flipped_chr${i}.vcf -- -f ${reference_fasta} -m flip -d
  bcftools +fixref ${bfile}_flipped_chr${i}.vcf -Ov -o ${bfile}_ref_unsorted_chr${i}.vcf -- -d -f ${reference_fasta} -i ${common_variants}
  bcftools sort ${bfile}_ref_unsorted_chr${i}.vcf -Ov -o ${bfile}_ref_dup_chr${i}.vcf
  bcftools norm ${bfile}_ref_dup_chr${i}.vcf --rm-dup all -o ${bfile}_ref_chr${i}.vcf
}


for i in {1..22}; do
  bcftools_fixref_sort_rm_dup ${bfile} /home/borisov/software/human_g1k_v37.fasta /home/borisov/software/common_all_20180423.vcf.gz
done
