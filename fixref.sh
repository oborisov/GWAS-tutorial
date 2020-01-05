
# Split data into single chromosomes and converting to vcf
#  | bcftools +fixref -- -d -f ${reference_fasta} -i ${common_variants}
split_chr () {
  bfile=$1
  plink --bfile ${study_dir}/${bfile} --chr ${i} --recode vcf --out ${study_dir}/${bfile}_chr${i}
}
for i in {1..26}; do
split_chr
done


## Processing autosomes
# coverting alleles to plus strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)
bfile=
working_dir=
echo '#!/bin/bash
bfile=$1
path_to_BCFTOOLS_PLUGINS=$2
reference_fasta=$3
common_variants=$4
chr=$5
export BCFTOOLS_PLUGINS=$2
echo "processing chromosome ${chr}"
bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
' > script.sh
for chr in {1..22}; do
sbatch --time=05:00:00 --cpus-per-task=8 script.sh \
    ${bfile} \
    /home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins/ \
    /home/borisov/software/human_g1k_v37.fasta \
    /home/borisov/software/common_all_20180423.vcf.gz \
    ${chr} >> ${working_dir}/jobs.list
done

squeue | grep -wFf ${working_dir}/jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${working_dir}/jobs.list > /dev/null; done
echo "fixref done"
