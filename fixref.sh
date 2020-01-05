## Processing autosomes
# Split data into single chromosomes and converting to vcf
# coverting alleles to plus strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)

## Assign path to bfile
bfile=
# the following run automatically
echo '#!/bin/bash
bfile=$1
path_to_BCFTOOLS_PLUGINS=$2
reference_fasta=$3
common_variants=$4
chr=$5
export BCFTOOLS_PLUGINS=$2
echo "processing chromosome ${chr}"
plink --bfile ${bfile} --chr ${chr} --recode vcf --out ${bfile}_chr${chr}
bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
' > ~/fixref.sh
for chr in {1..22}; do
sbatch --parsable --mem=16G --time=05:00:00 --cpus-per-task=8 ~/fixref.sh \
    ${bfile} \
    /home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins/ \
    /home/borisov/software/human_g1k_v37.fasta \
    /home/borisov/software/common_all_20180423.vcf.gz \
    ${chr} >> ~/jobs.list
done
squeue | grep -wFf ~/jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ~/jobs.list > /dev/null; done
echo "fixref done"
