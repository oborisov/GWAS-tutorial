%%bash
#convert data to reference strand to impute sequencing panel
## Check that FIDs or IIDs do not have underscores "_"
# Processing autosomes
# Split data into single chromosomes and converting to vcf
# coverting alleles to reference strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)
## Assign path to bfile
bfile=""
## Reassign default variables if run in other account / on other machine
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
