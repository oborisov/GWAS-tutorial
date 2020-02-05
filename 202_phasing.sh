%%bash
# Phasing with shapeit in 2 steps:
# 1. shapeit in "-check" mode to double check problematic variants (duplicates, strand missmatches)
# 2. shapeit in phasing mode
## Assign path to bfile and 1kg superpopulation
bfile=""
my_superpopulation=EUR # EUR, AMR, SAS, EAS, AFR
## Reassign default variables if run in other account / on other machine
dir1000genomes=/home/borisov/software/1000GP_Phase3/
n_threads=16
echo ${my_superpopulation} > ${bfile}_group.list
# the following run automatically
echo '#!/bin/bash
bfile=$1
dir1000genomes=$2
chr=$3
n_threads=$4

echo "phasing chromosome ${chr} (shapeit check)"

shapeit \
-check \
--input-vcf ${bfile}_ref_chr${chr}.vcf \
-M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt \
--output-log ${bfile}_ref_chr${chr}.alignment \
-R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz \
${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz \
${dir1000genomes}/1000GP_Phase3.sample \
--include-grp ${bfile}_group.list

echo "phasing chromosome ${chr} (shapeit phase)"

if [ -f ${bfile}_ref_chr${chr}.alignment.snp.strand.exclude ]; then
shapeit \
--thread ${n_threads} \
--input-vcf ${bfile}_ref_chr${chr}.vcf \
-M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt \
-O ${bfile}_ref_chr${chr}.phased \
-R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz \
${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz \
${dir1000genomes}/1000GP_Phase3.sample \
--include-grp ${bfile}_group.list \
--exclude-snp ${bfile}_ref_chr${chr}.alignment.snp.strand.exclude

else

shapeit \
--thread ${n_threads} \
--input-vcf ${bfile}_ref_chr${chr}.vcf \
-M ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt \
-O ${bfile}_ref_chr${chr}.phased \
-R ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz \
${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz \
${dir1000genomes}/1000GP_Phase3.sample \
--include-grp ${bfile}_group.list
fi

' > ${bfile}_phasing.sh

for chr in {1..22}; do
sbatch --parsable --partition=medium --time=24:00:00 --cpus-per-task ${n_threads} --mem=16G \
    ${bfile}_phasing.sh \
    ${bfile} \
    ${dir1000genomes} \
    ${chr} \
    ${n_threads} >> ${bfile}_jobs.list
    sleep 1
done

squeue | grep -wFf ${bfile}_jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done
echo "phasing is done"

