%%bash
# imputation with impute2
## Assign path to bfile and 1kg superpopulation
bfile=""
my_superpopulation=EUR # EUR, AMR, SAS, EAS, AFR
chunk_size=5000000
## Reassign default variables if run in other account / on other machine
dir1000genomes=/home/borisov/software/1000GP_Phase3/
n_threads=1
echo ${my_superpopulation} > ${bfile}_group.list
# the following run automatically
echo '#!/bin/bash
    bfile=$1
    dir1000genomes=$2
    chr=$3
    lower=$4
    upper=$5
    my_superpopulation=$6
    echo "imputing chromosome ${chr}"

    impute2 \
    -m ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt \
    -h ${dir1000genomes}/1000GP_Phase3_chr${chr}.hap.gz \
    -l ${dir1000genomes}/1000GP_Phase3_chr${chr}.legend.gz \
    -use_prephased_g \
    -known_haps_g ${bfile}_ref_chr${chr}.phased.haps \
    -int ${lower} ${upper} \
    -Ne 20000 \
    -o ${bfile}_ref_chr${chr}.phased.impute2.${upper} \
    -filt_rules_l "${my_superpopulation}<0.01" "TYPE==LOWCOV"
' > ${bfile}_imputation.sh


for chr in {1..22}; do
    chromosome_min_coord=`tail -n +2 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | head -1 | cut -d " " -f 1`
    chromosome_max_coord=`tail -n -1 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | cut -d " " -f 1`
    for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do
        # calculate upper chunk coordinate
        upper=`echo ${lower} + ${chunk_size} | bc`
        sbatch --parsable --partition=medium --time=24:00:00 --cpus-per-task ${n_threads} --mem=16G \
        ${bfile}_imputation.sh \
        ${bfile} \
        ${dir1000genomes} \
        ${chr} \
        ${lower} \
        ${upper} \
        ${my_superpopulation} >> ${bfile}_jobs.list
        sleep 1
    done
done

squeue | grep -wFf ${bfile}_jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done


echo "imputation is done"
