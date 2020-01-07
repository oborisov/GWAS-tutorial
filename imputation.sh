%%bash
## Assign path to bfile and 1kg superpopulation
bfile=
my_superpopulation=EUR # EUR, AMR, SAS, EAS, AFR
## Reassign default variables if run in other account / on other machine
dir1000genomes=/home/borisov/software/1000GP_Phase3/
chunk_size=5000000
n_threads=16
# the following run automatically
echo ${my_superpopulation} > ${bfile}_group.list
echo '#!/bin/bash
  bfile=$1
  dir1000genomes=$2
  chr=$3
  echo "imputing chromosome ${chr}"

  rm ${bfile}_ref_chr${chr}.phased.impute2 2>/dev/null
  rm ${bfile}_ref_chr${chr}.phased.impute2_info 2>/dev/null
  chromosome_min_coord=`tail -n +2 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | head -1 | cut -d " " -f 1`
  chromosome_max_coord=`tail -n -1 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | cut -d " " -f 1`
  
  for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do
    # calculate upper chunk coordinate
    upper=`echo ${lower} + ${chunk_size} | bc`
    
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

    done
    wait
  for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do
    # calculate upper chunk coordinate
    upper=`echo ${lower} + ${chunk_size} | bc`
    # concatenate .gen imputed files to whole chromosome
    cat ${bfile}_ref_chr${chr}.phased.impute2.${upper} >> \
    ${study_dir}/${bfile}_ref_chr${chr}.phased.impute2
    # concatenate _info files to whole chromosome
    if [ ! -f ${bfile}_ref_chr${chr}.phased.impute2_info ]; then
      head -1 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info > \
      ${bfile}_ref_chr${chr}.phased.impute2_info
    fi
    tail -n +2 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info >> \
    ${study_dir}/${bfile}_ref_chr${chr}.phased.impute2_info
    # remove temporary chunks files
    rm ${bfile}_ref_chr${chr}.phased.impute2.${upper}*

  done

  # compress .gen to .gen.gz
  gzip ${bfile}_ref_chr${chr}.phased.impute2_info
  gzip ${bfile}_ref_chr${chr}.phased.impute2
' > ${bfile}_imputation.sh


for chr in {1..22}; do
sbatch --parsable --partition=medium --time=24:00:00 --cpus-per-task ${n_threads} --mem=16G \
    ${bfile}_imputation.sh \
    ${bfile} \
    ${dir1000genomes} \
    ${chr} >> ${bfile}_jobs.list
    sleep 1
done

squeue | grep -wFf ${bfile}_jobs.list > /dev/null
while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done
echo "imputation done"
