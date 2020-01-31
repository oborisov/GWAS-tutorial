%%bash
job_name='snptest_'
sample="_geno02_mind02_geno002_mind002_norelated_pca_ref.phased.sample"
gen_prefix='geno02_mind02_geno002_mind002_norelated_pca_ref_chr${chr}.phased.impute2'
for chr in {1..22}; do
    eval gen=${gen_prefix}
    salloc --job-name ${job_name} srun --mem=1000M --partition=medium --time=24:00:00 --cpus-per-task=1 > ${gen}_imputedPC.out_slurm 2>&1 snptest \
    -data ${gen}.gz ${sample} \
    -pheno case_control \
    -frequentist 1 \
    -method expected \
    -cov_all_continuous \
    -hwe \
    -missing_code NA \
    -assume_chromosome $chr \
    -o ${gen}_imputedPC.out &
done
#sstatus=$(sacct --format="JobName%30, State" | grep snptest_LKG_2010 | awk '{print $2}' | sort | uniq -c | awk '{print $2}')
#if [ ${sstatus} == "RUNNING" ]; then echo eq; f
#while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done
