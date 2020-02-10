%%bash
### Computing lambda
bfile=""
# plink logistic with 4 pc
bfile="${bfile}_ref.phased.sample_concat"
plink --bfile ${bfile} \
--logistic \
--covar ${bfile}_eigen.eigenvec \
--covar-name PC1, PC2, PC3, PC4 \
--out ${bfile}_lambda
# select p values for phenotype only
awk '{if ($5 == "ADD") print $9}' ${bfile}_lambda.assoc.logistic > ${bfile}_lambda_P
# remove temp files
rm ${bfile}_lambda.assoc.logistic
# Computing lambda
Rscript -e 'bfile=commandArgs(T); dt=data.table(lambda_4pc=median(qchisq(1 - fread(paste0(bfile, "_lambda_P"), header=F)[[1]], 1), na.rm=T) / qchisq(0.5, 1)); system(paste0("rm ", paste0(bfile, "_lambda_P"))); fwrite(dt, paste0(bfile, "_lambda")); dt' ${bfile} > ${bfile}_lambda
#####################################################

%%R
#creating sample with PCs
bfile=""
# the rest will run automatically
shapeit_sample=fread(paste0(bfile, "_ref_chr1.phased.sample"), header=T)
eigenvec=fread(paste0(bfile, "_ref.phased.sample_concat_eigen.eigenvec"), header=T)
fam=fread(paste0(bfile, "_ref.phased.sample_concat.fam"), header=F)
colnames(eigenvec)[2]="ID_1"
eigenvec=eigenvec[,-1]
colnames(fam)[c(2,5:6)]=c("ID_1", "sex", "case_control")
fam=fam[,c(2,5:6)]
shapeit_sample_pc_fam=Reduce(merge, list(shapeit_sample, eigenvec, fam))
shapeit_sample_pc_fam=shapeit_sample_pc_fam[match(shapeit_sample$ID_1, shapeit_sample_pc_fam$ID_1)]
if (!all(shapeit_sample_pc_fam$ID_1 == shapeit_sample$ID_1, na.rm=T)) {
    stop("sample with PC does not match shapeit sample!")
}
shapeit_sample_pc_fam=shapeit_sample_pc_fam[,c("ID_1","ID_2", "missing", "sex", "case_control", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
shapeit_sample_pc_fam[, case_control:=case_control-1]
shapeit_sample_pc_fam=shapeit_sample_pc_fam[-1]
fwrite(shapeit_sample_pc_fam, paste0(bfile, "_ref.phased.sample"), sep=" ", na="NA", quote=F)
#####################################################

%%bash
bfile=""
# the rest will run automatically
job_suffix=$(echo $bfile | sed 's/.*\///')
job_name="snptest_${job_suffix}"
cat <(head -1 ${bfile}_ref.phased.sample) \
<(echo "0 0 0 D B C C C C C C C C C C") \
<(tail -n +2 ${bfile}_ref.phased.sample) > \
${bfile}_ref.phased.sample2
gen_prefix=${bfile}'_ref_chr${chr}.phased.impute2'
for chr in {1..22}; do
    eval gen=${gen_prefix}
    salloc --job-name ${job_name} --mem=1000M --partition=long --time=72:00:00 --cpus-per-task=1 > ${gen}.out_slurm 2>&1 srun snptest \
    -data ${gen}.gz ${bfile}_ref.phased.sample2 \
    -pheno case_control \
    -frequentist 1 \
    -method expected \
    -cov_names PC1 PC2 PC3 PC4 \
    -hwe \
    -missing_code NA \
    -assume_chromosome $chr \
    -o ${gen}_imputedPC.out.gz &
done


#sstatus=$(sacct --format="JobName%30, State" | grep snptest_LKG_2010 | awk '{print $2}' | sort | uniq -c | awk '{print $2}')
#if [ ${sstatus} == "RUNNING" ]; then echo eq; f
#while [ $? -ne 1 ]; do sleep 5; squeue | grep -wFf ${bfile}_jobs.list > /dev/null; done
