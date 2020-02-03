# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
# user can set maf, geno, mind by passing 2nd 3rd and 4th argument to maf_geno_mind_hwe
# minor allele frequency (maf): excluding variants with maf < 0.01
# SNP missingness (geno): excluding variants with missingness in > 2% of samples
# Hardy–Weinberg equilibrium (hwe): excluding variants deviating from HWE with P<1e-10 (for cases and controls) and P<1-e6 (for controls only)
# individual missingness (mind): excluding samples with > 2% of missing variants 
bfile=""
maf_geno_mind_hwe () {
  bfile=$1; maf_set=$2; geno_set=$3; mind_set=$4
  maf_set="${maf_set:-0.01}"
  geno_set="${geno_set:-0.02}"
  mind_set="${mind_set:-0.02}"
  echo $bfile $maf_set $geno_set $mind_set
  plink --bfile ${bfile} \
  --maf ${maf_set} \
  --geno 0.2 \
  --hwe 1e-6 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02
  plink --bfile ${bfile}_geno02 \
  --hwe 1e-10 include-nonctrl \
  --mind 0.2 \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02
  plink --bfile ${bfile}_geno02_mind02 \
  --geno ${geno_set} \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002
  plink --bfile ${bfile}_geno02_mind02_geno002 \
  --mind ${mind_set} \
  --allow-no-sex \
  --make-bed --out ${bfile}_geno02_mind02_geno002_mind002
}
maf_geno_mind_hwe ${bfile}
grep "pass filters and QC" ${bfile}_geno02_mind02_geno002_mind002.log
echo "REMOVED SAMPLES"
cat ${bfile}*irem
