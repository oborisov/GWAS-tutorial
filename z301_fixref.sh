%%bash
#convert data to reference strand to impute sequencing panel
## Check that FIDs or IIDs do not have underscores "_"
# Processing autosomes
# Split data into single chromosomes and converting to vcf
# coverting alleles to reference strand using +fixref bcftools
# sorting files by coordinates
# removing duplicated variants (bcftools)
## Assign path to bfile
bfile="/home/borisov/LUTO/Dutch/dcg_checkedsex_geno02_mind02_geno002_mind002_norelated_pca_pca_pca_pca"
## Reassign default variables if run in other account / on other machine
bcftools_plugins="/home/borisov/.conda/envs/borisov_env/bcftools-1.9/plugins"
reference_fasta="/home/borisov/software/human_g1k_v37.fasta"
# the following run automatically
for chr in {1..22}; do
plink --bfile ${bfile} --chr ${chr} --recode vcf --out ${bfile}_chr${chr}
export BCFTOOLS_PLUGINS=${bcftools_plugins}
salloc --job-name fixref_${bfile} \
srun bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | \
bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
done
echo "fixref is done"
