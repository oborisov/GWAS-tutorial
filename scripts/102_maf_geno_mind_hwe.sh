%%bash
## Standard QC for genotyping data: maf, missingness of samples and variants, hwe
## Make sure that sex is present, otherwise HWE will remove the majority of X-chromosome variants
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
bfile=""
maf_set=0.001
geno_set=0.02
mind_set=0.02
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
if [ -f ${bfile}_geno02_mind02.irem ]; then cat ${bfile}_geno02_mind02.irem; fi
if [ -f ${bfile}_geno02_mind02_geno002_mind002.irem ]; then cat ${bfile}_geno02_mind02_geno002_mind002.irem; fi



%%R
bfile=""
# intersection between checksex and samples with missingness
sexcheck=fread("/home/borisov/ARM/ARM_SRC_splitX_sexcheck.sexcheck")
sexcheck=sexcheck[SNPSEX == 0]
irem=rbindlist(lapply(list.files(path="/home/borisov/ARM/", pattern="irem", full.names=T), fread))
sexcheck_irem=merge(sexcheck, irem, by.x=c("FID", "IID"), by.y=c("V1", "V2"), all=T)
print(dim(irem))
sexcheck_irem
