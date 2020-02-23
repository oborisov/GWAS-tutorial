

  # qqman R script
  echo '
  options(bitmapType="cairo")
  library(qqman)
  library(QCEWAS)
  args=commandArgs(T)
  assoc_file=args[1]
  assoc=fread(assoc_file)
  assoc[CHR == "X", CHR := 23]
  assoc[, CHR := as.numeric(CHR)]
  assoc=assoc[!is.na(assoc$P)]
  assoc=assoc[P > 0]
  my_lambda=P_lambda(assoc$P)
  fwrite(data.table(filename=args, lambda=P_lambda(assoc$P)), paste0(args, "_lambda"), sep="\t")
  jpeg(paste0(assoc_file, "_manh.jpeg"), width = 12, height = 6, units = "in", res = 600)
  print(manhattan(assoc, chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T))
  title(main = assoc_file)
  dev.off()
  jpeg(paste0(assoc_file, "_qq.jpeg"), width = 6, height = 6, units = "in", res = 600)
  print(qq(assoc$P))
  title(main = assoc_file, sub = paste0("Lambda=", round(my_lambda,4)))
  dev.off()
  ' > ${bfile}.assoc.logistic.R
  Rscript ${bfile}.assoc.logistic.R ${bfile}.assoc.logistic_pc
}

# bcftools to use dbsnp ids
# | bcftools +fixref -- -d -f ${reference_fasta} -i ${common_variants} 
#     /home/borisov/software/common_all_20180423.vcf.gz \



bcftools +fixref ${bfile}_chr${chr}.vcf -Ov -o ${bfile}_flipped_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d
    bcftools +fixref ${bfile}_flipped_chr${chr}.vcf -Ov -o ${bfile}_ref_unsorted_chr${chr}.vcf -- -d -f ${reference_fasta} -i ${common_variants}
    bcftools sort ${bfile}_ref_unsorted_chr${chr}.vcf -Ov -o ${bfile}_ref_dup_chr${chr}.vcf
    bcftools norm ${bfile}_ref_dup_chr${chr}.vcf --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
 
 
### Relatedness
#  from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
king_relatedness () {
  bfile=$1
  awk '{print $2,$2,$3,$4,$5,$6}' ${bfile}.fam > \
  ${bfile}1.fam
  cp ${bfile}.fam ${bfile}.famBK
  cp ${bfile}1.fam ${bfile}.fam
  king -b ${bfile}.bed \
  --related --degree 2 --prefix ${bfile} --cpus 22
  if [ -f ${bfile}.kin0 ]; then
    first_option_cases=$(grep -wFf <(awk '{if ($6 == 2) print $1,$2}' ${bfile}.fam | tr " " "_") \
    <(awk '{if ($14 != "UN") print $1,$2}' ${bfile}.kin0 | tr " " "_") | wc -l)
    second_option_cases=$(grep -wFf <(awk '{if ($6 == 2) print $1,$2}' ${bfile}.fam | tr " " "_") \
    <(awk '{if ($14 != "UN") print $3,$4}' ${bfile}.kin0 | tr " " "_") | wc -l)
    if [ ${second_option_cases} -ge ${first_option_cases} ]; then
      plink --bfile ${bfile} \
      --remove <(awk '{if ($14 != "UN") print $1,$2}' ${bfile}.kin0) \
      --allow-no-sex \
      --make-bed --out ${bfile}_norelated
    else
      plink --bfile ${bfile} \
      --remove <(awk '{if ($14 != "UN") print $3,$4}' ${bfile}.kin0) \
      --allow-no-sex \
      --make-bed --out ${bfile}_norelated
    fi
  else
    for suff in bed bim fam; do mv ${bfile}.${suff} ${bfile}_norelated.${suff}; done
  fi
  awk '{print 0,$2,$3,$4,$5,$6}' ${bfile}_norelated.fam > \
  ${bfile}_norelated1.fam
  mv ${bfile}_norelated1.fam ${bfile}_norelated.fam
}


### Fixref

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
for chr in {1..22}; do
plink --bfile ${bfile} --chr ${chr} --recode vcf --out ${bfile}_chr${chr}
export BCFTOOLS_PLUGINS=${bcftools_plugins}
salloc --job-name fixref_${bfile} \
srun bcftools +fixref ${bfile}_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d | \
bcftools sort | bcftools norm --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
done
echo "fixref is done"
