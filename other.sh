

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

# | bcftools +fixref -- -d -f ${reference_fasta} -i ${common_variants}


bcftools +fixref ${bfile}_chr${chr}.vcf -Ov -o ${bfile}_flipped_chr${chr}.vcf -- -f ${reference_fasta} -m flip -d
    bcftools +fixref ${bfile}_flipped_chr${chr}.vcf -Ov -o ${bfile}_ref_unsorted_chr${chr}.vcf -- -d -f ${reference_fasta} -i ${common_variants}
    bcftools sort ${bfile}_ref_unsorted_chr${chr}.vcf -Ov -o ${bfile}_ref_dup_chr${chr}.vcf
    bcftools norm ${bfile}_ref_dup_chr${chr}.vcf --rm-dup all -o ${bfile}_ref_chr${chr}.vcf
 
