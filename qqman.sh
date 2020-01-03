options(bitmapType="cairo")
library(qqman)
library(RColorBrewer)
sumstats=fread("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic")
colnames(sumstats)=c("CHR", "BP", "SNP", "A1", "P")
sumstats[CHR == "X", CHR := 23]
sumstats[, CHR := as.numeric(CHR)]
sumstats=sumstats[!is.na(CHR)]
sumstats=sumstats[!is.na(P)]
sumstats=sumstats[P > 0]
jpeg("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic_manh.jpeg", width = 12, height = 6, units = "in", res = 200)
print(manhattan(sumstats, chr="CHR", bp="BP", p="P", snp="SNP", annotatePval = 1, annotateTop = T, col=brewer.pal(8, "Dark2")))
title(main = "ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic")
dev.off()
my_lambda=median(qchisq(1 - (sumstats[[5]]), 1))/qchisq(0.5,1)
jpeg("/home/borisov/ACE/ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic_qq.jpeg", width = 6, height = 6, units = "in", res = 200)
print(qq(sumstats[[5]]))
title(main = "ACE_HNR_full_QC_PC3-5.PHENO1.glm.logistic", sub = paste0("Lambda=", round(my_lambda,4)))
dev.off()
