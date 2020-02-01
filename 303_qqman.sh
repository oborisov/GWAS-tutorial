%%R
# set up 4 variables
path_to_wd="/home/borisov/"
pattern_of_snptest="_geno02_mind02_geno002_mind002_norelated_pca_ref_chr.*.phased.impute2_imputedPC.out$"
info_threshold=0.8
maf_threshold=0.01

## the rest should run automatically

# xcat default is bitmapType
# getOption("bitmapType")
# [1] "Xlib"
# but it requires X11 which is not present
# to overcome, use cairo graphics (https://github.com/easybuilders/easybuild-easyconfigs/issues/4918)
options(bitmapType="cairo")

# libraries
lapply(c("qqman", "QCEWAS", "RColorBrewer"), require, character.only = TRUE)

# reading snptest output into dt
snptest_out_list_files=list.files(path=path_to_wd, pattern=pattern_of_snptest, full.names=T)
snptest_out_dt=rbindlist(sapply(snptest_out_list_files, function(x) {
    fread(cmd=paste0("grep -v \\# ", x))[,c("rsid", "chromosome", "position", "info", "cases_maf", "controls_maf", "frequentist_add_pvalue")]
}, simplify=F))
colnames(snptest_out_dt)=c("SNP", "CHR", "BP", "info", "cases_maf", "controls_maf", "P")

# filtering snptest using info and maf
snptest_out_dt=snptest_out_dt[info > info_threshold & cases_maf > maf_threshold & controls_maf > maf_threshold]

# adjusting CHR and p
snptest_out_dt[CHR == "X", CHR := 23]
snptest_out_dt[, CHR := as.numeric(CHR)]
snptest_out_dt=snptest_out_dt[!is.na(snptest_out_dt$P) & P > 0]

# calculating lambda https://www.biostars.org/p/298847/
lambda <- median(qchisq(1 - snptest_out_dt$P, 1)) / qchisq(0.5, 1)

# producing Manhattan plot
path_to_wd="/home/borisov/nsCLP/LKG_2010/"
pattern_of_snptest="LKG_2010_cases_399_controls_1318_geno02_mind02_geno002_mind002_norelated_pca_ref_chr.*.phased.impute2_imputedPC.out$"

plot_prefix=paste0(path_to_wd, gsub("\\$", "", gsub("_chr\\.\\*", "", pattern_of_snptest)))

jpeg(paste0(plot_prefix, "_manh.jpeg"),
     width = 12, height = 6, units = "in", res = 600)
print(manhattan(snptest_out_dt, chr="CHR", bp="BP", p="P", snp="SNP",
                annotatePval = 1, annotateTop = T,
                col=brewer.pal(8, "Dark2")))
title(main = plot_prefix)
dev.off()

# producing Q-Q plot
jpeg(paste0(plot_prefix, "_qq.jpeg"),
     width = 6, height = 6, units = "in", res = 600)
print(qq(snptest_out_dt$P))
title(main = plot_prefix,
      sub = paste0("Lambda=", round(lambda,4)))
dev.off()


# ploting Manhattan and Q-Q plot to stdout
print(manhattan(snptest_out_dt, chr="CHR", bp="BP", p="P", snp="SNP",
                annotatePval = 1, annotateTop = T,
                col=brewer.pal(8, "Dark2")))
print(qq(snptest_out_dt$P))

