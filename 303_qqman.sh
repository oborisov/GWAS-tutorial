%%R
# set up 3 variables
bfile=""
info_threshold=0.8
maf_threshold=0.01

# the rest will run automatically

# libraries
lapply(c("qqman", "QCEWAS", "RColorBrewer"), require, character.only = TRUE)

# reading snptest output into dt
snptest_out_list_files=list.files(path=gsub("(.*/).*","\\1",bfile), pattern=paste0(gsub(".*/","",bfile), ".*.out.gz$"), full.names=T)
snptest_out_dt=rbindlist(sapply(snptest_out_list_files, function(x) {
    fread(cmd=paste0("zcat ", x, " | grep -v \\#"))[,c("rsid", "chromosome", "position", "info", "cases_maf", "controls_maf", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1", "alleleB", "alleleA", "all_total")]
}, simplify=F))
colnames(snptest_out_dt)=c("SNP", "CHR", "BP", "info", "cases_maf", "controls_maf", "P", "BETA", "SE", "effectAllele", "otherAllele", "NMISS")

# filtering snptest using info and maf
snptest_out_dt=snptest_out_dt[info > info_threshold & cases_maf > maf_threshold & controls_maf > maf_threshold]

# adjusting CHR and p
snptest_out_dt[CHR == "X", CHR := "23"]
snptest_out_dt[CHR == "0X", CHR := "23"]
snptest_out_dt[, CHR := as.numeric(CHR)]
snptest_out_dt=snptest_out_dt[!is.na(snptest_out_dt$P) & P > 0]

# calculating lambda https://www.biostars.org/p/298847/
lambda <- median(qchisq(1 - snptest_out_dt$P, 1)) / qchisq(0.5, 1)

# producing Manhattan plot
jpeg(paste0(bfile, "_manh.jpeg"),
     width = 12, height = 6, units = "in", res = 600)
print(manhattan(rbind(snptest_out_dt[P<5e-2], snptest_out_dt[P>5e-2][seq(1,nrow(snptest_out_dt[P>5e-2]),10)]), chr="CHR", bp="BP", p="P", snp="SNP",
                annotatePval = 1, annotateTop = T,
                col=brewer.pal(8, "Dark2")))
title(main = gsub(".*/(.*)_checkedsex.*","\\1",bfile))
dev.off()

# producing Q-Q plot
jpeg(paste0(bfile, "_qq.jpeg"),
     width = 6, height = 6, units = "in", res = 600)
print(qq(snptest_out_dt$P))
title(main = gsub(".*/(.*)_checkedsex.*","\\1",bfile),
      sub = paste0("Lambda=", round(lambda,4)))
dev.off()

# writing suggestive associations into a file
fwrite(snptest_out_dt, paste0(bfile, "_all.tsv"), sep="\t", na=NA, quote=F)
fwrite(snptest_out_dt[P < 1e-5], paste0(bfile, "_suggestive.tsv"), sep="\t", na=NA, quote=F)

system(paste0("gzip -f ", bfile, "_all.tsv"))
system(paste0("gzip -f ", bfile, "_suggestive.tsv"))

# ploting Manhattan and Q-Q plot to stdout
print(manhattan(snptest_out_dt, chr="CHR", bp="BP", p="P", snp="SNP",
                annotatePval = 1, annotateTop = T,
                col=brewer.pal(8, "Dark2")))
print(qq(snptest_out_dt$P))
