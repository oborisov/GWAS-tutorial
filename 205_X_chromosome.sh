### X ###

# rename 23 to X
for i in 23; do
sbatch --parsable --chdir $dir --mem=16G --time=4:00:00 --wrap="
bcftools annotate -Ov --rename-chrs /home/borisov/software/ucsc2nsembl.txt \
$dir/${bfile}_chr${i}.vcf > \
$dir/${bfile}_renamed_chrX.vcf;
" >> $dir/jobs.list
done
wait_fun
#mv $dir/${bfile}_renamed_chr${i}.vcf $dir/${bfile}_renamed_chrX.vcf

# force the unambigous alleles onto the forward strand
for i in 23; do
sbatch --parsable --chdir $dir --mem=16G --time=4:00:00 --wrap="
${export_BCFTOOLS_PLUGINS}
bcftools +fixref $dir/${bfile}_renamed_chrX.vcf \
-Ov -o $dir/${bfile}_renamed_flipped_chrX.vcf \
-- -f /home/borisov/software/human_g1k_v37.fasta -m flip -d
" >> $dir/jobs.list
done
wait_fun

# bgzip
sbatch --parsable --chdir $dir --mem=16G --time=4:00:00 --wrap="
bgzip -fc $dir/${bfile}_renamed_flipped_chrX.vcf > \
$dir/${bfile}_renamed_flipped_chrX.vcf.gz;
bcftools index $dir/${bfile}_renamed_flipped_chrX.vcf.gz
" >> $dir/jobs.list
wait_fun

# check sanity
echo 'library(data.table); library(ggplot2)
data.dist.txt_name=commandArgs(T)[1]
data.dist.txt=fread(data.dist.txt_name, header=F)
ggplot(data.dist.txt, aes(x=V2, y=V4))+
  aes(y = stat(data.dist.txt$V4 / max(data.dist.txt$V4))) +
  geom_density(stat = "identity") +
  scale_x_continuous(breaks=seq(0,1,0.1)) +
  scale_y_continuous(breaks=seq(0,1,0.1))
ggsave(paste0(data.dist.txt_name, ".pdf"), plot = last_plot())
' > $dir/sanity.R
sbatch --parsable --chdir $dir --mem=16G --partition=medium --time=23:00:00 --wrap="
${export_BCFTOOLS_PLUGINS}
bcftools annotate -c INFO/AF \
-a /home/borisov/software/1000GP_Phase3/vcf/af.vcf.gz \
$dir/${bfile}_renamed_flipped_chrX.vcf.gz | \
bcftools +af-dist | grep ^PROB > $dir/${bfile}_renamed_flipped_chrX.vcf_data.dist.txt;
Rscript sanity.R $dir/${bfile}_renamed_flipped_chrX.vcf_data.dist.txt
"

#fixref
sbatch --parsable --chdir $dir --mem=16G --time=4:00:00 --wrap="
${export_BCFTOOLS_PLUGINS}
bcftools +fixref $dir/${bfile}_renamed_flipped_chrX.vcf \
-Ov -o $dir/${bfile}_ref_chrX.vcf \
-- -d -f /home/borisov/software/human_g1k_v37.fasta \
-i /home/borisov/software/common_all_20180423.vcf.gz" >> $dir/jobs.list
wait_fun

# convert back to plink; check if fam is sample; write sex
sbatch --parsable --mem=16G --time=4:00:00 --chdir $dir --wrap="
plink --vcf $dir/${bfile}_ref_chrX.vcf --keep-allele-order \
--const-fid \
--make-bed --out $dir/${bfile}_ref_chrX;
" >> $dir/jobs.list
wait_fun
awk '{print $1,$2}' $dir/${bfile}_ref_chrX.fam | tr " " "." > $dir/new_fam_col1
awk '{print $1,$2}' $dir/${bfile}.fam | tr " " "." > $dir/old_fam_col1
#n_missm_lines=$(comm -3 $dir/old_fam_col1 $dir/new_fam_col1 | wc -l | bc)
#if [ $n_missm_lines -gt 0 ]; then
diff <(awk '{print $1,$2}' $dir/${bfile}.fam | tr " " "_") <(awk '{print $2}' $dir/${bfile}_ref_chrX.fam)
if [ ! $? ]; then
echo "Problem with X imputation, fam files do not match" > $dir/ERROR_MESSAGE
exit
fi
cp $dir/${bfile}.fam $dir/${bfile}_ref_chrX.fam

# shapeit -check
for i in X_PAR1 X_PAR2 X_NONPAR; do
sbatch --parsable --mem=16G --time=4:00:00 --chdir $dir --wrap="shapeit -check --chrX \
-B $dir/${bfile}_ref_chrX \
-M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
--output-log $dir/${bfile}_ref_chr${i}.alignment \
-R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
$dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
$dir1000genomes/1000GP_Phase3.sample \
--include-grp $dir/group.list" >> $dir/jobs.list
sleep 10
done
wait_fun

# shapeit phasing (-O or --output-max)  --time=4:00:00
for i in X_PAR1 X_PAR2 X_NONPAR; do
if [ -f $dir/${bfile}_ref_chr${i}.alignment.snp.strand.exclude ]
then
sbatch --parsable --partition=medium --time=23:59:59 --cpus-per-task 24 --mem=32G --chdir $dir --wrap="shapeit --chrX \
--thread 8 \
-B $dir/${bfile}_ref_chrX \
-M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
-O $dir/${bfile}_ref_chr${i}.phased \
-R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
$dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
$dir1000genomes/1000GP_Phase3.sample \
--include-grp $dir/group.list \
--exclude-snp $dir/${bfile}_ref_chr${i}.alignment.snp.strand.exclude" >> $dir/jobs.list
else
sbatch --parsable --partition=medium --time=23:59:59 --cpus-per-task 24 --mem=32G --chdir $dir --wrap="shapeit --chrX \
--thread 8 \
-B $dir/${bfile}_ref_chrX \
-M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
-O $dir/${bfile}_ref_chr${i}.phased \
-R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
$dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
$dir1000genomes/1000GP_Phase3.sample \
--include-grp $dir/group.list" >> $dir/jobs.list
fi
sleep 20
done
wait_fun


# impute2
#without "-filt_rules_l" => --partition=long + --time=80:00:00
for i in X_PAR1 X_PAR2 X_NONPAR; do
for j in $jseq; do
upper=`echo ${chunk_size}*$j | bc`; lower=`echo $upper - ${chunk_size} | bc`
sbatch  --parsable --mem=8G --time=4:50:00 --chdir $dir --wrap="impute2 -chrX \
-m $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
-h $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
-l $dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
-use_prephased_g \
-known_haps_g $dir/${bfile}_ref_chr${i}.phased.haps \
-sample_g $dir/${bfile}_ref_chr${i}.phased.sample \
-int $lower $upper -Ne 20000 \
-o $dir/${bfile}_ref_chr${i}.phased.impute2.${upper} \
$my_filt_rules 'TYPE==LOWCOV'" >> $dir/jobs.list
sleep 2
done
done
wait_fun

# concatenating
#head -1 $dir/${bfile}_ref_chr${i}.phased.impute2.${upper}_info > \
for i in X_PAR1 X_PAR2 X_NONPAR; do
nl=$(ls $dir/${bfile}_ref_chr${i}.phased.impute2.*0000*_info 2> /dev/null | wc -l)
if [ $nl -gt 0 ]; then
  rm $dir/${bfile}_ref_chr${i}.phased.impute2
  head -1 `ls $dir/${bfile}_ref_chr${i}.phased.impute2.*0000*_info | head -1` > \
  $dir/${bfile}_ref_chr${i}.phased.impute2_info
  for j in $jseq; do
  upper=$(echo ${chunk_size}*$j | bc); lower=$(echo $upper - ${chunk_size} | bc)
  cat $dir/${bfile}_ref_chr${i}.phased.impute2.${upper} >> \
  $dir/${bfile}_ref_chr${i}.phased.impute2
  tail -n +2 $dir/${bfile}_ref_chr${i}.phased.impute2.${upper}_info >> \
  $dir/${bfile}_ref_chr${i}.phased.impute2_info
  done
fi
done

# concatenating X_PAR1 X_PAR2 X_NONPAR
rm $dir/${bfile}_ref_chrX.phased.impute2
head -1 `ls $dir/${bfile}_ref_chr${i}.phased.impute2_info | head -1` > \
$dir/${bfile}_ref_chrX.phased.impute2_info
for i in X_PAR1 X_PAR2 X_NONPAR; do
cat $dir/${bfile}_ref_chr${i}.phased.impute2 >> \
$dir/${bfile}_ref_chrX.phased.impute2
tail -n +2 $dir/${bfile}_ref_chr${i}.phased.impute2_info >> \
$dir/${bfile}_ref_chrX.phased.impute2_info
done

# create sample
sed 's/plink_pheno/phenotype/' $dir/${bfile}_ref_chrX_NONPAR.phased.sample > \
$dir/${bfile}_ref_chrX.phased.sample_pheno

# removing X_PAR1 X_PAR2 X_NONPAR
for i in X_PAR1 X_PAR2 X_NONPAR; do
rm $dir/${bfile}_ref_chr${i}.phased.impute2*
done

for i in X; do
awk '{ if ($7 > 0.3) print $2}' $dir/${bfile}_ref_chr${i}.phased.impute2_info > \
$dir/${bfile}_ref_chr${i}.phased.impute2_info0.3_snp
done

rm $dir/${bfile}_ref_chrX.phased.impute2_info_snp
for i in X; do
awk '{ print $2,$7}' $dir/${bfile}_ref_chr${i}.phased.impute2_info >> \
$dir/${bfile}_ref_chrX.phased.impute2_info_snp
done


### Post_filtering
# plink --extract
for i in X; do
sbatch --parsable --mem=16G --partition=medium --time=23:59:59 --chdir $dir --wrap="
plink --gen $dir/${bfile}_ref_chr${i}.phased.impute2 \
--sample $dir/${bfile}_ref_chr${i}.phased.sample_pheno \
--extract $dir/${bfile}_ref_chr${i}.phased.impute2_info_snp \
--oxford-single-chr 23 \
--make-bed --out $dir/${bfile}_ref_chr${i}.phased.impute2
"  >> $dir/jobs.list
done
wait_fun

rm $dir/*col1
