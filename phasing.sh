# shapeit in "-check" mode to double check problematic variants (duplicates, strand missmatches)
shapeit_check () {
  shapeit -check \
  --input-vcf $dir/${bfile}_ref_chr${i}.vcf \
  -M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
  --output-log $dir/${bfile}_ref_chr${i}.alignment \
  -R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
  $dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
  $dir1000genomes/1000GP_Phase3.sample \
  --include-grp $dir/group.list" >> $dir/jobs.list
}


sleep 1

for i in {1..22}; do
done
wait_fun

# shapeit phasing (-O or --output-max) --ntasks 1 --ntasks-per-node 1
for i in {1..22}; do
if [ -f $dir/${bfile}_ref_chr${i}.alignment.snp.strand.exclude ]
then
sbatch --parsable --chdir $dir --partition=medium --time=23:59:59 --cpus-per-task 16 --mem=16G --wrap="
shapeit \
--thread 24 \
--input-vcf $dir/${bfile}_ref_chr${i}.vcf \
-M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
-O $dir/${bfile}_ref_chr${i}.phased \
-R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
$dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
$dir1000genomes/1000GP_Phase3.sample \
--include-grp $dir/group.list \
--exclude-snp $dir/${bfile}_ref_chr${i}.alignment.snp.strand.exclude
" >> $dir/jobs.list
else
sbatch --parsable --chdir $dir --partition=medium --time=23:59:59 --cpus-per-task 16 --mem=16G --wrap="
shapeit \
--thread 24 \
--input-vcf $dir/${bfile}_ref_chr${i}.vcf \
-M $dir1000genomes/genetic_map_chr${i}_combined_b37.txt \
-O $dir/${bfile}_ref_chr${i}.phased \
-R $dir1000genomes/1000GP_Phase3_chr${i}.hap.gz \
$dir1000genomes/1000GP_Phase3_chr${i}.legend.gz \
$dir1000genomes/1000GP_Phase3.sample \
--include-grp $dir/group.list
" >> $dir/jobs.list
fi
sleep 0.2
done
wait_fun

rm $dir/jobs.list
