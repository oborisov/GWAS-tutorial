pedfile_prefix=""
sed -n 100000p ${pedfile_prefix}.ped

... are not in hg19
https://www.ncbi.nlm.nih.gov/snp/
GRCh37.p13 chr 
 is located at  BP in GRCh37
Trying to liftover chr:- from hg18 to hg19:
chr:-
which corresponds to hg19

### Genotypes are presented in hg18 !
### Need to convert to hg19

### Liftover 

binary http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/

Hg18toHg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

Hg38toHg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

%%bash
suffix="map" # or bim
plink_file_prefix=""
chain_file="/home/borisov/software/hg38ToHg19.over.chain.gz" # /home/borisov/software/hg18ToHg19.over.chain.gz
# dos2unix "${pedfile_prefix}.${suffix}" # optional, check with vim
/home/borisov/software/liftOver \
<(awk 'OFS="\t"{print "chr",$1,$4,$4+1,$2}' \
"${plink_file_prefix}.${suffix}" | \
sed 's/\t//' | \
sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' | \
sed 's/chr25/chrXY/g' | sed 's/chr26/chrM/g') \
${chain_file} \
"${plink_file_prefix}"_hg19 \
"${plink_file_prefix}"_unMapped
if [ ${suffix} == map ]; then plink_input="file"; elif [ ${suffix} == bim ]; then plink_input="bfile"; fi
plink --${plink_input} "${plink_file_prefix}" \
--update-map <(awk '{print $4,$2}' "${plink_file_prefix}"_hg19) \
--exclude <(grep -v \# "${plink_file_prefix}"_unMapped | awk '{print $4}') \
--make-bed --out "${plink_file_prefix}_hg19"


for suffix in bed bim fam; do
    mv ${plink_file_prefix}.${suffix} ${plink_file_prefix}_NOThg19.${suffix} 
    mv ${plink_file_prefix}_hg19.${suffix} ${plink_file_prefix}.${suffix} 
done
