# from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
# make sure that family ids are not identical (e.g., "0" for all samples)
# if needed, change FIDs to IIDs
cp ${bfile}.fam \
${bfile}.famBK
awk '{print $2,$2,$3,$4,$5,$6}' \
${bfile}.famBK > \
${bfile}.fam

%%bash
king_relatedness () {
  bfile=$1
  king -b ${bfile}.bed \
  --related \
  --degree 2 \
  --prefix ${bfile} \
  --cpus 22
}
bfile="/home/borisov/LUTO/Dutch/dcg_checkedsex_geno02_mind02_geno002_mind002"
king_relatedness ${bfile}

## if there are no relatives, copy files with "_norelated" suffix
cp ${bfile}.bed ${bfile}_norelated.bed
cp ${bfile}.bim ${bfile}_norelated.bim
cp ${bfile}.fam ${bfile}_norelated.fam
