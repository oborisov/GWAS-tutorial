# from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
# make sure that family ids are not identical (e.g., "0" for all samples)
# if needed, change FIDs to IIDs
%%bash
bfile=""
cp ${bfile}.fam \
${bfile}.famBK
awk '{print $2,$2,$3,$4,$5,$6}' \
${bfile}.famBK > \
${bfile}.fam

%%bash
bfile=""
salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 \
srun king -b ${bfile}.bed \
--related \
--degree 2 \
--prefix ${bfile} \
--cpus 20

# if FIDs were changed to IIDs, change back famBK file
mv ${bfile}.famBK ${bfile}.fam 

# remove relatives
%%bash
bfile=""
plink --bfile ${bfile} \
--remove <(awk '{print $1,$2}' ${bfile}.kin0 | tail -n +2) \
--make-bed --out ${bfile}_norelated

## if there are no relatives, copy files with "_norelated" suffix
%%bash
bfile=""
cp ${bfile}.bed ${bfile}_norelated.bed
cp ${bfile}.bim ${bfile}_norelated.bim
cp ${bfile}.fam ${bfile}_norelated.fam
