%%bash
# Checking relatedness using king
# from http://people.virginia.edu/~wc9c/KING/manual.html
# make sure that family ids are not identical (e.g., "0" for all samples), if needed, change FIDs to IIDs
%%bash
bfile=""
salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 \
srun king -b ${bfile}.bed \
--related \
--degree 2 \
--prefix ${bfile} \
--cpus 20

# if FIDs were changed to IIDs, change back famBK file
%%bash; bfile=""; cp ${bfile}.fam ${bfile}.famBK; awk '{print $2,$2,$3,$4,$5,$6}' ${bfile}.famBK > ${bfile}.fam
mv ${bfile}.famBK ${bfile}.fam 

%%bash
# remove relatives
bfile=""
plink --bfile ${bfile} \
--remove <(awk '{print $1,$2}' ${bfile}.kin0 | tail -n +2) \
--make-bed --out ${bfile}_norelated

%%bash
## if there are no relatives, copy files with "_norelated" suffix
bfile=""
cp ${bfile}.bed ${bfile}_norelated.bed
cp ${bfile}.bim ${bfile}_norelated.bim
cp ${bfile}.fam ${bfile}_norelated.fam
