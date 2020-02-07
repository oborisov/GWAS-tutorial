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

# script to update fids
%%bash
bfile="/home/borisov/LUTO/Dutch/dcg_checkedsex_geno02_mind02_geno002_mind002"
awk '{print $1,$2,$2,$2}' ${bfile}.fam > ${bfile}_update_ids_temp
plink --bfile ${bfile} \
--update-ids ${bfile}_update_ids_temp \
--make-bed --out ${bfile}_updids
rm ${bfile}_update_ids_temp

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
