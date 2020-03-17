%%bash
# Checking relatedness using king
# from http://people.virginia.edu/~wc9c/KING/manual.html
# make sure that family ids are not identical (e.g., "0" for all samples), if needed, change FIDs to IIDs
bfile=""
salloc --mem=16000M --time=5:00:00 --cpus-per-task=20 \
srun king -b ${bfile}.bed \
--related \
--degree 3 \
--prefix ${bfile} \
--cpus 20

%%bash
# call rate
bfile=""
plink --bfile ${bfile} \
--keep <(cat <(awk '{print $1,$2}' ${bfile}.kin0 | tail -n +2) <(awk '{print $3,$4}' ${bfile}.kin0 | tail -n +2) ) \
--missing --out ${bfile}_CR
cat ${bfile}_CR.imiss



%%bash
# remove relatives
bfile=""
plink --bfile ${bfile} \
--remove <(awk '{print $1,$2}' ${bfile}.kin0 | tail -n +2) \
--make-bed --out ${bfile}_norelated
cat ${bfile}.kin0

%%bash
## if there are no relatives, copy files with "_norelated" suffix
bfile=""
plink --bfile ${bfile} \
--make-bed --out ${bfile}_norelated

%%bash
# script to update fids
bfile=""
awk '{print $1,$2,$2,$2}' ${bfile}.fam > ${bfile}_update_ids_temp
plink --bfile ${bfile} \
--update-ids ${bfile}_update_ids_temp \
--make-bed --out ${bfile}_updids
rm ${bfile}_update_ids_temp

