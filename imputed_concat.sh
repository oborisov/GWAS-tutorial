%%bash
## Assign path to bfile
bfile=
## Reassign default variables if run in other account / on other machine
dir1000genomes=/home/borisov/software/1000GP_Phase3/
chunk_size=5000000
imputed_concat () {
    bfile=$1
    chr=$2
    chunk_size=$3
    dir1000genomes=$4
    chromosome_min_coord=`tail -n +2 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | head -1 | cut -d " " -f 1`
    chromosome_max_coord=`tail -n -1 ${dir1000genomes}/genetic_map_chr${chr}_combined_b37.txt | cut -d " " -f 1`
    # concatenate .gen imputed files to whole chromosome
    for lower in `seq ${chromosome_min_coord} ${chunk_size} ${chromosome_max_coord}`; do
        upper=`echo ${lower} + ${chunk_size} | bc`
        echo "processing chunk called ${bfile}_ref_chr${chr}.phased.impute2.${upper}"
        cat ${bfile}_ref_chr${chr}.phased.impute2.${upper} >> \
        ${bfile}_ref_chr${chr}.phased.impute2
        # concatenate _info files to whole chromosome
        if [ ! -f ${bfile}_ref_chr${chr}.phased.impute2_info ]; then
            head -1 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info > \
            ${bfile}_ref_chr${chr}.phased.impute2_info
        fi
        tail -n +2 ${bfile}_ref_chr${chr}.phased.impute2.${upper}_info >> \
        ${bfile}_ref_chr${chr}.phased.impute2_info
        # remove temporary chunks files
#        rm ${bfile}_ref_chr${chr}.phased.impute2.${upper}*
    done
    # compress .gen to .gen.gz
    echo "gzipping ${bfile}_ref_chr${chr}.phased.impute2"
    gzip -f ${bfile}_ref_chr${chr}.phased.impute2_info &
    gzip -f ${bfile}_ref_chr${chr}.phased.impute2 &

}

for chr in {1..22}; do
    rm ${bfile}_ref_chr${chr}.phased.impute2 2>/dev/null
    rm ${bfile}_ref_chr${chr}.phased.impute2_info 2>/dev/null
    echo "concatenating chromosome ${chr}"
        imputed_concat \
        ${bfile} \
        ${chr} \
        ${chunk_size} \
        ${dir1000genomes} &
done
wait
echo "concatenation is done"
