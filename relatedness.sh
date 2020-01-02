# from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
# make sure that family ids are not identical (e.g., "0" for all samples)
king_relatedness () {
  bfile=$1
  king -b ${bfile}.bed \
  --related \
  --degree 2 \
  --prefix ${bfile} \
  --cpus 22
}
