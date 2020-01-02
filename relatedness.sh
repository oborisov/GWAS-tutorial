# from http://people.virginia.edu/~wc9c/KING/manual.html
# the script looks and filter out 2nd degree relatives
king_relatedness () {
  king -b ${bfile}.bed \
  --related \
  --degree 2 \
  --prefix ${bfile} \
  --cpus 22
}
