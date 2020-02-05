## Converting the chip

Strand: TOP  

Strand files:
Script: https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh  
TOP: https://www.well.ox.ac.uk/~wrayner/strand/index.html  
SRC: https://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html  

%%bash
bfile=""
strand_file="GSAMD-24v1-0_20011747_A1"
bash /home/borisov/software/update_build.sh \
${bfile} \
/home/borisov/software/${strand_file}-b37.strand \
${bfile}_SRC
