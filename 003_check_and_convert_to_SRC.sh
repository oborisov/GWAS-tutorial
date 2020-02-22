### Checking (and converting) chip strand
http://mccarthy.well.ox.ac.uk/chipendium/ui/
GSAMD-24v1-0_20011747_A1
Strand alignment TOP
Probability of strand match 100 

%%bash
bfile=""
strand_file="GSAMD-24v1-0_20011747_A1"
strand_file=${strand_file}-b37 # add .Source here if SRC
bash /home/borisov/software/update_build.sh \
${bfile} \
/home/borisov/software/manifest_strand/${strand_file}.strand \
${bfile}_SRC
