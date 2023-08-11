# script to extract GO-terms from protein-level PANNZER homology search and generate gene-level output
# James Ord 25/07/23

cut -f1 GO.out | grep "CAAJ" | awk 'BEGIN{FS=OFS="."} NF--' > gnames1.temp
cut -f1 GO.out | grep -v "CAAJ" | tail -n +2 | awk 'BEGIN{FS=OFS="."} NF--' > gnames2.temp
cat gnames1.temp gnames2.temp > gnames3.temp
echo 'gid' | cat - gnames3.temp > gnames4.temp
cut -f2,3,4 GO.out > GO2.temp
paste gnames4.temp GO2.temp | sort | uniq > GO_genelevel.txt
rm *.temp