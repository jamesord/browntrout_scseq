# script to extract descriptions from protein-level PANNZER homology search and generate gene-level output
# James Ord 25/07/23

cut -f1 DE.out | grep "CAAJ" | awk 'BEGIN{FS=OFS="."} NF--' > gnames1.temp
cut -f1 DE.out | grep -v "CAAJ" | tail -n +2 | awk 'BEGIN{FS=OFS="."} NF--' > gnames2.temp
cat gnames1.temp gnames2.temp > gnames3.temp
echo 'gid' | cat - gnames3.temp > gnames4.temp
cut -f9 DE.out > DE2.temp
paste gnames4.temp DE2.temp | sort | uniq > DE_genelevel.temp
cut -f1 DE_genelevel.temp | uniq > gnames_uniq.temp

for i in `cat gnames_uniq.temp`;do
fgrep -w ${i} DE_genelevel.temp | cut -f2 | paste -s -d";" - > descs${i}.temp
echo ${i} | paste - descs${i}.temp > descs2${i}.temp
done

cat descs2*.temp > DE_genelevel.txt

rm *.temp