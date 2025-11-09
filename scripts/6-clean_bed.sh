#########################################################################
# File Name: clean_bed.sh
# Last Modified: 2025-11-09
#########################################################################

input_cell=$1
names=`ls $input_cell | cut -f 1 -d '.'`
for name in $names
do
    bedtools sort -i $name.rawbed > $name.2.bed
    bedtools merge -i $name.2.bed -c 4 -o max > $name.3.bed
    mv $name.3.bed $name.bed
    rm $name.2.bed
done
