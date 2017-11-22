#!/bin/bash -l

if [ $# -lt 7 ]; then
cat <<EOF
Required arguments:

geneId=\$1
rank=\$2
h2=\$3
n_causal_snp=\$4
n_causal_tis=\$5
rseed=\$6
tempdata=\$7
EOF

    exit 1
fi

geneId=$1
rank=$2
h2=$3
n_causal_snp=$4
n_causal_tis=$5
rseed=$6
tempdata=$7

[ -d $(dirname $tempdata) ] || mkdir -p $(dirname $tempdata)

cisk=1000
genefile=data/genes.txt 

geneInfo=$(zless $genefile | tail -n+2 | head -n ${geneId} | tail -n1)

printf "$geneInfo\n" > $tempdata.gene

ensg=$(echo $geneInfo | awk '{ print $1 }')
chr=$(echo $geneInfo | awk '{ print $2 }')
tss=$(echo $geneInfo | awk '{ print $3 }')

lb=$(echo $tss | awk -vcis=${cisk}000 '{ print ($1 > cis) ? ($1 - cis) : 0; }')
ub=$(echo $tss | awk -vcis=${cisk}000 '{ print ($1 + cis); }')

## Subset genotype matrix
if [ -f $tempdata.bed ]; then
    rm $tempdata.bed $tempdata.bim $tempdata.fam
    echo "Warning: overwritting PLINK files"
fi

./bin/plink --bfile genotypes/chr${chr} --make-bed \
    --chr $chr --from-bp $lb --to-bp $ub \
    --out $tempdata

# ./bin/plink --bfile $tempdata --recode vcf --out $tempdata
# gzip $tempdata.vcf

[ -f $tempdata.bed ] || exit 1

## Generate simulation data
./sim.data.R $rank $h2 $n_causal_snp $n_causal_tis $rseed $tempdata

[ -f $tempdata.yfull.txt.gz ] || exit 1
[ -f $tempdata.y.txt.gz ] || exit 1
[ -f $tempdata.eta.txt.gz ] || exit 1
[ -f $tempdata.causal.txt.gz ] || exit 1
[ -f $tempdata.settings.txt.gz ] || exit 1

echo "Successfully generated : $tempdata"
