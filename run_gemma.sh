#!/bin/bash -l

if [ $# -lt 4 ]; then
    exit 1
fi

BED=$1        # e.g., BED=scratch/temp
Y=$2          # e.g., Y=scratch/temp.y.txt.gz
output=$3     # e.g., temp.bslmm.gz
pve_output=$4 # e.g., PVE output

[ -f $BED.bed ] || exit 1
[ -f $Y ] || exit 1
[ -f $output ] || exit 1
[ -f $pve_output ] || exit 1

tempdir=$5

echo "Start $0 $@"

if [ -z $tempdir ]; then
    outdir=$(dirname $output)
    [ -d $outdir ] || mkdir -p $outdir
    tempdir=$(mktemp -d "$outdir/temp.XXXXXX")
else
    [ -d $tempdir ] || mkdir -p $tempdir
fi

# 1. construct temporary BED for each column of Y
ncol=$(zcat $Y | head -n 1 | wc -w)

for((j=1; j<=$ncol; ++j)); do

    paste -d' ' <(cat $BED.fam | awk '{ print $1 FS $2 }') \
        <(zcat $Y | cut -f $j) > $tempdir/pheno.${j}.txt
    
    ./bin/plink --bfile $BED --pheno $tempdir/pheno.${j}.txt --out $tempdir/data.${j} \
        --make-bed --keep-allele-order
done

# 2. run GEMMA within working directory

back=$(pwd)
gemma=$back/bin/gemma

burin=10000
mcmc=100000

[ -f $output ] && rm -f $output
[ -f $pve_output ] && rm -f $pve_output

printf "chr\trs\tloc\tn.mis\tdense\tsparse\tpip\ttis\n" | gzip > $output

for((j=1; j<=$ncol; ++j)); do

    cd $tempdir

    $gemma -bfile data.${j} -bslmm 1 -o bslmm1.${j} -w $burin -s $mcmc || exit 1

    cd ${back}

    Rscript util.col.stat.R $tempdir/output/bslmm1.${j}.hyp.txt \
        $tempdir/output/bslmm1.${j}.summary.txt || exit 1

    cat $tempdir/output/bslmm1.${j}.summary.txt | \
        awk -vtis=${j} -F'\t' '{ print $0 FS tis }' | \
        gzip >> $pve_output

    cat $tempdir/output/bslmm1.${j}.param.txt | \
        awk -F'\t' -vtis=${j} 'NR > 1 { print $0 FS tis }' | \
        gzip >> $output

    printf "Finished : $j / ${ncol}\n\n\n"

done

[ -d $tempdir ] && rm -r $tempdir

[ -f ${output} ] || exit 1

printf "successfully copied results: $0 $@ --> $output and $pve_output"
