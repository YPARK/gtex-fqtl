#!/bin/bash -l

# Some methods do not know how to handle missing values.
# So, we just subset simulated full data

if [ $# -lt 3 ]; then
cat <<EOF

datahdr=\$1
n=\$2
rseed=\$3

EOF

    exit 1
fi

printf "[%s] take subsets\n\n" "$(date)"

datahdr=$1
n=$2
rseed=$3

cat <<EOF

datahdr=$1
n=$2
rseed=$3

EOF

get_seeded_random () {
    seed="$1"
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
        </dev/zero 2>/dev/null
}

[ -f $datahdr.fam ] || exit 1

ntot=$(cat $datahdr.fam| wc -l)

[ ${ntot} -ge $n ] || exit 1

cat $datahdr.fam | awk '{ print NR FS $0  }' | \
    shuf -n${n} --random-source=<(get_seeded_random ${rseed}) | \
    sort -k1n > $datahdr.random.$n.ind

cat $datahdr.random.$n.ind | cut -d' ' -f 2- > $datahdr.random.$n.fam

./bin/plink --bfile $datahdr --make-bed \
    --keep-fam $datahdr.random.$n.fam \
    --out $datahdr.random.$n || \
    exit 1

zcat $datahdr.yfull.txt.gz | \
    awk -vROWS=$(cat $datahdr.random.$n.ind | awk '{ if(NR > 1) printf ","; printf $1 }') \
    -f util_subset_rows.awk | \
    gzip > $datahdr.random.$n.y.txt.gz

[ -f $datahdr.random.$n.y.txt.gz ] || exit 1

[ $(zcat $datahdr.random.$n.y.txt.gz | wc -l) -eq ${n} ] || exit 1

printf "[%s] subsets taken n = %d.\n\n" "$(date)" ${n}

