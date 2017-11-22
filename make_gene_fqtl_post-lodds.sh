#!/bin/bash -l

if [ $# -lt 1 ]; then
    exit 1
fi

fqtl_hdr=$1
snp_out=$fqtl_hdr.snp-null.txt.gz
tis_out=$fqtl_hdr.tis-null.txt.gz

[ -f $snp_out ] || \
    cat $fqtl_hdr.snp.lodds.txt.gz | gzip -d | \
    awk -F'\t' '(NR == 1) { for(j=1; j<=NF; ++j) maxval[j] = $j; }
NR > 1 { for(j=1; j<=NF; ++j) if($j>maxval[j]) maxval[j] = $j; }
END { for(j in maxval) print j FS maxval[j] }' | sort -k1n | \
    gzip > $snp_out

[ -f $tis_out ] || \
    cat $fqtl_hdr.tis.lodds.txt.gz | gzip -d | \
    awk -F'\t' '(NR == 1) { for(j=1; j<=NF; ++j) maxval[j] = $j; }
NR > 1 { for(j=1; j<=NF; ++j) if($j>maxval[j]) maxval[j] = $j; }
END { for(j in maxval) print j FS maxval[j] }' | sort -k1n | \
    gzip > $tis_out
