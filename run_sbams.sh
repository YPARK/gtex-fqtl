#!/bin/bash -l

datahdr=scratch/temp.random.123
tempdir=scratch/temp.sbams

[ -d $tempdir ] || mkdir -p $tempdir

tempdata=$(mktemp "$tempdir/temp.sbams.XXXXX")

# 1. data formatting
n=$(zcat $datahdr.y.txt.gz | wc -l)
m=$(zcat $datahdr.y.txt.gz | head -n1 | wc -w)

[ -f $tempdata.mcmc.dat ] && rm -f $tempdata.mcmc.dat

R --vanilla <<EOF

plink.hdr = '${datahdr}'
y.file = '${datahdr}.y.txt.gz'

options(stringsAsFactors = FALSE)
library(fqtl)

y = read.table(y.file)
plink = read.plink(plink.hdr)

xx = data.frame(type = 'covariate', name = plink[['BIM']][, 2],
                data = t(scale(plink[['BED']], center = TRUE, scale = FALSE)))

yy = data.frame(type = 'response', name = 1:dim(y)[2],
                data = t(scale(y)))

out = rbind(yy, xx)

write.table(out, file = '$tempdata.mcmc.dat', row.names = FALSE, 
col.names = FALSE, quote = FALSE)
EOF

# 2. simply re-use hyper-parameter settings of Wen (2014)

cat > $tempdata.hyp << EOF
0.05 0.20
0.10 0.40
0.20 0.80
0.40 1.60
EOF

# 3. We cannot run this method for large number of tissues

./bin/sbams_mvlr -d $tempdata.mcmc.dat \
    -g $tempdata.hyp \
    -mcmc -b 10 -r 10 \
    -o $tempdata.output.txt
    
