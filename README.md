# Factored QTL analysis on GTEx v6p data

Results of Sarkar and Park _et al._ (2017) submitted.

# Results

All the tissue and SNP effect sizes can be bound in

```
result/stat/chr1/50/combined.txt.gz

(...)

result/stat/chr22/50/combined.txt.gz
```

Note that individual by tissue gene expression
matrix `Y` was regressed on cis-regulatory genotype matrix `X` with
factored effect size matrix.

```
Y ~ X * sum_k (Theta[, k] * t(Theta[, k])) + confounders + errors
```

Each row contains 

1. `ENSEMBL.ID` : unique gene ID
2. `Chromosome` : chromosome name (`1` to `22`)
3. `TSS` : transcription start site (provided by GTEx v6p)
4. `Tissue.idx` : comma-separate tissue indexes
5. `Tissue.names` : comma-separate tissue names
6. `Tissue.theta` : comma-separate tissue effect sizes
7. `Tissue.se` : comma-separate tissue effect size standard errors
8. `Tissue.lodds` : comma-separate tissue PIP log-odds
9. `SNP.names` : comma-separate SNP names
10. `SNP.theta` : comma-separate SNP effect sizes
11. `SNP.se` : comma-separate SNP effect size standard errors
12. `SNP.lodds` : comma-separate SNP PIP log-odds
13. `k` : factor index
14. `pip` : posterior inclusion probability cutoff

