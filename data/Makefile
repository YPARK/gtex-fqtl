
all: $(foreach tis, $(shell cat tissues.txt), expr/$(tis).expr.txt.gz expr/$(tis).genes.txt.gz) \
  $(foreach tis, $(shell cat tissues.txt), covariates/$(tis).covariates.txt.gz) \
  genes.txt \
  individuals/matched.txt

EXPR := ../eqtl_data/GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices.tar.gz
COV := ../eqtl_data/GTEx_Analysis_2015-01-12_eQTLInputFiles_covariates.tar.gz 
GENE := ../eqtl_data/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt

genes.txt: $(GENE)
	cp $< $@

expr/%.expr.txt.gz: 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	tar xzvf $(EXPR) -C $(dir $@) $*.expr.txt
	gzip expr/$*.expr.txt

expr/%.genes.txt.gz: expr/%.expr.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | cut -f 1 | gzip > $@

covariates/%.covariates.txt.gz: 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	tar xzvf $(COV) -C $(dir $@) $*.covariates.txt
	gzip covariates/$*.covariates.txt

individuals/matched.txt: individuals/tissues.ind.txt individuals/geno.ind.txt
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./match.ind.R $^ $@

individuals/geno.ind.txt:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat ../genotypes/chr1.fam | awk '{ print $$1 }' | awk -F'-' '{ print $$1 FS $$2 }' > $@

individuals/tissues.ind.txt: $(foreach tis, $(shell cat tissues.txt), individuals/$(tis).ind.txt)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ > $@
	rm $^

individuals/%.ind.txt: expr/%.expr.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | head -n1 | cut -f2- | tr '\t' '\n' | awk -F'\t' '{ print $$0 FS NR FS "$*" }' > $@

