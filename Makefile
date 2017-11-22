all:

################################################################
DATA: data/tissues.txt
	$(MAKE) -C data

EXPR := eqtl_data/GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices.tar.gz

CHR := $(shell seq 1 22)

data/tissues.txt:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	tar tzf $(EXPR) | awk '{ gsub(".expr.txt", "", $$0); print $$0 }' > $@

# 1. Take GTEx v6 data

FQTL-STEP0: data/fqtl.genes.txt jobs/fqtl/data.job.txt.gz scratch/data/covariates

# create matched covariates
scratch/data/covariates:
	Rscript --vanilla make.covariates.R $@

data/fqtl.genes.txt: data/genes.txt
	cat $< | tail -n +2 | awk '$$2 ~ /^[0-9]+$$/' | sort -k2n -k3n | awk -F'\t' '{ print NR FS $$0 }' > $@

jobs/fqtl/data.job.txt.gz: $(foreach g, $(shell cat data/fqtl.genes.txt | cut -f1), jobs/fqtl/data-$(g)-job)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	for f in jobs/fqtl/data-*-job; do cat $$f | gzip >> $@; done
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N FQTL.DATA.$* -binding "linear:1" -q short -l h_vmem=1g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/fqtl/data-%-job:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d scratch/data/$*/ ] || mkdir -p scratch/data/$*/
	echo ./make_gene_expr_data.sh $* scratch/data/$*/ > $@


FQTL-STEP1: $(foreach chr, $(CHR), data/fqtl-$(chr)-valid.genes.txt)

data/fqtl-%-valid.genes.txt: data/fqtl.genes.txt
	cat data/fqtl.genes.txt | awk -vChr=$* -F'\t' '$$3 == Chr { if(system("[ -f scratch/data/" $$1 "/tissues.txt.gz ]") == 0) { "zcat scratch/data/" $$1 "/tissues.txt.gz | wc -l" | getline ntis; print $$0 FS ntis; } }' > $@

################################################################
# 2. Confounder correction with sparse regression and fit factored QTL
# models.
FQTL-STEP2: $(foreach chr, $(CHR), jobs/fqtl/resid-$(chr)-job.txt.gz)

FQTL-STEP2-resubmit: $(foreach chr, $(CHR), jobs/fqtl/resid-$(chr)-job-resubmit.txt.gz)

FQTL-STEP2-clean: 
	rm -f $(foreach chr, $(CHR), jobs/fqtl/resid-$(chr)-job.txt.gz)
	rm -f $(foreach chr, $(CHR), jobs/fqtl/resid-$(chr)-job-resubmit.txt.gz)

# for genes with mininum 10 tissues
jobs/fqtl/resid-%-job.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" | gzip > $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vRDIR=result/fqtl 'system("mkdir -p " RDIR "/" $$1 "; [ ! -f " RDIR "/" $$1 "/fqtl.resid.txt.gz ]") == 0 { print "./make.gene.resid.R" FS $$1 FS "FALSE" FS (RDIR "/" $$1 "/fqtl.resid.txt.gz") FS (RDIR "/" $$1 "/fqtl.pve.txt.gz") }' | gzip > $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vRDIR=result/fqtl-std 'system("mkdir -p " RDIR "/" $$1 "; [ ! -f " RDIR "/" $$1 "/fqtl.resid.txt.gz ]" ) == 0 { print "./make.gene.resid.R" FS $$1 FS "TRUE" FS (RDIR "/" $$1 "/fqtl.resid.txt.gz") FS (RDIR "/" $$1 "/fqtl.pve.txt.gz") }' | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Resid.$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/fqtl/resid-%-job-resubmit.txt.gz: jobs/fqtl/resid-%-job.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Resid.$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

################################################################
# 3. Fit FQTL models
FQTL-STEP3: $(foreach chr, $(CHR), jobs/fqtl/fit-$(chr)-job.txt.gz)

FQTL-STEP3-resubmit: $(foreach chr, $(CHR), jobs/fqtl/fit-$(chr)-job-resubmit.txt.gz)

# for genes with mininum 10 tissues
FACTORS := 1 10 50
Rseed := 1001

jobs/fqtl/fit-%-job.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vFACTORS="$(FACTORS)" -vRDIR=result/fqtl-std/ 'BEGIN{ split(FACTORS, factors); } { rdir = (RDIR "/" $$1 "/"); if(system("[ -f " rdir "/fqtl.resid.txt.gz ]") == 0) for(f in factors) { k = factors[f]; system("mkdir -p " rdir "/" k "/"); if(system("[ -f " rdir "/" k "/fqtl.tis.lodds.txt.gz ]") == 0) break; print "./make.gene.fqtl.R" FS $$1 FS k FS "TRUE" FS (rdir "/fqtl.resid.txt.gz") FS (rdir "/" k "/" fqtl) } }' | gzip >> $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vFACTORS="$(FACTORS)" -vRDIR=result/fqtl-perm/ -vDDIR=result/fqtl-std/ 'BEGIN{ split(FACTORS, factors); } { rdir = (RDIR "/" $$1 "/"); ddir = (DDIR "/" $$1 "/"); if(system("[ -f " ddir "/fqtl.resid.txt.gz ]") == 0) for(f in factors) { k = factors[f]; system("mkdir -p " rdir "/" k "/"); if(system("[ -f " rdir "/" k "/fqtl.tis.lodds.txt.gz ]") == 0) break; print "./make.gene.fqtl.R" FS $$1 FS k FS "TRUE" FS (ddir "/fqtl.resid.txt.gz") FS (rdir "/" k "/" fqtl) FS $(Rseed) } }' | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N FQTL.$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/fqtl/fit-%-job-resubmit.txt.gz: jobs/fqtl/fit-%-job.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'system("[ ! -f " $$NF ".tis.lodds.txt.gz ]") == 0 && system("[ -f " $$(NF-1) " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N FQTL.$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

queue-fqtl-step3:
	for j in jobs/fqtl/fit-*-job.txt.gz; do [ $$(zcat $$j | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $$j | wc -l) -N FQTL.$$(basename $$j) -binding "linear:1" -q long -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $$j; done

# check NA results
check-%-fqtl:
	for g in $$(cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }'); do for k in $(FACTORS); do if [ -f result/fqtl/$$g/$$k/fqtl.tis.lodds.txt.gz ] && [ $$(zcat result/fqtl/$$g/$$k/fqtl.tis.lodds.txt.gz | awk '/NA/' | wc -l) -gt 0 ] ; then printf "%s deleted!\n\n" "result/fqtl/$$g/$$k/fqtl.*.txt.gz"; rm result/fqtl/$$g/$$k/fqtl.*.txt.gz; fi; done; done


################################################################
# 4. Post-processing -- compute basic statistics and combine results
FQTL-STEP4: $(foreach chr, $(CHR), jobs/fqtl/stat-$(chr)-job.txt.gz)

jobs/fqtl/stat-%-job.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vFACTORS="$(FACTORS)" -vRDIR=result/fqtl-std/ 'BEGIN{ split(FACTORS, factors); } { rdir = (RDIR "/" $$1 "/"); if(system("[ -f " rdir "/fqtl.resid.txt.gz ]") == 0) for (f in factors) { k = factors[f]; system("mkdir -p " rdir "/" k "/"); hdr = rdir "/" k "/fqtl"; if(system("[ -f " hdr ".tis.lodds.txt.gz ]") == 0)  { print "./make.gene.fqtl.post.R" FS hdr FS ("scratch/data/" $$1) FS (hdr ".s1.txt.gz") FS (hdr ".s2.txt.gz") FS (hdr ".combined.txt.gz"); } } }' | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip >> $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vFACTORS="$(FACTORS)" -vRDIR=result/fqtl-perm/ -vDDIR=result/fqtl-std/ 'BEGIN{ split(FACTORS, factors); } { rdir = (RDIR "/" $$1 "/"); if(system("[ -f " DDIR "/" $$1 "/fqtl.resid.txt.gz ]") == 0) for (f in factors) { k = factors[f]; system("mkdir -p " rdir "/" k "/"); hdr = rdir "/" k "/fqtl"; if(system("[ -f " hdr ".tis.lodds.txt.gz ]") == 0)  { print "./make_gene_fqtl_post-lodds.sh" FS hdr; } } }' | awk 'system("[ ! -f " $$NF ".tis-null.txt.gz ]") == 0' | gzip >> $@
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vFACTORS="$(FACTORS)" -vRDIR=result/fqtl-std/ -vDDIR=result/fqtl-std/ 'BEGIN{ split(FACTORS, factors); } { rdir = (RDIR "/" $$1 "/"); if(system("[ -f " DDIR "/" $$1 "/fqtl.resid.txt.gz ]") == 0) for (f in factors) { k = factors[f]; system("mkdir -p " rdir "/" k "/"); hdr = rdir "/" k "/fqtl"; if(system("[ -f " hdr ".tis.lodds.txt.gz ]") == 0)  { print "./make_gene_fqtl_post-lodds.sh" FS hdr; } } }' | awk 'system("[ ! -f " $$NF ".tis-null.txt.gz ]") == 0' | gzip >> $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N STAT.$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@


################################################################
# 5. Merge them all
FQTL-STEP5: $(foreach chr, $(CHR), $(foreach k, $(FACTORS), result/stat/chr$(chr)/$(k)/s1.txt.gz result/stat/chr$(chr)/$(k)/s2.txt.gz result/stat/chr$(chr)/$(k)/combined.txt.gz result/stat/chr$(chr)/$(k)/snp-null.txt.gz result/stat/chr$(chr)/$(k)/tis-null.txt.gz result/perm/chr$(chr)/$(k)/snp-null.txt.gz result/perm/chr$(chr)/$(k)/tis-null.txt.gz))

# % = chr/rank/stat-type
result/stat/chr%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	cat data/fqtl-$(shell echo $* | awk -F'/' '{ print $$1 }')-valid.genes.txt | awk  -vRANK=$(shell echo $* | awk -F'/' '{ print $$2 }') -vSTAT=$(shell echo $* | awk -F'/' '{ print $$3 }') '$$NF >= 10 { print "result/fqtl-std/" $$1 "/" RANK "/fqtl." STAT ".txt.gz" }' | awk 'system("[ -f " $$1 " ]") == 0 { system("cat " $$1 " >> $@") }'

result/perm/chr%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	cat data/fqtl-$(shell echo $* | awk -F'/' '{ print $$1 }')-valid.genes.txt | awk  -vRANK=$(shell echo $* | awk -F'/' '{ print $$2 }') -vSTAT=$(shell echo $* | awk -F'/' '{ print $$3 }') '$$NF >= 10 { print "result/fqtl-perm/" $$1 "/" RANK "/fqtl." STAT ".txt.gz" }' | awk 'system("[ -f " $$1 " ]") == 0 { system("cat " $$1 " >> $@") }'


################################################################
# 6. Show basic statistics
FIG := result/figures

FQTL-STEP6: $(foreach k, $(FACTORS), $(FIG)/fdr-$(k).pdf \
  $(foreach pip, 50 90 95, $(FIG)/factor-s2-$(k)-$(pip).pdf) \
  $(foreach pip, 50 90 95, $(FIG)/tis-pair-$(k)-$(pip).pdf))

$(FIG)/fdr-%.pdf: ./figure4.fdr.R
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./figure4.fdr.R $* $@

# % = $(k)-$(pip)
$(FIG)/factor-s2-%.pdf: ./figure5.factor.stat.R
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./figure5.factor.stat.R $(shell echo $* | sed 's/-/ /g') $(foreach s, 1 2, $(FIG)/factor-s$(s)-$*.pdf)

# % = $(k)-$(pip)
$(FIG)/tis-pair-%.pdf: ./figure6.tissue.stat.R
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./figure6.tissue.stat.R $(shell echo $* | sed 's/-/ /g') $(foreach s, hist_g hist_t pair, $(FIG)/tis-$(s)-$*.pdf)


######################
# summary-based TWAS #
######################

FQTL-STEP7: $(foreach chr, $(CHR), jobs/fqtl-stwas/stwas-IGAP-$(chr)-job.txt.gz)

FQTL-STEP7-long: $(foreach chr, $(CHR), jobs/fqtl-stwas/stwas-IGAP-$(chr)-job-long.txt.gz)

# for genes with mininum 10 tissues
STWAS_FACTORS := 50

# % = IGAP-chr21
jobs/fqtl-stwas/stwas-%-job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p result/stwas/fqtl-std/$(shell echo $* | sed 's/-/\//')
	cat data/fqtl-$(shell echo $* | awk -F'-' '{ print $$2 }')-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vGWAS=$(shell echo $* | awk -F'-' '{ print $$1 }') -vDIR=result/stwas/fqtl-std/$(shell echo $* | sed 's/-/\//') 'system("[ -f result/fqtl-std/" $$1 "/$(STWAS_FACTORS)/fqtl.resid.txt.gz ]") == 0 { print "./make.gene.fqtl.stwas.R result/fqtl-std/" $$1 "/$(STWAS_FACTORS)/fqtl" FS "scratch/data/" $$1 FS GWAS FS DIR "/" $$1 ".stwas.gz" }' | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N STWAS-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/fqtl-stwas/stwas-%-job-long.txt.gz: jobs/fqtl-stwas/stwas-%-job.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0 && ' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N STWAS-$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@


##############################################
# Comparison with other methods on real data #
##############################################

################################################################

##############
# sparse QTL #
##############

SQTL-STEP1: $(foreach chr, $(CHR), jobs/sqtl/sqtl-$(chr)-job.txt.gz)

SQTL-STEP1-resubmit: $(foreach chr, $(CHR), jobs/sqtl/sqtl-$(chr)-job-resubmit.txt.gz)

# % = $(chr)
jobs/sqtl/sqtl-%-job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p result/sqtl/$*
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vRDIR=result/sqtl/$* '{ print "./make.gene.sqtl.R" FS $$1 FS "TRUE" FS ("result/fqtl-std/" $$1 "/fqtl.resid.txt.gz") FS (RDIR "/" $$1 "/sqtl") }' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N SQTL-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/sqtl/sqtl-%-job-resubmit.txt.gz: jobs/sqtl/sqtl-%-job.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF ".effect.lodds.txt.gz ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N SQTL-$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

SQTL-STEP2: $(foreach chr, $(CHR), jobs/sqtl-stwas/sqtl-IGAP-$(chr)-job.txt.gz)

SQTL-STEP2-long: $(foreach chr, $(CHR), jobs/sqtl-stwas/sqtl-IGAP-$(chr)-job-long.txt.gz)

# % = IGAP-$(chr)
jobs/sqtl-stwas/sqtl-%-job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p result/stwas/sqtl/$(shell echo $* | sed 's/-/\//')
	cat data/fqtl-$(shell echo $* | awk -F'-' '{ print $$2 }')-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vGWAS=$(shell echo $* | awk -F'-' '{ print $$1 }') -vRDIR=result/sqtl/$(shell echo $* | awk -F'-' '{ print $$2 }') -vDIR=result/stwas/sqtl/$(shell echo $* | sed 's/-/\//') '{ print "./make.gene.sqtl.stwas.R " FS (RDIR "/" $$1 "/sqtl") FS ("scratch/data/" $$1 "/") FS "IGAP" FS (DIR "/" $$1 ".stwas.gz") }' | gzip > $@
#@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N SST-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/sqtl-stwas/sqtl-%-job-long.txt.gz: jobs/sqtl-stwas/sqtl-%-job.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N SST-$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@





################################################################

#########
# BSLMM #
#########

BSLMM-STEP1: $(foreach chr, $(CHR), jobs/bslmm/bslmm-$(chr)-job.txt.gz)

BSLMM-STEP1-check: $(foreach chr, $(CHR), check/bslmm-$(chr)-check.txt.gz)

BSLMM-STEP1-resubmit: $(foreach chr, $(CHR), jobs/bslmm/bslmm-$(chr)-job-resubmit.txt.gz)

# % = $(chr)
jobs/bslmm/bslmm-%-job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p result/bslmm/$*
	@mkdir -p /broad/hptmp/ypp/bslmm/$*
	cat data/fqtl-$*-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vRDIR=result/bslmm/$* -vTEMP=/broad/hptmp/ypp/bslmm/$* '{ print "./run_gemma.sh scratch/data/" $$1 "/plink" FS ("result/fqtl-std/" $$1 "/fqtl.resid.txt.gz") FS (RDIR "/bslmm-" $$1 "-weights.gz") FS (RDIR "/bslmm-" $$1 "-pve.gz") FS (TEMP "/" $$1) }' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N BSLMM-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/bslmm/bslmm-%-job-resubmit.txt.gz: jobs/bslmm/bslmm-%-job.txt.gz
	zcat $< | awk 'system("[ -f " $$3 " ]") == 0 && system("[ ! -f " $$(NF - 2) " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N BSLMM-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

check/bslmm-%-check.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	ls -1 result/bslmm/$*/bslmm-*-weights.gz | xargs -I file sh -c "[ \$$(zcat file | wc -l) -eq 1 ] && rm file && echo removed file" | gzip > $@

BSLMM-STEP2: $(foreach chr, $(CHR), jobs/bslmm-stwas/bslmm-IGAP-$(chr)-job.txt.gz)

BSLMM-STEP2-long: $(foreach chr, $(CHR), jobs/bslmm-stwas/bslmm-IGAP-$(chr)-job-long.txt.gz)

# % = IGAP-$(chr)
jobs/bslmm-stwas/bslmm-%-job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p result/stwas/bslmm/$(shell echo $* | sed 's/-/\//')
	cat data/fqtl-$(shell echo $* | awk -F'-' '{ print $$2 }')-valid.genes.txt | awk '$$NF >= 10 { print $$1 }' | awk -vGWAS=$(shell echo $* | awk -F'-' '{ print $$1 }') -vRDIR=result/bslmm/$(shell echo $* | awk -F'-' '{ print $$2 }') -vDIR=result/stwas/bslmm/$(shell echo $* | sed 's/-/\//') '{ print "./make.gene.bslmm.stwas.R " FS (RDIR "/bslmm-" $$1 "-weights.gz") FS ("scratch/data/" $$1 "/") FS "IGAP" FS (DIR "/" $$1 ".stwas.gz") }' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N BST-$* -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/bslmm-stwas/bslmm-%-job-long.txt.gz: jobs/bslmm-stwas/bslmm-%-job.txt.gz
	zcat $< | awk 'system("[ ! -f " $$NF " ]") == 0' | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N BST-$* -binding "linear:1" -q long -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@


################################################################

MERGE-STWAS: $(foreach m, bslmm sqtl fqtl-std, $(foreach chr, $(CHR), result/stat/stwas/$(m)/IGAP/$(chr).stwas.combined.txt.gz)) $(foreach chr, $(CHR), result/stat/stwas/fqtl-std/IGAP/chr$(chr).stwas.matched.txt.gz)

# % = $(method)/$(gwas)/$(chr)
result/stat/stwas/%.stwas.combined.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	ls -1 result/stwas/$*/*.stwas.gz | xargs cat >> $@

# % = $(chr)
result/stat/stwas/fqtl-std/IGAP/chr%.stwas.matched.txt.gz: result/stat/stwas/fqtl-std/IGAP/%.stwas.combined.txt.gz result/stat/chr%/50/combined.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.matched.stwas.R $^ $@


###############
# simulations #
###############

################################################################
# randomly select 100 genes
result/simulation/selected.genes:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@shuf -i1-$(shell cat data/genes.txt | wc -l) | head -n100 > $@

simdata := scratch/simulation/data
glmnet := scratch/simulation/result/glmnet
phenix := scratch/simulation/result/phenix
bslmm := scratch/simulation/result/bslmm
sqtl := scratch/simulation/result/sqtl
fqtl := scratch/simulation/result/fqtl
metasoft := scratch/simulation/result/metasoft
mashr := scratch/simulation/result/mashr
zfqtl := scratch/simulation/result/zfqtl

Rank := 1 2 3
H2 := 0.05 0.1 0.15 0.2 0.3
Ncausal := 1 3 5
Ntissue := 5 10 20 30
Rseed := 1001 # 1003 1667

SIM-DATA-JOB: result/simulation/selected.genes $(foreach g, $(shell cat result/simulation/selected.genes), jobs/data-$(g)-jobs.txt.gz)
	echo $@

SIM-JOB: $(foreach g, $(shell cat result/simulation/selected.genes), jobs/sim-$(g)-jobs.txt.gz)
	echo $@

SIM-QUEUE-JOB: $(foreach g, $(shell cat result/simulation/selected.genes), jobs/sim-$(g)-jobs.txt.gz)
	for f in $^; do [ $$(zcat $$f | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $$f | wc -l) -N SIM.$$(basename $$f) -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $$f; done

SIM-EVAL-JOB: $(foreach g, $(shell cat result/simulation/selected.genes), jobs/eval-$(g)-jobs.txt.gz)
	echo $@

SIM-RANK-JOB: $(foreach g, $(shell cat result/simulation/selected.genes), jobs/rank-$(g)-jobs.txt.gz)
	echo $@

SIM-RANK-POST-JOB: $(foreach g, $(shell cat result/simulation/selected.genes), jobs/rcollect-$(g)-jobs.txt.gz)
	echo $@

# generate jobs for each gene
jobs/data-%-jobs.txt.gz: $(foreach r, $(Rank), \
                         $(foreach h, $(H2), \
                         $(foreach c, $(Ncausal), \
                         $(foreach t, $(Ntissue), \
                         $(foreach s, $(Rseed), jobs/data-%-${r}-$(h)-$(c)-$(t)-$(s)-sub-jobs)))))
	@[ -d log/ ] || mkdir -p log/
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N GTEx.DATA.$* \
           -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/sim-%-jobs.txt.gz: $(foreach r, $(Rank), \
                         $(foreach h, $(H2), \
                         $(foreach c, $(Ncausal), \
                         $(foreach t, $(Ntissue), \
                         $(foreach s, $(Rseed), jobs/sim-%-${r}-$(h)-$(c)-$(t)-$(s)-sub-jobs)))))
	@[ -d log/ ] || mkdir -p log/
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N GTEx.SIM.$* \
           -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/eval-%-jobs.txt.gz: $(foreach r, $(Rank), \
                         $(foreach h, $(H2), \
                         $(foreach c, $(Ncausal), \
                         $(foreach t, $(Ntissue), \
                         $(foreach s, $(Rseed), jobs/eval-%-${r}-$(h)-$(c)-$(t)-$(s)-sub-jobs)))))
	@[ -d log/ ] || mkdir -p log/
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N Eval.$* \
           -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/rank-%-jobs.txt.gz: $(foreach r, $(Rank), \
                         $(foreach h, $(H2), \
                         $(foreach c, $(Ncausal), \
                         $(foreach t, $(Ntissue), \
                         $(foreach s, $(Rseed), jobs/rank-%-${r}-$(h)-$(c)-$(t)-$(s)-sub-jobs)))))
	@[ -d log/ ] || mkdir -p log/
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N GTEx.rank.$* \
           -binding "linear:1" -q short -l h_vmem=2g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@

jobs/rcollect-%-jobs.txt.gz: $(foreach r, $(Rank), \
                         $(foreach h, $(H2), \
                         $(foreach c, $(Ncausal), \
                         $(foreach t, $(Ntissue), \
                         $(foreach s, $(Rseed), jobs/rcollect-%-${r}-$(h)-$(c)-$(t)-$(s)-sub-jobs)))))
	@[ -d log/ ] || mkdir -p log/
	@cat $^ | gzip > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N GTEx.rank.$* \
           -binding "linear:1" -q short -l h_vmem=1g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@


# % = $(gene)-$(rank)-$(h2)-$(nc)-$(nt)-$(rs)
jobs/data-%-sub-jobs:
	@printf "" > $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p $(simdata)/$(shell echo $* | sed 's/-/\//g')
	@[ -f $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.bed ] || \
          echo ./sim_data.sh $(shell echo $* | sed 's/-/ /g') \
               $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp > $@

# % = $(gene)-$(rank)-$(h2)-$(nc)-$(nt)-$(rs)
jobs/sim-%-sub-jobs:
	@printf "" > $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p $(fqtl)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(sqtl)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(glmnet)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(phenix)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(bslmm)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(metasoft)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(mashr)/$(shell echo $* | sed 's/-/\//g')
	@mkdir -p $(zfqtl)/$(shell echo $* | sed 's/-/\//g')

	@[ -f $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.pve.gz ] || \
           echo ./run.sim.fqtl.R $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.snp.gz \
           $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.pve.gz >> $@

	@[ -f $(sqtl)/$(shell echo $* | sed 's/-/\//g')-sqtl.pve.gz ] || \
           echo ./run.sim.sqtl.R $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(sqtl)/$(shell echo $* | sed 's/-/\//g')-sqtl.snp.gz \
           $(sqtl)/$(shell echo $* | sed 's/-/\//g')-sqtl.pve.gz >> $@

	@[ -f $(metasoft)/$(shell echo $* | sed 's/-/\//g')-metasoft.tis.gz ] || \
           echo ./run_metasoft.sh $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(metasoft)/$(shell echo $* | sed 's/-/\//g')-metasoft.snp.gz \
           $(metasoft)/$(shell echo $* | sed 's/-/\//g')-metasoft.tis.gz >> $@

	@[ -f $(mashr)/$(shell echo $* | sed 's/-/\//g')-mashr.tis.gz ] || \
           echo ./run_mashr.sh $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(mashr)/$(shell echo $* | sed 's/-/\//g')-mashr.snp.gz \
           $(mashr)/$(shell echo $* | sed 's/-/\//g')-mashr.tis.gz >> $@

	@[ -f $(zfqtl)/$(shell echo $* | sed 's/-/\//g')-zfqtl.tis.gz ] || \
           echo ./run_sim_zfqtl.sh $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(zfqtl)/$(shell echo $* | sed 's/-/\//g')-zfqtl.snp.gz \
           $(zfqtl)/$(shell echo $* | sed 's/-/\//g')-zfqtl.pve.gz >> $@

	@[ -f $(bslmm)/$(shell echo $* | sed 's/-/\//g')-bslmm.hyper.gz ] || \
           echo ./run_gemma.sh $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(bslmm)/$(shell echo $* | sed 's/-/\//g')-bslmm.snp.gz \
           $(bslmm)/$(shell echo $* | sed 's/-/\//g')-bslmm.hyper.gz >> $@

	@[ -f $(glmnet)/$(shell echo $* | sed 's/-/\//g')-glmnet.pve.gz ] || \
           echo ./run.glmnet.R $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(glmnet)/$(shell echo $* | sed 's/-/\//g')-glmnet.snp.gz \
           $(glmnet)/$(shell echo $* | sed 's/-/\//g')-glmnet.pve.gz >> $@

	@[ -f $(phenix)/$(shell echo $* | sed 's/-/\//g')-phenix.pve.gz ] || \
           echo ./run.sim.phenix.R $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(phenix)/$(shell echo $* | sed 's/-/\//g')-phenix.snp.gz \
           $(phenix)/$(shell echo $* | sed 's/-/\//g')-phenix.pve.gz >> $@

# % = $(gene)-$(rank)-$(h2)-$(nc)-$(nt)-$(rs)
jobs/eval-%-sub-jobs:
	@printf "" > $@
	@[ -d result/simulation/power ] || mkdir -p result/simulation/power
	if [ -f $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.pve.gz ] && [ -f $(sqtl)/$(shell echo $* | sed 's/-/\//g')-sqtl.pve.gz ] && [ -f $(metasoft)/$(shell echo $* | sed 's/-/\//g')-metasoft.tis.gz ] && [ -f $(mashr)/$(shell echo $* | sed 's/-/\//g')-mashr.tis.gz ] && [ -f $(zfqtl)/$(shell echo $* | sed 's/-/\//g')-zfqtl.pve.gz ] && [ -f $(bslmm)/$(shell echo $* | sed 's/-/\//g')-bslmm.hyper.gz ] && [ -f $(glmnet)/$(shell echo $* | sed 's/-/\//g')-glmnet.pve.gz ] && [ -f $(phenix)/$(shell echo $* | sed 's/-/\//g')-phenix.pve.gz ] ; then [ -f result/simulation/power/$*-power.gz ] || echo ./run.sim.evaluation.R $(shell echo $* | sed 's/-/\//g') result/simulation/power/$*-power.gz > $@; fi

jobs/rank-%-sub-jobs:
	@printf "" > $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p $(fqtl)/$(shell echo $* | sed 's/-/\//g')

	@[ -f $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.lodds.gz ] || \
           echo ./run.sim.fqtl-rank.R $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp \
           $(simdata)/$(shell echo $* | sed 's/-/\//g')-temp.y.txt.gz \
           $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.lodds.gz >> $@

jobs/rcollect-%-sub-jobs:
	@printf "" > $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d result/simulation/rank ] || mkdir -p result/simulation/rank
	@[ -f $(fqtl)/$(shell echo $* | sed 's/-/\//g')-fqtl.lodds.gz ] && \
	[ -f result/simulation/rank/$*-rank.gz ] || \
	echo ./run.sim.eval.fqtl-rank.R $(shell echo $* | sed 's/-/\//g') \
	     result/simulation/rank/$*-rank.gz > $@

result/simulation/power.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	for f in result/simulation/power/*-power.gz; do cat $$f >> $@; done

result/simulation/fqtl-power.txt.gz: result/simulation/power.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat result/simulation/power.txt.gz | awk '/fqtl/' | gzip > $@

result/simulation/rank.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	for f in result/simulation/rank/*-rank.gz; do cat $$f >> $@; done
