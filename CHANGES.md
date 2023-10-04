v0.3.1
======
+ add CSQ parsing js code: https://github.com/brentp/slivar/wiki/CSQ
+ expose sample.GT
+ update impact order to include new effects from VEP (#157)

v0.3.0
======
+ bugfix for 0.2.8 regression where SNPs present in gnotate file would be instead annotated as missing #149.
  This happened when annotating a chromosome *after* a chromosome present in the query set that was not in the gnotate file

v0.2.9
======
+ bugfix for 0.2.8 regression where indels present in gnotate file would be instead annotated as missing #149

v0.2.8
======
+ [internal] use same zip library for make-gnotate as gnotate. this should
  improve speed for reading many small chromosomes.
+ don't quit on CSQ/ANN/BCSQ fields that don't have enough information to parse (#122)
+ fix segregating_dominant_x

v0.2.7
======
+ [tsv] fix for empty VCFs

v0.2.6
======
+ allow specifying `QUAL` and `ID` to `-i` argument to `slivar tsv` (#104 #105)

v0.2.5
======
+ it was previously not possible to adjust the order of impacts with
  `SLIVAR_IMPACTFUL_ORDER`. This is now fixed. (#97)

v0.2.4
======
+ fix long-standing bug (#27) that would cause sigsegv in some cases

v0.2.3
======
+ add `feature_fusion` to impact order list (#92)
+ [compound-hets] write summary even if no variants found

v0.2.2
======
+ fix bug with '.' in ALT field (caused message about incorrect number of alts in some cases. thanks Batsal for reporting)
+ slivar compound-hets: fix bug with parents specified in ped file but absent from VCF (#79)
+ slivar compound-hets: add "intergenic_region" to list of impacts that are skipped by default
+ add SLIVAR_NO_REPORT_ALL to prevent reporting variants for familys with no affected samples
+ fix bug when creating huge zip files with make-gnotate (#86)
+ **NOTE**: change default min depth in slivar-functions.js to 6. (was 0)
+ don't hard-code tmp directory to /tmp (use $TMPDIR)
+ [compound-hets] dont fail if no usable variants were found, just issue warning.
  this can happen for small chroms or regions when slivar was parallelized.

v0.2.1
======
+ add `INFO.highest_impact_order` so user can do things like:
  ```
  --info "INFO.highest_impact_order == ImpactOrder.missense"
  ```
  to get only missense mutations. or

  ```
  --info "INFO.highest_impact_order <= ImpactOrder.missense"
  ```
  to get variants with an impact of missense or higher. (note that lower values have higher impact)
+ `tsv`: bugfix!! previously, `slivar tsv` used a `< 0` check where it should hvae used `<= 0` so
   candidates with a parental sample at index 0 in the samples of the vcf would have missing information
   for depth, GQ, and AB. (#78 . Thanks @amwenger for finding the problem and its source).
+ `expr`: always annotate variant and INFO even if no sample (or trio,group,fam) expressions are given so one can use slivar just to get the 'impactful' annotation.
+ remove duplicates from default-order.txt

v0.2.0
======
+ clear values from expressions for family, group, and sample. was previously only done for
  trios. so if you relied on something like:
  ```
  --group-expr "expr1: ..." \
  --group-expr "expr2:('expr1' in INFO) && ..."\
  ```
  then the expr2 would still see expr1 in INFO even if the preceding expression did not pass but it
  has passed for a previous variant. Thanks Amelia Wallace for finding and reporting a simple test-case.
+ add /opt/slivar/slivar-functions.js to docker image (#75)


v0.1.13
=======
+ fix bug with VCFs with >256 info fields (#74). Thanks @liserjrqlxue for making a great test-case.

v0.1.12
=======
+ slivar expr: handle ploidy > 2 (updated hts-nim) (#55 thanks to @markw3lsh for reporting and providing a test-case)
+ slivar: optionally support strings in format (sample) fields. default is to ignore these for performance reasons, but users can set the environment variable `SLIVAR_FORMAT_STRINGS` to anything to force `slivar` to populate the javascript objects with any strings in the VCF. e.g. `SLIVAR_FORMAT_STRINGS=1 slivar expr ...` and use `unset SLIVAR_FORMAT_STRINGS` to remove. (#37)
+ slivar compound-hets: **fix bug with 2 affect sibs sharing same compound-het previously only 1 sib would be reported**

v0.1.11
=======
+ gnotate: ~10-20% speed improvement by inlining cmp in binarysearch
+ slivar expr: fix counts in summary table when some families in the cohort have no affected samples
+ slivar tsv: improve readability when 1 or both parents are missing
+ slivar tsv: fix comma separation with 1 or both parents missing (thanks @brwnj for reporting)
+ slivar ddc: when looking at all chromosomes, skip X, Y
+ slivar expr: add --exclude option which takes a bed file of exclude regions
+ slivar compound-hets: default will now output compound-het pairs that include synonymous variants 
+ pslivar: fix tsv output from pslivar (which allows running slivar in parallel). more info [here](https://github.com/brentp/slivar/wiki/parallel-slivar)
+ add `INFO.genic` boolean to complement `INFO.impactful`. By default this includes all impacts included by `impactful` along with
  synonymous and other exonic, but non-protein-altering variants (but does not include `intronic`. (more info [here](https://github.com/brentp/slivar/wiki/impactful)

v0.1.10
=======
+ [slivar ddc] fix bug when selecting sample filters.
+ [slivar ddc] UI fix: remove button when [x] is clicked.


v0.1.9
======
+ new command `slivar ddc`: https://github.com/brentp/slivar/wiki/data-driven-cutoffs
+ fix for slivar compound-hets reporting orphaned variants in some cases (thanks Steve Boyden for reporting)
+ update to htslib >= 1.10
+ change recommended setup for rare disease to use family-based approach, rather than trios as it is more
  flexible and now well-tested. (https://github.com/brentp/slivar/wiki/rare-disease#full-analysis-for-trios-with-unaffected-parents)
+ add `VCF` object available in javascript which gets populated with `VCF.CSQ == ["CONSEQUENCE", "CODONS","AMINO_ACIDS", "GENE", ...]` and VCF.ANN, VCF.BCSQ if any
  or all of those annotations are available in the VCF (see #3)

v0.1.8
======
+ fix bug when subsetting vcfs (ped file contained subset of samples in vcf). Thanks Matt for reporting.

v0.1.7
======
+ slivar make-gnotate will error with warning if field is not a float or int
+ compound-het provides a way to ignore some impacts (intron, non_coding, etc) and sets sane defaults. this removes
  most candidates from WGS which were predominantly pairs of intronic or non-coding variants.
+ compound-het: remove -f flag, this is now discovered by default
+ only samples with parents in the vcf (not just ped) are counted as trios
+ expose sample.phenotype attribute which is a string taken directly from the pedigree file
+ given an error on an INFO expression, don't bail, just report and continue
+ fix for VCFs without index AND without contigs in header (#44)
+ support for file sizes > 4.2GB and support for gnotate files with chr prefix annotating query files without
+ fix gnotate bug where some long alleles were not annotated


v0.1.6
======
+ add INFO.impactful (boolean) if CSQ/BCSQ/ANN is present. this value will be true of any of those annotation are high enough
+ the `impactful` flag is automaticaly added to the output VCF for any variant passing through slivar and meeting the criteria (https://github.com/brentp/slivar/wiki/impactful)
+ better checks for length of AD field
+ fix for empty groups (#38)
+ add --family-expr for family expressions like: 
    `fam.every(function(s) { return s.het == s.affected && s.hom_ref == !s.affected && s.GQ > 5 })`
+ slivar tsv now outputs a sortable column for highest-impact. it makes use of a default
  list of impact orderings from: https://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html and
  supplemented with any type seen in bcftools or snpEff.
+ bug fix with gnotate (-g) when annotating a file with a empty chromsome followed by a non-empty (e.g. if `chr16_random` is mixed in with canonical chromosomes)

v0.1.5
======
+ !!MAJOR: fix bug with `hom_alt` alias for `alts == 2`. Since v0.1.3, `hom_alt` would not get set to false for samples, variants after
  being set once.
+ update [wiki](https://github.com/brentp/slivar/wiki/rare-disease#full-analysis-for-trios-with-unaffected-parents) to simplify rare disease analysis.
+ tune javascript in [my.js](https://raw.githubusercontent.com/brentp/slivar/master/js/my.js) and rare-disease expressions so that many analyses will run
  better than twice as fast.
+ [expr] add `--skip-non-variable` which further improves speed by not evaluating javascript for a trio if all members are hom-ref or unknown (but can the same variant will
  still be evlauated for trios with at least 1 sample that has a variant at that site.
+ [compound-het] major speed increase for large cohorts. 
+ [compound-het] remove -i/--index option as this is now discovered from the header for the CSQ/BCSQ field.

v0.1.4
======
+ better error messages on bad VCF
+ [tsv] add --csq-column to `slivar tsv` to allow extracting extra CSQ fields
+ general usability improvements
+ [compound-het]: support singleton (kid only) or duo's (kid and 1 parent) to compound het. this will give many false positives because
  the variants can't be phased, but the number can still be quite small given sufficient filtering on population allele frequency.
  requires --allow-non-trios argument

v0.1.3
======
+ add new sample attributes `hom_ref`, `het`, `hom_alt`, `unknown` which are synonyms of `alts == 0`, `alts == 1`, `alts == 2`, `alts == -1` respectively.
+ remove `slivar gnotate` sub-command. the same functionality can be had from `slivar expr` with `-g` and `--info`.

v0.1.2
======
+ new sub-command `slivar tsv` to output a tab-separated value file from a filtered VCF for final processing
  this command also allows adding columns based on the gene the variant is in, such as description of gene
  function, pLI, etc.
+ fix bug where slivar expr would segfault if labels were re-used (now quits with error message)
+ better help in `slivar make-gnotate`
+ improve and document output format for `duo-del` https://github.com/brentp/slivar/wiki/finding-deletions-in-parent-child-duos
+ add a prelude function with hasSample(INFO, key, sample) to see if sample has been added to info field by previous filter.
+ allow outputting summary table to file with `SLIVAR_SUMMARY_FILE` environment variable
+ fix for make-gnotate with multiple files (still requires first file to contain all chromosomes)
+ allow more kinds of pedigree files (#30)
+ change format of slivar compound het field: adds an id that uniquely identifies the pair of variants.

v0.1.1
======
+ better checks on missing parents/kids
+ new sub-command `duo-del` that uses non-transmission in parent-child duos to find *de novo* structural deletions.
  e.g. a cluster of sites where kid is 0/0 and parent is 1/1 would be a candidate.
  see: https://github.com/brentp/slivar/wiki/finding-de-novo-deletions-in-parent-child-duos

v0.1.0
======
+ allow accessing INFO via `variant.INFO` as well as `INFO`
+ better error message on missing gnotate file
+ add --sample-expr to `slivar expr` sub-command to allow applying an expression to each sample:
   `--sample-expr "hi_quality:sample.DP > 10 && sample.GQ > 20"
+ expose samples in javascript in `$S` object, e.g. `$S["sampleABC"].DP > 10 && $S["sampleXYZ"].DP > 10`
+ better checking of matching between ped and vcf => don't fail on missing kids

v0.0.9
======
+ fix bug without --pass-only. slivar always behaved as if --pass-only was used.
+ [compound-hets] fix for multi-sample VCFs with `--sample-field`

v0.0.8
======
+ use random file names so concurrent slivar processes don't clobber files.
+ more informative error on bad js expression
+ fix for empty groups (#20)
+ fix bug when later expressions depended on previous ones.


v0.0.7
======
+ [expr] allow --region to be a bed file
+ [compound-hets] new command for finding compound-hets including where 1 side is a de novo (see: https://github.com/brentp/slivar/wiki/rare-disease#compound-heterozygotes)
+ default to send output to stdout for streaming
+ fix bug in some -g gnotate files
