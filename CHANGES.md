v0.1.5 (dev)
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
