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
