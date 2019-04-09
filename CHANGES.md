v0.1.0 (dev)
============
+ allow accessing INFO via `variant.INFO` as well as `INFO`
+ better error message on missing gnotate file
+ add --sample-expr to `slivar expr` sub-command to allow applying an expression to each sample:
   `--sample-expr "hi_quality:sample.DP > 10 && sample.GQ > 20"
+ expose samples in javascript in `$S` object, e.g. `$S["sampleABC"].DP > 10 && $S["sampleXYZ"].DP > 10`

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
