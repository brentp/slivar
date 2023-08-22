# slivar: filter/annotate variants in VCF/BCF format with simple expressions [![Build Status](https://travis-ci.com/brentp/slivar.svg?branch=master)](https://travis-ci.com/brentp/slivar)
[![downloads](https://anaconda.org/bioconda/slivar/badges/downloads.svg)](https://anaconda.org/bioconda/slivar)


If you use `slivar`, please cite [the paper](https://www.nature.com/articles/s41525-021-00227-3)

slivar is a set of command-line tools that enables rapid querying and filtering of VCF files. 
It facilitates operations on trios and [groups](#groups) and allows arbitrary expressions using simple javascript.

#### use-cases for `slivar`

+ annotate variants with [gnomad](https://gnomad.broadinstitute.org/) allele frequencies from combined exomes + whole genomes at > 30K variants/second using only a 1.5GB compressed annotation file.
+ call *denovo* variants with a simple expression that uses *mom*, *dad*, *kid* labels that is applied to each trio in a cohort (as inferred from a pedigree file).
  `kid.het && mom.hom_ref && dad.hom_ref && kid.DP > 10 && mom.DP > 10 && dad.DP > 10`
+ define and filter on arbitrary groups with labels. For example, 7 sets of samples each with 1 normal and 3 tumor time-points:
  `normal.AD[0] = 0 && tumor1.AB  < tumor2.AB && tumor2.AB < tumor3.AB`
+ filter variants with simple expressions:
  `variant.call_rate > 0.9 && variant.FILTER == "PASS" && INFO.AC < 22 && variant.num_hom_alt == 0`
+ see [using slivar for rare disease research](https://github.com/brentp/slivar/wiki/rare-disease)

![slivar logo](https://user-images.githubusercontent.com/1739/85797163-55738480-b6f8-11ea-8486-389de3e492e6.png)

slivar has sub-commands:
+ [expr](#expr): filter and/or annotate with INFO, trio, sample, group expressions
+ [make-gnotate](#make-gnotate): make a compressed zip file of annotations for use by slivar
+ [compound-hets](#compound-het): true compound hets using phase-by-inheritance within gene annotations

# Table of Contents

* [Installation](#installation)
* [QuickStart](#quickstart)
* [Commands](#commands)
 * [expr](#expr)
    * [trio](#trio)
    * [Family Expressions](#family-expressions)
    * [Groups](#groups)
 * [compound-het](#compound-het)
 * [tsv](#tsv)
 * [duo-del](#duo-del)
 * [ddc](#data-driven-cutoffs)
* [Attributes](#attributes)
* [How it works](#how-it-works)
* [Gnotation Files](#gnotation-files)



## Installation

get the latest binary from: https://github.com/brentp/slivar/releases/latest

`slivar_static` does not depend on any libraries and should work on any 64 bit linux system.        

`slivar_shared` will require libhts.so (from [htslib](https://htslib.org)) to be in the usual places or in a directory indicated in `LD_LIBRARY_PATH`.

or use via docker from: [brentp/slivar:latest](https://hub.docker.com/r/brentp/slivar)

## QuickStart

To get started quickly, grab a static binary [for the latest release](https://github.com/brentp/slivar/releases/latest) and
then follow [this example](https://github.com/brentp/slivar/wiki/rare-disease#full-analysis-for-trios-with-unaffected-parents)

So for hg38:

```
vcf=/path/to/your/vcf.vcf.gz
ped=/path/to/your/pedigree.ped
wget https://github.com/brentp/slivar/releases/download/v0.2.8/slivar
chmod +x ./slivar
wget https://raw.githubusercontent.com/brentp/slivar/master/js/slivar-functions.js
wget https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip

# example command
./slivar expr --js slivar-functions.js -g gnomad.hg38.genomes.v3.fix.zip \
	--vcf $vcf --ped $ped \
	--info "INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS'" \
	--trio "example_denovo:denovo(kid, dad, mom)" \
	--family-expr "denovo:fam.every(segregating_denovo)" \
	--trio "custom:kid.het && mom.het && dad.het && kid.GQ > 20 && mom.GQ > 20 && dad.GQ > 20" \
	--pass-only
```

The pedigree format is explained [here](https://github.com/brentp/slivar/wiki/pedigree-file)

## Commands

### expr

`expr` allows filtering on (abstracted) trios and groups. For example, given a VCF (and ped/fam file) with
100 trios, `slivar` will apply an expression with `kid`, `mom`, `dad` identifiers to each trio that it automatically
extracts.

`expr` can also be used, for example to annotate with population allele frequencies from a `gnotate` file without
any sample filtering. See [the wiki](https://github.com/brentp/slivar/wiki/gnotate) for more detail and [the gnotate](#gnotation-files)
section for gnotation files that we distribute for `slivar`.

`expr` commands are quite fast, but can be parallelized using [pslivar](https://github.com/brentp/slivar/wiki/parallel-slivar).

#### trio

when --trio is used, `slivar` finds all trios in a VCF, PED pair and let's the user specify an expression with indentifiers
of `kid`, `mom`, `dad` that is applied to each possible trio. For example, a simple expression to call
*de novo* variants:

```javascript
variant.FILTER == 'PASS' && \                         # 
variant.call_rate > 0.95 && \                         # genotype must be known for most of cohort.
INFO.gnomad_af < 0.001 && \                           # rare in gnomad (must be in INFO [but see below])
kid.het && mom.hom_ref && dad.hom_ref && \            # also unknown
kid.DP > 7 && mom.DP > 7 && dad.DP > 7 && \           # sufficient depth in all
(mom.AD[1] + dad.AD[1]) == 0                          # no evidence for alternate in the parents
```

This requires passing variants that are rare in gnomad that have the expected genotypes and do
not have any alternate evidence in the parents. If there are 200 trios in the `ped::vcf` given, then this expression
will be tested on each of those 200 trios.

When trios are not sufficient, use [Family Expressions](#family-expressions) which allow more heterogeneous
family structures.

The expressions are javascript so the user can make these as complex as needed.


```bash
slivar expr \
   --pass-only \ # output only variants that pass one of the filters (default is to output all variants)
   --vcf $vcf \
   --ped $ped \
   # compressed zip that allows fast annotation so that `gnomad_af` is available in the expressions below.
   --gnotate $gnomad_af.zip \ 
   # any valid javascript is allowed in a file here. provide functions to be used below.
   --js js/slivar-functions.js \ 
   --out-vcf annotated.bcf \
   # this filter is applied before the trio filters and can speed evaluation if it is stringent.
   --info "variant.call_rate > 0.9" \ 
   --trio "denovo:kid.het && mom.hom_ref && dad.hom_ref \
                   && kid.AB > 0.25 && kid.AB < 0.75 \
                   && (mom.AD[1] + dad.AD[1]) == 0 \
                   && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 \
                   && kid.DP >= 12 && mom.DP >= 12 && dad.DP >= 12" \
   --trio "informative:kid.GQ > 20 && dad.GQ > 20 && mom.GQ > 20 && kid.alts == 1 && \
           ((mom.alts == 1 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 1))" \
   --trio "recessive:trio_autosomal_recessive(kid, mom, dad)"

```

Note that `slivar` does not give direct access to the genotypes, instead exposing 
`hom_ref`, `het`, `hom_alt` and `unknown` or via `alts` where 0 is homozygous reference, 1 is heterozygous, 2 is
homozygous alternate and -1 when the genotype is unknown. It is recommended to **decompose** a VCF before sending to `slivar`

Here it is assumed that `trio_autosomal_recessive` is defined in `slivar-functions.js`; an example implementation of that
and other useful functions is provided [here](https://github.com/brentp/slivar/blob/4856c503a15f2647270a2ac24e4e1b1455208e34/js/slivar-functions.js).
Note that it's often better to use --family-expr instead as it's more flexible than trio expressions.


#### Family Expressions

Trios are a nice abstraction for cohorts consisting of only trios, but for more general uses, there is `--family-expr`
for example, given either a duo, or a quartet, we can find variants present only in affected samples with:

```		
 --family-expr "aff_only:fam.every(function(s) { return s.het == s.affected && s.hom_ref == !s.affected && s.GQ > 5 })"
```

Note that this does not explicitly check for transmission or non-transmission between parents and off-spring
so it is less transparent than the `trio` mode, but more flexible.


#### Groups

A `trio` is a special-case of a `group` that can be inferred from a pedigree. For more specialized use-cases, a `group` can be
specified. For example we could, instead of  using `--trio`, use a `group` file like:
```
#kid	mom	dad
sample1	sample2	sample3
sample4	sample5	sample6
sample7	sample8	sample9
```

Where, here we have specified 3 trios below a header with their "labels". This can be accomplished using `--trio`, but we can
for example specify quartets like this:

```
#kid	mom	dad	sibling
sample1	sample2	sample3	sample10
sample4	sample5	sample6	sample11
sample7	sample8	sample9	sample12
```

where `sample10` will be available as "sibling" in the first family and an expression like:
```bash
kid.alts == 1 && mom.alts == 0 && dad.alts == 0 and sibling.alts == 0
```
could be specified and it would automatically be applied to each of the 3 families.

Another example could be looking at somatic variants with 3 samples, each with a normal and 4 time-points of a tumor:
```
#normal	tumor1	tumor2	tumor3	tumor4
ss1	ss8	ss9	ss10	ss11
ss2	ss12	ss13	ss14	ss15	
ss3	ss16	ss17	ss18	ss19	
```

where, again each row is a sample and the ID's (starting with "ss") will be injected for each sample to allow a single
expression like:
```bash
normal.hom_ref && normal.DP > 10 \
  && tumor1.AB > 0 \
  && tumor1.AB < tumor2.AB \
  && tumor2.AB < tumor3.AB \
  && tumor3.AB < tumor4.AB
```

to find a somatic variant that has increasing frequency (AB is allele balance) along the tumor time-points.
More detail on groups is provided [here](https://github.com/brentp/slivar/wiki/groups-in-slivar)

#### Sample Expressions

Users can specify a boolean expression that is tested against each `sample` using e.g.:

```
--sample-expr "hi_quality:sample.DP && sample.GQ > 10"
```

Each sample that passes this expression will be have its sample id appended to the INFO field of `hi_quality` which
is added to the output VCF.


#### make-gnotate

Users can make their own `gnotate` files like:

```bash
slivar make-gnotate --prefix gnomad \
    --field AF_popmax:gnomad_popmax_af \
    --field nhomalt:gnomad_num_homalt \
    gnomad.exomes.r2.1.sites.vcf.gz gnomad.genomes.r2.1.sites.vcf.gz
```

this will pull `AF_popmax` and `nhomalt` from the INFO field and put them into `gnomad.zip` as `gnomad_popmax_af` and `gnomad_num_homalt` respectively.
The resulting zip file will contain the union of values seen in the exome and genomes files with the maximum value for any intersection.
Note that the names (`gnomad_popmax_af` and `gnomad_num_homalt` in this case) should be chosen carefully as those will be the names added to the INFO of any file to be annotated with the resulting `gnomad.zip`

More information on `make-gnotate` is [in the wiki](https://github.com/brentp/slivar/wiki/make-gnotate)

### compound-het

This command is used to find compound heterozygous variants (with phasing-by-inheritance) in trios.
It is used after filtering to rare(-ish) heterozygotes.

See a full description of use [here](https://github.com/brentp/slivar/wiki/compound-heterozygotes)

**NOTE** that by default, this command limits to a subset of impacts; this is
adjustable with the `--skip` flag. See more on the
[wiki](https://github.com/brentp/slivar/wiki/compound-heterozygotes)

### tsv

This command is used to convert a filtered and annotated VCF to a TSV (tab-separated value file) for final 
examination. An example use is:

```
slivar tsv -p $ped \
    -s denovo -s x_recessive \
    -c CSQ \
    -i gnomad_popmax_af -i gnomad_nhomalt \
    -g gene_desc.txt -g clinvar_gene_desc.txt \
    $vcf > final.tsv
```

where `denovo` and `x_recessive` indicate the INFO fields that contain lists of samples (as added by slivar) that should be extracted.
and `gnomad_popmax_af` and `gnomad_nhomalt` are pulled from the INFO field. 
The `-c` arugment (CSQ) tells `slivar` that it can get gene, transcript and impact information from the CSQ field in the INFO.
And the `-g` arguments are tab-delimited files of gene -> description where the description is added to the text output for quick inspection.
Run `slivar tsv` without any arguments for examples on how to create these for pLI and clinvar.

Also see the [wiki](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF)

## duo-del

slivar duo-del finds structural deletions in parent-child duos using non-transmission of alleles. this 
can work to find deletions in exome data using genotypes, thereby avoiding the problems associated with
depth-based CNV calling in exomes.

see: https://github.com/brentp/slivar/wiki/finding-deletions-in-parent-child-duos

## Data Driven Cutoffs

`slivar ddc` is a tool to discover data-driven cutoffs from a VCF and pedigree information.
It generates an interative VCF so a user can see how **mendelian violation and transmissions**
are effected by varying cutoffs for values in the INFO and FORMAT fields.

See [the wiki](https://github.com/brentp/slivar/wiki/data-driven-cutoffs) for more details.

## Attributes

 + anything in the INFO is available as e.g. INFO.DP
 + INFO.impactful which, if CSQ (VEP), BCSQ (bcftools), or ANN (snpEff) is present indicates if the highest impact is "impactful". see [wiki](https://github.com/brentp/slivar/wiki/impactful) and `INFO.genic` which includes other gene impacts like `synonymous`. Also `INFO.highest_impact_order` explained in the wiki
 + variant consequences such as in INFO.CSQ can be parsed and used as object as described [here](https://github.com/brentp/slivar/wiki/CSQ)
 + if FORMAT.AB is not present, it is added so one can filter with kid.AB > 0.25 && kid.AB < 0.75
 + variant attributes are: `CHROM`, `POS`, `start`, `end`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`,
                           `is_multiallelic`
 + calculated variant attributes include: `aaf`, `hwe_score`, `call_rate`, `num_hom_ref`, `num_het`, `num_hom_alt`, `num_unknown`

 + numeric and flag sample attributes (via `kid`, `mom`, `dad`) included in the FORMAT. available as e.g. `kid.AD[1]`, `mom.DP`, etc.
 + if the environment variable `SLIVAR_FORMAT_STRINGS` is not empty, then string sample fields will be available. these are not populated
   by default as they are used less often and impact performance.
 + sample attributes for `hom_ref`, `het`, `hom_alt`, `unknown` which are synonums for `sample.alts` of 0, 1, 2, -1 respectively.
 + sample attributes from the ped for `affected`, `phenotype`, `sex`, `id` are available as, e.g. kid.sex.
   phenotype is a string taken directly from the pedigree file while affected is a boolean.
 + sample relations are available as `mom`, `dad`, `kids`. `mom` and `dad` will be undefined if not available and kids will be an empty array.
 + a `VCF` object contains `CSQ`, `BCSQ`, `ANN` if those are present in the header (from VEP, BCFTOOLS, SnpEFF). The content is a list indicating
   the order of entries in the field e.g. `["CONSEQUENCE", "CODONS","AMINO_ACIDS", "GENE", ...]`

## How it works

 `slivar` embeds the [duktape javascript engine](https://duktape.org/) to allow the user to specify expressions.
 For each variant, each trio (and each sample), it fills the appropriate `attributes`. This can be intensive for
 VCFs with many samples, but this is done **as efficiently as possible** such that `slivar` can evaluate 10's of
 thousand of variants per second even with dozens of trios.

## Summary Table

slivar outputs a summary table with rows of samples and columns of expression where each value
indicates the number of variants that passed the expression in each sample. By default, this goes to STDOUT
but if the environment variable `SLIVAR_SUMMARY_FILE` is set, `slivar` will write the summary to that file
instead.

## Gnotation Files

Users can create their own gnotation files with `slivar make-gnotate`, but we provide:

+ gnomad for hg37 with AF popmax, numhomalts (total and controls only) [here](https://s3.amazonaws.com/slivar/gnomad.hg37.zip)
+ gnomad for hg38 (v3) genomes [here](https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip)


+ **lifted** gnomad exomes+genomes for hg38 with AF popmax, numhomalts (updated in release v0.1.2) [here](https://s3.amazonaws.com/slivar/gnomad.hg38.v2.zip)
<!--
+ gnomad genomes (71,702 samples) for hg38 with AF popmax, numhomalts (updated in release v0.1.7) [here](https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.zip)
-->
+ spliceai scores (maximum value of the 4 scores in spliceai) [here](https://s3.amazonaws.com/slivar/spliceai.hg37.zip)

+ [topmed allele frequencies (via dbsnp)](https://slivar.s3.amazonaws.com/topmed.hg38.dbsnp.151.zip) these can be used with `INFO.topmed_af`. Useful when analyzing data in hg38 because [some variants in hg38 are not visible in GRCh37](https://twitter.com/brent_p/status/1139540523364917248)

The available fields can be seen with, for example:

```
$ unzip -l gnomad.hg38.v2.zip | grep -oP "gnotate-[^.]+" | sort -u
gnotate-gnomad_nhomalt
gnotate-gnomad_nhomalt_controls
gnotate-gnomad_popmax_af
gnotate-gnomad_popmax_af_controls
gnotate-variant
```

indicating that `INFO.gnomad_nhomalt`, `INFO.gnomad_nhomalt_controls`, `INFO.gnomad_popmax_af` and `INFO.gnomad_popmax_af_controls` will be
the fields after they are added to the INFO.

