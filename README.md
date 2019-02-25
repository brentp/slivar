# slivar: filter/annotate variants in VCF/BCF format with simple expressions

slivar has sub-commands:
+ [expr](#expr): trio and group expressions and filtering
+ [gnotate](#gnotate): rapidly annotate a VCF/BCF with gnomad
+ [make-gnotate](#gnotate): make a compressed zip file of annotations for use by slivar
+ filter: filter a VCF with an expression

# Table of Contents

* [Commands](#commands)
 * [expr](#expr)
    * [trio](#trio)
    * [Groups](#groups)
 * [Gnotate](#gnotate)
* [Installation](#installation)
* [Attributes](#attributes)
* [How it works](#how-it-works)


## Commands

### expr

`expr` allows filtering on (abstracted) trios and groups. For example, given a VCF (and ped/fam file) with
100 trios, `slivar` will apply an expression with `kid`, `mom`, `dad` identifiers to each trio that it automatically
extracts.

#### trio

when --trio is used, `slivar` finds all trios in a VCF, PED pair and let's the user specify an expression with indentifiers
of `kid`, `mom`, `dad` that is applied to each possible trio. For example, a simple expression to call
*de novo* variants:

```javascript
variant.FILTER == 'PASS' && \                         # 
variant.call_rate > 0.95 && \                         # genotype must be known for most of cohort.
INFO.gnomad_af < 0.001 && \                           # rare in gnomad (must be in INFO [but see below])
kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && \  # alts are 0:hom_ref, 1:het, 2:hom_alt, -1:unknown
kid.DP > 7 && mom.DP > 7 && dad.DP > 7 && \           # sufficient depth in all
(mom.AD[1] + dad.AD[1]) == 0                          # no evidence for alternate in the parents
```

This requires passing variants that are rare in gnomad that have the expected genotypes and do
not have any alternate evidence in the parents. If there are 200 trios in the `ped::vcf` given, then this expression
will be tested on each of those 200 trios.

The expressions are javascript so the user can make these as complex as needed.


```
slivar expr \
   --pass-only \ # output only variants that pass one of the filters (default is to output all variants)
   --vcf $vcf \
   --ped $ped \
   --gnotate $gnomad_af.zip \ # a compressed gnomad zip that allows fast annotation so that `gnomad_af` is available in the expressions below.
   --load functions.js \ # any valid javascript is allowed here.
   --out-vcf annotated.bcf \
   --info "variant.call_rate > 0.9" \ # this filter is applied before the trio filters and can speed evaluation if it is stringent.
   --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 \
                   && kid.AB > 0.25 && kid.AB < 0.75 \
                   && (mom.AD[1] + dad.AD[1]) == 0 \
                   && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 \
                   && kid.DP >= 12 && mom.DP >= 12 && dad.DP >= 12" \
   --trio "informative:kid.GQ > 20 && dad.GQ > 20 && mom.GQ > 20 && kid.alts == 1 && ((mom.alts == 1 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 1))" \
   --trio "recessive:recessive_func(kid, mom, dad)"

```

Note that `slivar` does not give direct access to the genotypes, instead exposing `alts` where 0 is homozygous reference, 1 is heterozygous, 2 is
homozygous alternate and -1 when the genotype is unknown. It is recommended to **decompose** a VCF before sending to `slivar`

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
```
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
```
normal.alts == 0 && normal.DP > 10 \
  && tumor1.AB > 0 \
  && tumor1.AB < tumor2.AB \
  && tumor2.AB < tumor3.AB \
  && tumor3.AB < tumor4.AB
```

to find a somatic variant that has increasing frequency (AB is allele balance) along the tumor time-points.


### Gnotate

More extensive documentation and justification for `gnotate` are [here](https://github.com/brentp/slivar/wiki/gnotate)

This uses a compressed, reduced representation of a single value pulled from a (population VCF) along with a boolean that indicates a
non-pass filter. This can, for example, reduce the 600+ GB of data for the **whole genome and exome** from gnomad to a 1.5GB file
distributed [here](https://s3.amazonaws.com/gemini-annotations/gnomad-2.1.zip).
The zip file encodes the popmax_AF (whichever is higher between whole genome and exome) and the presence of FILTER for every variant
in gnomad.  It can annotate at faster than 10K variants per second.

slivar gnotate --vcf $input_vcf -o $output_bcf --threads 3 -g encoded.zip

## Installation

get the latest binary from: https://github.com/brentp/slivar/releases/latest
This will require libhts.so (from [htslib](https://htslib.org)) to be in the usual places or in a directory indicated in `LD_LIBRARY_PATH`.

or use via docker from: [brentp/slivar:latest](https://hub.docker.com/r/brentp/slivar)

## Attributes

 + anything in the INFO is available as e.g. INFO.CSQ
 + if FORMAT.AB is not present, it is added so one can filter with kid.AB > 0.25 && kid.AB < 0.75
 + variant attributes are: `CHROM`, `POS`, `start`, `end`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`
 + calculated variant attributes include: `aaf`, `hwe_score`, `call_rate`, `num_hom_ref`, `num_het`, `num_hom_alt`, `num_unknown`

 + sample attributes (via `kid`, `mom`, `dad`) include in the FORMAT. available as e.g. kid.AD[1]
 + sample attributes from the ped for `affected`, `sex`, `id` are available as, e.g. kid.sex.
 + sample relations are available as `mom`, `dad`, `kids`. `mom` and `dad` will be undefined if not available and kids will be an empty array.

## How it works

 `slivar` embeds the [duktape javascript engine](https://duktape.org/) to allow the user to specify expressions.
 For each variant, each trio (and each sample), it fills the appropriate `attributes`. This can be intensive for
 VCFs with many samples, but this is done **as efficiently as possible** such that `slivar` can evaluate 10's of
 thousand of variants per second even with dozens of trios.
