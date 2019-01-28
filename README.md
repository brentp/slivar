# slivar: filter/annotate variants in VCF/BCF format with simple expressions

slivar finds all trios in a VCF, PED pair and let's the user specify an expression with indentifiers
of `kid`, `mom`, `dad` that is applied to each possible trio. For example, a simple expression to call
*de novo* variants:

```javascript
variant.FILTER == 'PASS' && \                         # 
variant.call_rate > 0.95 && \                         # genotype must be known for most of cohort.
INFO.gnomad_af < 0.001 && \                           # rare in gnomad (must be in INFO)
kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && \  # alts are 0:hom_ref, 1:het, 2:hom_alt, -1:unknown
kid.DP > 7 && mom.DP > 7 && dad.DP > 7 && \           # sufficient depth in all
(mom.AD[1] + dad.AD[1]) == 0                          # no evidence for alternate in the parents
```

This requires passing variants that are rare in gnomad that have the expected genotypes and do
not have any alternate evidence in the parents. If there are 200 trios in the `ped::vcf` given, then this expression
will be tested on each of those 200 trios.

The expressions are javascript so the user can make these as complex as needed.


```
slivar \
   --pass-only \ # output only variants that pass one of the filters (default is to output all variants)
   --vcf $vcf \
   --ped $ped \
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


## Installation

get the latest binary from: https://github.com/brentp/slivar/releases/latest

or use via docker from: [brentp/slivar:latest](https://cloud.docker.com/repository/docker/brentp/slivar)

### Attributes

 + anything in the INFO is available as e.g. INFO.CSQ
 + if FORMAT.AB is not present, it is added so one can filter with kid.AB > 0.25 && kid.AB < 0.75
 + variant attributes are: `CHROM`, `POS`, `start`, `end`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`
 + calculated variant attributes include: `aaf`, `hwe_score`, `call_rate`, `num_hom_ref`, `num_het`, `num_hom_alt`, `num_unknown`

 + sample attributes (via `kid`, `mom`, `dad`) include in the FORMAT. available as e.g. kid.AD[1]
 + sample attributes from the ped for `affected`, `sex` are available as, e.g. kid.sex.


## How it works

 `slivar` embeds the [duktape javascript engine](https://duktape.org/) to allow the user to specify expressions.
 For each variant, each trio (and each sample), it fills the appropriate `attributes`. This can be intensive for
 VCFs with many samples, but this is done **as efficiently as possible** such that `slivar` can evaluate 10's of
 thousand of variants per second even with dozens of trios.
