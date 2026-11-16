# javascript

The javascript engine is JavaScriptCore (JSC), through https://github.com/brentp/jsc-nim. On linux
that downloads a prebuilt, statically-linkable JSC from Bun's WebKit fork at compile time.

Calls across the JSC C API are relatively expensive (~150-300ns each), while JIT-compiled javascript
is very fast, so the design minimizes crossings:

+ nim decodes INFO and FORMAT values into C buffers that javascript sees as typed arrays
  (`Buf` in evaluator.nim). Nothing is copied into javascript objects.
+ `INFO`, `variant` and each sample (`kid`, `mom`, `sample`, ...) are javascript `Proxy` objects
  (see `runtime.js`) that read from those typed arrays when an attribute is accessed. They
  also implement the "strict" behavior: accessing a missing attribute throws `unknown attribute:x`
  unless `SLIVAR_NO_ATTRIBUTE_CHECKING` is set.
+ strings (CHROM, REF, INFO strings, ...) are written to a shared utf-8 arena and decoded lazily.
+ per variant, nim makes at most two calls: `$sl.info()` for the `--info` expression and
  `$sl.samples()` which loops over all trios/groups/families and expressions in javascript
  and writes the indexes of matches into a shared `Int32Array`.

# gnotate

the encoding used by `gnotate` can handle ref+alt length of up to 14 bases. For variants
longer than that, it still adds an entry to the compressed version, but with an empty ref
and alt allele. The `find` function will see there is an entry in there with missing ref, alt
and it knows to then look into the `longs` datastructure (which contains < 0.1% of the data)
to check for large alleles.

Both the `longs` and the `encs` are sorted and use binary search. The `enc` is a simple
encoding to uint64 such that ordering is maintained. This setup allows 5-7 million searches
per second and should not be a bottleneck as such.

The data is stored in a .zip file with 1 directory per chromosome and a list of FILTERs since
`slivar` annotates with the gnomad filters as well.

It defaults to using `popmax_AF` but this can sometimes be empty even when no FILTER is set.
It might be useful to look at `popmax_AC` and make sure that value is large enough. This could
get hairy. Ideally, instead of isolating this in utils/vkgnomad.nim, it would be nice to be
able to specify this behavior via the javascript expressions. For now, probably a few simple
exceptions suffice.

# filling sample and info fields

slivar decodes every FORMAT field for every sample. If variant `i` has `GT:GQ:DP:AD:PL` but
variant `i+1` has only `GT:GQ:AD`, then `DP` and `PL` will be unavailable. Each header id has a
count in `FMT_N` / `INFO_N` (0 means absent for this variant); nim resets only the ids that were set
for the previous variant.

# compression

I have to use the nim/zip/zipfiles for decompression and https://github.com/brentp/nim-minizip for
compression (in utils/vkgnomad) because I get problems using `minizip` for compression and problems
with `zlib` for decompression. The compilation flags include `-d:useLibzipSrc` so there is no
dependency on libz.so.
