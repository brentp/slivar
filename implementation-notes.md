# javascript

The javascript stuff uses the `duktape` javascript engine. This is quite fast, but the
syntax is somewhat onerous. There is a (zero-cost) wrapper in `duko.nim` to simplify this.
That wrapper is used in slivar to make it simpler to set/get values to the javascript objects.

Currently all values in the INFO and FORMAT fields are sent to the javascript objects. There
could be some simple regular expressions to check which fields are used an only fill those
if it improved performance substantially.

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
It might be useful to look at popmax_AC and make sure that value is large enough. This could
get hairy. Ideally, instead of isolating this in utils/vkgnomad.nim, it would be nice to be
able to specify this behavior via the javascript expressions. For now, probably a few simple
exceptions suffice.

# compression

I have to use the nim/zip/zipfiles for decompression and https://github.com/brentp/nim-minizip for
compression (in utils/vkgnomad) because I get problems using `minizip` for compression and problems
with `zlib` for decompression. The compilation flags include `-d:useLibzipSrc` so there is no
dependency on libz.so.
