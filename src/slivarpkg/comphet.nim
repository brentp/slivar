import hts/vcf
import pedfile
import ./counter
import ./tsv
import strformat
import strutils
import algorithm
import argparse
import tables
import sets

template inherited(kid: pedfile.Sample, a:openarray[int8]): bool =
  # this avoids checks as they are already done in is_compound het
  a[kid.i] == 1 and (a[kid.mom.i] == 1 or a[kid.dad.i] == 1)

template denovo(kid: pedfile.Sample, a:openarray[int8]): bool =
  # this avoids checks as they are already done in is_compound het
  a[kid.i] == 1 and (a[kid.mom.i] == 0 and a[kid.dad.i] == 0)

proc compound_denovo(kid: pedfile.Sample, a:openarray[int8], b:openarray[int8]): bool {.inline.} =
  ## this should only be callled from is_compount_het as it avoids re-checking stuff already in there.
  ## it checks if there is transmitted het and 1
  if kid.denovo(a) and kid.inherited(b): return true
  if kid.denovo(b) and kid.inherited(a): return true
  return false

proc is_compound_het*(kid: pedfile.Sample, a:openarray[int8], b:openarray[int8], allow_non_trios: bool=false): bool {.inline.} =
  # check if the kid (with mom and dad) has either a classic compound het
  # where each parent transmits 1 het to the kid or where 1 het is transmitted
  # and the other is de novo.
  if a[kid.i] != 1 or b[kid.i] != 1: return false
  let dad_nil = kid.dad == nil or kid.dad.i == -1
  let mom_nil = kid.mom == nil or kid.mom.i == -1

  if allow_non_trios:
    if dad_nil and mom_nil:
      return true

    if dad_nil or mom_nil:
      var p = if mom_nil: kid.dad else: kid.mom
      # parent must be het at 1 and only 1 site. assume the missing parent is
      # het at the other site.
      return (a[p.i] == 0 and b[p.i] == 1) or (a[p.i] == 1 and b[p.i] == 0)


  doAssert kid.dad != nil and kid.mom != nil and kid.dad.i != -1 and kid.mom.i != -1, $(kid, kid.dad, kid.mom)

  if a[kid.dad.i] == -1 or a[kid.mom.i] == -1: return false
  if b[kid.dad.i] == -1 or b[kid.mom.i] == -1: return false

  # if only a single het between the parents, only possiblity is a de novo.
  if a[kid.dad.i] + a[kid.mom.i] + b[kid.dad.i] + b[kid.mom.i] == 1:
    return compound_denovo(kid, a, b)

  if a[kid.dad.i] + a[kid.mom.i] != 1: return false

  if b[kid.dad.i] + b[kid.mom.i] != 1: return false
  # now we have 1 het parent at each site and kid het at both sites.
  # check that different parents have the different alleles.
  return a[kid.dad.i] + b[kid.dad.i] == 1 and a[kid.mom.i] + b[kid.mom.i] == 1

proc kids(samples:seq[Sample], allow_non_trios:bool): seq[Sample] =
  for s in samples:
    if allow_non_trios and (s.kids.len == 0):
      # NOTE we skip parents here on purpose.
      result.add(s)
    elif s.dad != nil and s.mom != nil:
      result.add(s)

proc key(v:Variant): string {.inline.} =
  var balts = join(v.ALT, ",")
  &"{$v.CHROM}/{v.start+1}/{v.REF}/{balts}"

proc add_comphet(a:Variant, b:Variant, gene: string, sample:string, id:int) =
  var s: string = ""
  #sample/gene/id/chrom/pos/ref/alt
  discard a.info.get("slivar_comphet", s)
  if s.len > 1:
    # htslib doesn't allow setting an existing field again. so we delete and re-add.
    doAssert a.info.delete("slivar_comphet") == Status.OK
    s &= ","
  # NOTE: if this format changes, must change slivar tsv as well
  s &= &"{sample}/{gene}/{id}/{b.key}"
  doAssert a.info.set("slivar_comphet", s) == Status.OK

proc get_samples(v:Variant, sample_fields: seq[string], samples: var HashSet[string]): bool =
  samples.init(4)
  var tmp: string
  for sf in sample_fields:
      if v.info.get(sf, tmp) != Status.OK: continue
      for s in tmp.split(seps={','}):
        samples.incl(s)
  return sample_fields.len == 0 or samples.len > 0


proc write_compound_hets(ovcf:VCF, kids:seq[Sample], tbl:TableRef[string, seq[Variant]], sample_fields: seq[string], comphet_id: var int, allow_non_trios: bool, cnt: var Counter): int =
  # sample_fields is a seq like @['lenient_ar', 'lenient_denovo'] where each is
  # a key in the INFO field containing a (comma-delimited list of samples)
  var
    x: seq[int32]
    a: seq[int8]
    b: seq[int8]

  var foundKeys = initHashSet[string]()
  # need to track variants in this because 1 variant can be a CH with multiple genes/transcripts.
  var found: seq[Variant]

  for gene, variants in tbl.mpairs:
    var variants = variants
    if variants.len < 2: continue

    # cache so we don't re-get the same alts repeatedly
    # especially important because getting the b alts is inside
    # of the kids loop.
    var cache = newSeq[seq[int8]](variants.len)
    cache[0] = variants[0].format.genotypes(x).alts

    # for each variant, just grab a hashset of samples 1 time.
    var sample_cache = newSeq[HashSet[string]](variants.len)
    discard variants[0].get_samples(sample_fields, sample_cache[0])

    for ai in 1..variants.high:
      if not variants[ai].get_samples(sample_fields, sample_cache[ai]):
        continue
      let asamples = sample_cache[ai]

      a = variants[ai].format.genotypes(x).alts
      cache[ai] = a
      for kid in kids:
        if sample_fields.len > 0 and kid.id notin asamples: continue
        # quick checks to rule out this variant.
        if a[kid.i] != 1: continue
        if not allow_non_trios:
          var parent_alleles = a[kid.mom.i] + a[kid.dad.i]
          if parent_alleles != 0 and parent_alleles != 1: continue

        for bi in 0..<ai:
          doAssert ai != bi
          b = cache[bi]

          let bsamples = sample_cache[bi]
          #if not variants[bi].get_samples(sample_fields, bsamples): continue
          if sample_fields.len > 0 and kid.id notin bsamples: continue
          if not is_compound_het(kid, a, b, allow_non_trios): continue

          if variants[ai].key notin foundKeys:
            found.add(variants[ai])
            foundKeys.incl(variants[ai].key)
          if variants[bi].key notin foundKeys:
            found.add(variants[bi])
            foundKeys.incl(variants[bi].key)

          var ab_key = variants[ai].key & "||" & variants[bi].key & "||" & kid.id
          if ab_key in foundKeys:
            continue
          foundKeys.incl(ab_key)

          cnt.inc(@[kid.id, kid.id], "compound-het")
          variants[ai].add_comphet(variants[bi], gene, kid.id, comphet_id)
          variants[bi].add_comphet(variants[ai], gene, kid.id, comphet_id)
          comphet_id += 1

  sort(found, proc (a:Variant, b:Variant): int =
    if a.start == b.start:
      return a.REF.len - b.REF.len
    return int((a.start - b.start)))

  for v in found:
    doAssert ovcf.write_variant(v)
    result.inc

proc main*(dropfirst:bool=false) =
  var p = newParser("slivar compound-hets"):
    help("find compound-hets in trios from pre-filtered variants")
    option("-v", "--vcf", default="/dev/stdin", help="input VCF")
    option("-s", "--sample-field", multiple=true, help="optional INFO field(s) that contains list of samples (kids) that have passed previous filters.\ncan be specified multiple times. this is needed for multi-family VCFs")
    option("-p", "--ped", default="", help="required ped file describing the trios in the VCF")
    option("--skip", default="splice_region,intergenic_region,intron,non_coding_transcript,non_coding,upstream_gene,downstream_gene,non_coding_transcript_exon,NMD_transcript,5_prime_UTR,3_prime_UTR", help="skip variants with these impacts (comma-separated)")
    option("-o", "--out-vcf", default="/dev/stdout", help="path to output VCF/BCF")
    flag("-a", "--allow-non-trios", help="allow samples with one or both parent unspecified. if this mode is used, any pair of heterozygotes co-occuring in the same gene, sample will be reported for samples without both parents that don't have kids. if a single parent is present some additional filtering is done.")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "compound-hets":
    argv = argv[1..argv.high]

  let opts = p.parse(argv)
  if opts.help:
    quit 0
  if opts.ped == "":
    echo p.help
    quit "--ped is required"


  var
    samples = parse_ped(opts.ped)
    ovcf:VCF
    ivcf:VCF
    skip_impacts = opts.skip.toLowerAscii.strip.split(',')

  block:
    var toAdd:seq[string]
    for sk in skip_impacts:
      if sk.endsWith("_variant"): toAdd.add(sk[0..sk.high-8])
      else: toAdd.add(sk & "_variant")
      if toAdd.len > 0 and toAdd[toAdd.high] in skip_impacts:
        discard toAdd.pop
    skip_impacts.add(toAdd)


  if not open(ivcf, opts.vcf, threads=2):
    quit "couldn't open vcf:" & opts.vcf

  samples = samples.match(ivcf)

  var kids = samples.kids(opts.allow_non_trios)
  if kids.len == 0:
    quit &"[slivar] no trios found for {opts.ped} with {opts.vcf}"
  var counter = initCounter(kids)

  var gene_fields: seq[GeneIndexes]
  for f in ["ANN", "CSQ", "BCSQ", getEnv("CSQ_FIELD")]:
    if f == "": continue
    try:
      var gf:GeneIndexes
      if ivcf.set_csq_fields(f, gf).len > 0:
        gene_fields.add(gf)
      # add this to the field names so we can clear it as needed
    except KeyError:
      continue

  doAssert gene_fields.len > 0, &"[slivar] error! no gene-like field found in {opts.vcf}"

  if not open(ovcf, opts.out_vcf, mode="w"):
    quit "couldn't open output vcf"

  # NOTE: if this format changes, must change slivar tsv as well
  if ivcf.header.add_info("slivar_comphet", ".", "String", "compound hets called by slivar. format is sample/gene/id/chrom/pos/ref/alt where id is a unique integer indicating the compound-het pair.") != Status.OK:
      quit "error adding field to header"
  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header

  var
    last_rid = -1
    tbl: TableRef[string,seq[Variant]]
    csqs: string = ""
    ncsqs = 0
    nwritten = 0
    comphet_id:int = 0


  var imp_counts = initCountTable[string]()
  var nvariants = -1
  for v in ivcf:
    nvariants += 1
    if v.rid != last_rid:
      if last_rid != -1:
        if nvariants > 5000 and tbl.len == 0:
          stderr.write_line &"[slivar] warning: no genes found for chromosome preceding rid:{v.rid} check your gene annotations have succeeded and that you are pulling the correct sample field(s)"
        nwritten += ovcf.write_compound_hets(kids, tbl, opts.sample_field, comphet_id, opts.allow_non_trios, counter)

      tbl = newTable[string, seq[Variant]]()
      last_rid = v.rid

    for g in gene_fields:
      if v.info.get(g.csq_field, csqs) != Status.OK or csqs.len == 0:
        continue

      # if there are variants without any of the sample fields, we can skip them.
      if opts.sample_field.len > 0:
        var has_any = false
        var samples = ""
        for f in opts.sample_field:
          if v.info.get(f, samples) == Status.OK and samples.len > 0:
            has_any = true
            break
        if not has_any: continue
      else:
        stderr.write_line "[slivar compound-het] WARNING: no sample fields specified, using only genotypes to call compound hets"


      var seen = initHashSet[string]()
      for csq in csqs.split(','):
        var fields = csq.split('|')
        if max(g.gene, g.consequence) >= fields.len:
          stderr.write_line &"[slivar] warning! {g.csq_field} has a CSQ of {csq} which is incomplete. skipping {v.CHROM}:{v.start + 1}({v.REF}/{v.ALT[0]})"
          continue
        var gene = fields[g.gene]
        if gene == "": continue
        if gene in seen: continue
        var impact_ok = false
        for imp in fields[g.consequence].split('&'):
          if imp.toLowerAscii notin skip_impacts:
            impact_ok = true
            break
          imp_counts.inc(imp)
        if not impact_ok: continue
        if gene notin tbl:
          tbl[gene] = @[v.copy()]
        else:
          tbl[gene].add(v.copy())
        ncsqs.inc
        seen.incl(gene)

  if last_rid != -1:
    nwritten += ovcf.write_compound_hets(kids, tbl, opts.sample_field, comphet_id, opts.allow_non_trios, counter)

  ovcf.close()
  ivcf.close()
  stderr.write_line &"[slivar compound-hets] wrote {nwritten} variants that were part of a compound het."

  var summaryPath = getEnv("SLIVAR_SUMMARY_FILE")

  if summaryPath == "":
    stderr.write_line counter.tostring(kids)
  else:
    var fh: File
    if not open(fh, summaryPath, fmWrite):
      quit "[slivar] couldn't open summary file:" & summaryPath
    fh.write(counter.tostring(kids))
    fh.close()
    stderr.write_line "[slivar] wrote summary table to:" & summaryPath

  if ncsqs == 0:
    stderr.write_line &"[slvar compound-hets] WARNING!!! no variants had any of the requested fields; unable to call compound hets!"
    quit 0

when isMainModule:
  import unittest
  suite "compound hets":
    test "simple":
      var kid = Sample(id:"kid", i:0)
      kid.dad = Sample(id:"dad", i: 1)
      kid.mom = Sample(id:"mom", i: 2)

      check is_compound_het(kid, [1'i8, 0, 1], [1'i8, 1, 0])

      check not is_compound_het(kid, [1'i8, 0, 1], [1'i8, 0, 1])
      check not is_compound_het(kid, [1'i8, 1, 1], [1'i8, 0, 1])

      check not is_compound_het(kid, [0'i8, 0, 1], [1'i8, 1, 0])

    test "with denovo":
      var kid = Sample(id:"kid", i:0)
      kid.dad = Sample(id:"dad", i: 1)
      kid.mom = Sample(id:"mom", i: 2)

      # first is de novo, 2nd is inherited
      check is_compound_het(kid, [1'i8, 0, 0], [1'i8, 1, 0])
      check not is_compound_het(kid, [1'i8, 1, 0], [1'i8, 1, 0])

      check is_compound_het(kid, [1'i8, 0, 1], [1'i8, 0, 0])


    test "without mom & dad":
      var kid = Sample(id:"kid", i:0)
      check is_compound_het(kid, [1'i8], [1'i8], allow_non_trios=true)


    test "kid with mom and dad without index":

      var mom = Sample(id:"mom", i: -1)
      var dad = Sample(id:"dad", i: -1)
      var kid = Sample(id:"kid", i: 0)
      kid.dad = dad
      kid.mom = mom

      check is_compound_het(kid, [1'i8], [1'i8], allow_non_trios=true)


