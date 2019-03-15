import hts/vcf
import bpbiopkg/pedfile
import strformat
import algorithm
import argparse
import tables
import sets

proc is_compound_het*(kid: pedfile.Sample, a:openarray[int8], b:openarray[int8]): bool {.inline.} =
  doAssert kid.dad != nil and kid.mom != nil
  if a[kid.i] != 1 or b[kid.i] != 1: return false

  if a[kid.dad.i] == -1 or a[kid.mom.i] == -1: return false
  if a[kid.dad.i] + a[kid.mom.i] != 1: return false

  if b[kid.dad.i] == -1 or b[kid.mom.i] == -1: return false
  if b[kid.dad.i] + b[kid.mom.i] != 1: return false
  # now we have 1 het parent at each site and kid het at both sites.
  # check that different parents have the different alleles.
  return a[kid.dad.i] + b[kid.dad.i] == 1 and a[kid.mom.i] + b[kid.mom.i] == 1

proc kids(samples:seq[Sample]): seq[Sample] =
  for s in samples:
    if s.dad == nil or s.mom == nil: continue
    result.add(s)

proc key(v:Variant): string {.inline.} =
  var balts = join(v.ALT, ",")
  &"{$v.CHROM}/{v.start+1}/{v.REF}/{balts}"

proc add_comphet(a:Variant, b:Variant, gene: string) =
  var s: string = ""
  discard a.info.get("slivar_comphet", s)
  if s.len > 1:
    # htslib doesn't allow setting an existing field again. so we delete and re-add.
    doAssert a.info.delete("slivar_comphet") == Status.OK
    s &= ","
  s &= b.key & "/" & gene
  doAssert a.info.set("slivar_comphet", s) == Status.OK

proc write_compound_hets(ovcf:VCF, kids:seq[Sample], tbl:TableRef[string, seq[Variant]]) =
  var
    x: seq[int32]

  var foundKeys = initHashSet[string]()
  # need to track variants in this because 1 variant can be a CH with multiple genes/transcripts.
  var found: seq[Variant]

  for gene, variants in tbl.mpairs:
    var variants = variants
    if variants.len < 2: continue

    for ai in 1..variants.high:
      var a = variants[ai].format.genotypes(x).alts
      for kid in kids:
        # quick checks to rule out this variant.
        if a[kid.i] != 1: continue
        if a[kid.mom.i] + a[kid.dad.i] != 1: continue

        for bi in 0..<ai:
          doAssert ai != bi
          var b = variants[bi].format.genotypes(x).alts
          if not is_compound_het(kid, a, b): continue

          if variants[ai].key notin foundKeys:
            found.add(variants[ai])
            foundKeys.incl(variants[ai].key)
          if variants[bi].key notin foundKeys:
            found.add(variants[bi])
            foundKeys.incl(variants[bi].key)

          variants[ai].add_comphet(variants[bi], gene)
          variants[bi].add_comphet(variants[ai], gene)

  sort(found, proc (a:Variant, b:Variant): int =
    if a.start == b.start:
      return a.REF.len - b.REF.len
    return a.start - b.start)

  for v in found:
    doAssert ovcf.write_variant(v)

proc main*(dropfirst:bool=false) =
  var p = newParser("slivar compound-hets"):
    help("find compound-hets in trios from pre-filtered variants")
    option("-v", "--vcf", default="/dev/stdin", help="input VCF")
    option("-p", "--ped", default="", help="required ped file describing the trios in the VCF")
    option("-f", "--field", default="BCSQ", help="INFO field containing the gene name")
    option("-i", "--index", default="2", help="(1-based) index of the gene-name in the field after splitting on '|'")

  let opts = p.parse()
  if opts.ped == "":
    echo p.help
    quit "--ped is required"

  var
    samples = parse_ped(opts.ped)
    ovcf:VCF
    ivcf:VCF

  if not open(ivcf, opts.vcf, threads=2):
    quit "couldn't open vcf:" & opts.vcf

  samples = samples.match(ivcf)
  var kids = samples.kids
  if kids.len == 0:
    quit &"[slivar] no trios found for {opts.ped} with {opts.vcf}"

  if not open(ovcf, "/dev/stdout", mode="w"):
    quit "couldn't open output vcf"

  if ivcf.header.add_info("slivar_comphet", ".", "String", "compound hets called by slivar. format is chrom/pos/ref/alt/gene") != Status.OK:
      quit "error adding field to header"
  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header

  var
    last_rid = -1
    tbl: TableRef[string,seq[Variant]]
    csqs: string = ""
    index = parseInt(opts.index) - 1

  for v in ivcf:
    if v.rid != last_rid:
      if last_rid != -1:
        ovcf.write_compound_hets(kids, tbl)

      tbl = newTable[string, seq[Variant]]()
      last_rid = v.rid

    if v.info.get(opts.field, csqs) != Status.OK or csqs.len == 0:
      continue

    var seen = initHashSet[string]()
    for csq in csqs.split(","):
      var fields = csq.split("|")
      var gene = fields[index]
      if gene == "": continue
      if gene in seen: continue
      if gene notin tbl:
        tbl[gene] = @[v.copy()]
      else:
        tbl[gene].add(v.copy())
      seen.incl(gene)

  if last_rid != -1:
    ovcf.write_compound_hets(kids, tbl)

  ovcf.close()
  ivcf.close()

when isMainModule:
  import unittest
  suite "compound hets":
    test "simple":
      var kid = Sample(id:"kid")
      kid.mom = Sample(id:"mom")
      kid.dad = Sample(id:"dad")
      kid.i = 0
      kid.dad.i = 1
      kid.mom.i = 2

      check is_compound_het(kid, [1'i8, 0, 1], [1'i8, 1, 0])

      check not is_compound_het(kid, [1'i8, 0, 1], [1'i8, 0, 1])
      check not is_compound_het(kid, [1'i8, 1, 1], [1'i8, 0, 1])

      check not is_compound_het(kid, [0'i8, 0, 1], [1'i8, 1, 0])
