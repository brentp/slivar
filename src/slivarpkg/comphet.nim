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

proc add_ch(a:Variant, b:Variant) =
  var s: string = ""
  discard a.info.get("slivar_comphet", s)
  if s.len > 1:
    echo "found"
    s &= ","
  else:
    echo "not found"

  var balts = join(b.ALT, ",")
  var toadd = &"{$b.CHROM}/{b.start+1}/{b.REF}/{balts}"
  s &= toadd
  doAssert a.info.set("slivar_comphet", s) == Status.OK

proc write_compound_hets(ovcf:VCF, kids:seq[Sample], tbl:TableRef[string, seq[Variant]]) =
  var
    x: seq[int32]

  for gene, variants in tbl.mpairs:
    var variants = variants
    var found = initHashSet[int]()
    if variants.len < 2: continue
    shallow(variants)

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
          found.incl(ai)
          found.incl(bi)
          var av = variants[ai]
          var bv = variants[bi]
          av.add_ch(bv)
          bv.add_ch(av)
    if found.len == 0: continue
    var ordered = newSeq[int]()
    for i in found:
      ordered.add(i)
    sort(ordered, system.cmp)
    for o in ordered:
      doAssert ovcf.write_variant(variants[o])

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

  if ivcf.header.add_info("slivar_comphet", ".", "String", "compound hets called by slivar") != Status.OK:
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

    for csq in csqs.split(","):
      var fields = csq.split("|")
      var gene = fields[index]
      if gene == "": continue
      if gene notin tbl:
        tbl[gene] = @[v.copy()]
      else:
        tbl[gene].add(v.copy())
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

  main()
