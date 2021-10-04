# data-driven cutoffs
import times
import random
import hts/vcf
import json
import lapper
import algorithm
import pedfile
import tables
import argparse
import strutils
import strformat
import os

proc `$$`(k:float32): string {.inline.} =
  result = &"{k:.3f}"
  while result.len > 1 and result[result.high] == '0':
    result = result.strip(chars={'0'}, leading=false)
  if result[result.high] == '.':
    result = result.strip(chars={'.'}, leading=false)

const tmpl_html* = staticRead("ddc.html")

type Trio = object
  sample_id: string
  tbl: Table[string, seq[float32]]
  #kid_alts: seq[int8]
  #dad_alts: seq[int8]
  #mom_alts: seq[int8]
  variant_idxs: seq[uint32]

  violations: seq[bool]

type VariantInfo = object
  float_tbl: Table[string, seq[float32]]
  bool_tbl: Table[string, seq[bool]]
  variant_lengths: seq[int]
  filters: seq[string]
  violations: seq[bool]

proc tojson(tbl: Table[string, seq[float32]]): string =
  result = newStringOfCap(16384)
  result.add('{')
  for k, vals in tbl:
    result.add(k)
    result.add(":[")
    for v in vals:
      result.add($$v)
      result.add(',')
    if len(vals) == 0: result.add(' ')
    result[result.high] = ']'
    result.add(",\n")
  if result.len > 1:
    result[result.high-1] = '\n'
    result[result.high] = '}'
  else:
    stderr.write_line "[slivar ddc] no values found"
    result = "{}"


proc tojson(t:Trio): string =
  result = newStringOfCap(16384)
  result.add(&"""{{sample_id: "{t.sample_id}",
tbl: {t.tbl.tojson},
violations: {%t.violations},
variant_idxs: {%t.variant_idxs}
}}""")

proc tojson(ts:seq[Trio]): string =
  result = newStringOfCap(16384*128)
  result.add('[')
  for t in ts:
    result.add(t.tojson)
    result.add(',')
  result[result.high] = ']'

proc tojson(v:VariantInfo): string =
  result = newStringOfCap(16384)
  result.add(&"""{{float_tbl: {v.float_tbl.tojson},
bool_tbl: {%v.bool_tbl},
variant_lengths: {%v.variant_lengths},
filters: {%v.filters},
violations: {%v.violations}
}}""")

proc trio_kids(samples: seq[Sample]): seq[Sample] =
  ## return all samples that have a mom and dad.
  result = newSeqOfCap[Sample](16)
  for sample in samples:
    if sample.mom == nil or sample.mom.i == -1 or sample.dad == nil or sample.dad.i == -1: continue
    result.add(sample)

proc getAB(v:Variant): seq[float32] =
  if v.format.get("AB", result) != Status.OK:

    var ad: seq[int32]
    if v.format.get("AD", ad) != Status.OK:
      return
    result = newSeq[float32](v.n_samples)
    for i in 0..<v.n_samples:
      result[i] = ad[2*i+1].float32 / max(1, ad[2*i+1] + ad[2*i]).float32
  for ab in result.mitems:
    if ab < 0: ab = 0
    if ab > 1: ab = 1

  # flip so that the alllele balance for hets is always < 0.5
  #for i, ab in result.mpairs:
  #  if ab == 0'f32 or ab == 1'f32: continue
  #  if ab > 0.5: ab = 1-ab

proc getINFOf32(v:Variant, f:string): float32 =
  if f == "QUAL":
    return v.QUAL.float32
  var rr = newSeqUninitialized[float32](1)
  var st = v.info.get(f, rr)
  if st == Status.OK:
    doAssert rr.len == 1, &"[slivar] require only fields with single values got {rr.len} in field {f}"
    return rr[0]
  if st != Status.UnexpectedType: return
  var i: seq[int32]
  st = v.info.get(f, i)
  if st != Status.OK: return
  doAssert i.len == 1, &"[slivar] require only fields with single values got {rr.len} in field {f}"
  return i[0].float32

proc getf32(v:Variant, f:string): seq[float32] =
  if f == "AB": return v.getAB
  var st = v.format.get(f, result)
  if st == Status.OK: return
  if st != Status.UnexpectedType: return
  var i: seq[int32]
  st = v.format.get(f, i)
  if st == Status.NotFound:
    # get GQ from PL
    if v.format.get("PL", i) != Status.OK: return
    var tmp = newSeq[int32](int(i.len/3))
    for idx in 0..<tmp.len:
      var pls = @[i[3*idx], i[3*idx+1], i[3*idx+2]]
      sort(pls)
      tmp[idx] = pls[1] - pls[0]
    i = tmp

  result.setLen(i.len)
  for k, val in i:
    # TODO: handle missing
    result[k] = val.float32

proc violation(kid:Sample, alts: seq[int8], allele_balances: seq[float32]): bool =
  result = alts[kid.i] >= 1 and alts[kid.mom.i] == 0 and alts[kid.dad.i] == 0
  if not result and alts[kid.i] == 2:
      result = [alts[kid.mom.i], alts[kid.dad.i]] in [[0'i8, 1], [1'i8, 0]]
  if not result and alts[kid.i] == 0:
      result = [alts[kid.mom.i], alts[kid.dad.i]] in [[2'i8, 1], [1'i8, 2]]

  if allele_balances.len == 0 or result == false: return
  if allele_balances[kid.mom.i] > 0 or allele_balances[kid.dad.i] > 0:
    raise newException(ValueError, "non zero allele balance for parents")

proc inherited(kid:Sample, alts: seq[int8], allele_balances: seq[float32]): bool =
  result = alts[kid.i] == 1 and [alts[kid.mom.i], alts[kid.dad.i]] in [[0'i8, 1], [1'i8, 0]]
  if result == false or allele_balances.len == 0: return
  if alts[kid.dad.i] == 0 and allele_balances[kid.dad.i] > 0:
    raise newException(ValueError, "non zero allele balance for parents")
  if alts[kid.mom.i] == 0 and allele_balances[kid.mom.i] > 0:
    raise newException(ValueError, "non zero allele balance for parents")

proc get_variant_length(v:Variant): int =
  var length = int(v.ALT[0].len - v.REF.len)
  if v.ALT[0][0] == '<':
    var lengths: seq[int32]
    if v.info.get("SVLEN", lengths) == Status.OK:
      length = lengths[0]
    else:
      length = int(v.stop - v.start - 1)
      var svt:string
      if v.info.get("SVTYPE", svt) == Status.OK and svt == "DEL":
        length = -length
  result = length

proc check*[T: VariantInfo|seq[Trio]](ivcf:VCF, fields: seq[string], ftype:BCF_HEADER_TYPE, infos: var T): seq[string] =
  for f in fields:
    try:
      let hr = ivcf.header.get(f, ftype)
      when T is VariantInfo:

        case hr["Type"]
        of "Float", "Integer":
          infos.float_tbl[f] = newSeqOfCap[float32](65536)
        of "Flag":
          infos.bool_tbl[f] = newSeqOfCap[bool](65536)
        else:
          echo "unknown type for:", f, " ", hr

      else:

        case hr["Type"]
        of "Float", "Integer":
          for s in infos.mitems:
            s.tbl[f] = newSeqOfCap[float32](65536)
        else:
          quit "only float and integer types for supported for format fields. got:" & $hr

    except KeyError:
      when T is VariantInfo:
        if f == "QUAL":
          infos.float_tbl[f] = newSeqOfCap[float32](65536)
        else:
          quit &"requested info field {f} not found in header"
      else:
        if f == "AB":
          for s in infos.mitems:
            s.tbl[f] = newSeqOfCap[float32](65536)
        else:
          quit &"requested info field {f} not found in header"
  return fields

proc ddc_main*(dropfirst:bool=false) =
  var p = newParser("slivar ddc"):
    #option("-x", help="haploid (x) chromosome", default="chrX")
    option("--chrom", help="limit to this chromosome only. use '-3' for all chromosomes (in the case of exome data)", default="chr15")
    option("--info-fields", help="comma-delimited list of info fields")
    option("--fmt-fields", help="comma-delimited list of sample fields")
    option("--html", default="slivar-ddc.html", help="path to output file")
    arg("vcf")
    arg("ped")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "ddc":
    argv = argv[1..argv.high]
  if len(argv) == 0: argv = @["--help"]

  let opts = p.parse(argv)
  if opts.help:
    quit 0

  var ivcf:VCF
  if not open(ivcf, opts.vcf, threads=2):
    quit "could not open vcf:" & opts.vcf

  var samples = parse_ped(opts.ped).match(ivcf)
  var kids = samples.trio_kids
  const max_trios = 4
  if len(kids) > max_trios:
    stderr.write_line &"[slivar] sub-sampling to {max_trios} random kids"
    randomize()
    shuffle(kids)
    kids = kids[0..<max_trios]
    kids.sort(proc(a, b:Sample): int = a.i - b.i)
    samples = newSeq[Sample]()
    for k in kids:
      samples.add(k)
      if k.dad notin samples:
        samples.add(k.dad)
      if k.mom notin samples:
        samples.add(k.mom)

  var output_infos = VariantInfo(float_tbl: initTable[string, seq[float32]](), bool_tbl: initTable[string, seq[bool]]())
  var output_trios = newSeq[Trio](kids.len)
  for i, o in output_trios.mpairs:
    o.sample_id = kids[i].id
    o.tbl = initTable[string, seq[float32]]()

  if opts.info_fields != "":
    discard ivcf.check(opts.info_fields.strip().split(','), BCF_HEADER_TYPE.BCF_HL_INFO, output_infos)
  var fmt_fields = ivcf.check(opts.fmt_fields.strip().split(','), BCF_HEADER_TYPE.BCF_HL_FMT, output_trios)

  var x: seq[int32]

  var ab: seq[float32]

  var variant_idx = 0
  var skip_xy = opts.chrom in ["*", "-3"]
  for v in ivcf.query(opts.chrom):
    if skip_xy:
      let CHROM = v.CHROM
      if CHROM[CHROM.high] in "XY": continue

    var fmts = newSeq[seq[float32]](fmt_fields.len)
    for i, f in fmt_fields:
      fmts[i] = v.getf32(f)
    shallow(fmts)
    if "AB" in fmt_fields:
      ab = fmts[fmt_fields.find("AB")]

    var alts = v.format.genotypes(x).alts
    shallow(alts)
    var any_used = false
    var any_violation = false

    for i, kid in kids:
      var vio: bool
      var inh: bool
      try:
        vio = kid.violation(alts, ab)
        inh = kid.inherited(alts, ab)
      except:
        continue
      if not (vio or inh): continue
      any_used = true
      if vio: any_violation = true

      var tr = output_trios[i]
      tr.variant_idxs.add((variant_idx).uint32)
      tr.violations.add(vio)
      #tr.alts.add(alts[kid.i])
      #tr.dad_alts.add(alts[kid.dad.i])
      #tr.mom_alts.add(alts[kid.mom.i])

      for k, kid_seq in tr.tbl.mpairs:

        let fi = fmt_fields.find(k)
        let fmt = fmts[fi]
        if fmt.len == 0: continue


        if k == "GQ":
          kid_seq.add(min(min(fmt[kid.i], fmt[kid.dad.i]), fmt[kid.mom.i]))
        else:
          kid_seq.add(fmt[kid.i])

      output_trios[i] = tr

    if not any_used: continue

    output_infos.filters.add(v.FILTER)
    output_infos.variant_lengths.add(v.get_variant_length)
    output_infos.violations.add(any_violation)
    for k, bseq in output_infos.bool_tbl.mpairs:
      bseq.add(v.info.has_flag(k))
    for k, fseq in output_infos.float_tbl.mpairs:
      fseq.add(v.getINFOf32(k))
    variant_idx.inc

  var html = tmpl_html.replace("<VARIANT_JSON>", output_infos.tojson)
  html = html.replace("<TRIO_JSON>", output_trios.tojson)


  var fh:File
  if not open(fh, opts.html, fmWrite):
    quit "couldn't open output html file:" & opts.html
  fh.write(html)
  stderr.write_line "[slivar] wrote output to:" & opts.html
  fh.close


  stderr.write_line "number of infos:", output_infos.variant_lengths.len
  for tr in output_trios:
    stderr.write_line "trio:", tr.variant_idxs.len

when isMainModule:
  ddc_main()
