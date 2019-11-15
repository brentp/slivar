# data-driven cutoffs
import times
import random
import hts/vcf
import math
import json
import lapper
import algorithm
import ./pedfile
import tables
import ./evaluator
import argparse
import ./utils
import strutils
import strformat
import regex
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
  kid_tbl: Table[string, seq[float32]]
  kid_alts: seq[int8]
  dad_tbl: Table[string, seq[float32]]
  dad_alts: seq[int8]
  mom_tbl: Table[string, seq[float32]]
  mom_alts: seq[int8]
  variant_idxs: seq[uint32]

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
    result[result.high] = ']'
    result.add(",\n")
  result[result.high-1] = '\n'
  result[result.high] = '}'

proc tojson(t:Trio): string =
  result = newStringOfCap(16384)
  result.add(&"""{{sample_id: "{t.sample_id}",
kid_tbl: {t.kid_tbl.tojson},
dad_tbl: {t.dad_tbl.tojson},
mom_tbl: {t.mom_tbl.tojson},
kid_alts: {%t.kid_alts},
dad_alts: {%t.dad_alts},
mom_alts: {%t.mom_alts},
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
  if v.format.get("AB", result) == Status.OK: return

  var ad: seq[int32]
  if v.format.get("AD", ad) != Status.OK:
    return
  result = newSeq[float32](v.n_samples)
  for i in 0..<v.n_samples:
    result[i] = ad[2*i+1].float32 / max(1, ad[2*i+1] + ad[2*i]).float32

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
  if st != Status.OK: return
  result.setLen(i.len)
  for k, v in i:
    result[k] = v.float32

proc violation(kid:Sample, alts: seq[int8]): bool =
  return alts[kid.i] >= 1 and alts[kid.mom.i] == 0 and alts[kid.dad.i] == 0

proc inherited(kid:Sample, alts: seq[int8]): bool =
  return alts[kid.i] == 1 and [alts[kid.mom.i], alts[kid.dad.i]] in [[0'i8, 1], [1'i8, 0]]

proc get_variant_length(v:Variant): int =
  var length = int(v.ALT[0].len - v.REF.len)
  if v.ALT[0][0] == '<':
    var lengths: seq[int32]
    if v.info.get("SVLEN", lengths) == Status.OK:
      length = lengths[0]
    else:
      length = v.stop - v.start - 1
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
          echo "unknown type for", f, " ", hr

      else:

        case hr["Type"]
        of "Float", "Integer":
          for s in infos.mitems:
            s.kid_tbl[f] = newSeqOfCap[float32](65536)
            s.mom_tbl[f] = newSeqOfCap[float32](65536)
            s.dad_tbl[f] = newSeqOfCap[float32](65536)
        else:
          quit "only float and integer types for supported for format fields. got:" & $hr

    except KeyError:
      when T is VariantInfo:
        quit &"requested info field {f} not found in header"
      else:
        if f == "AB":
          for s in infos.mitems:
            s.kid_tbl[f] = newSeqOfCap[float32](65536)
            s.mom_tbl[f] = newSeqOfCap[float32](65536)
            s.dad_tbl[f] = newSeqOfCap[float32](65536)
        else:
          quit &"requested info field {f} not found in header"
  return fields

proc ddc_main*() =
  var p = newParser("slivar ddc"):
    #option("-x", help="haploid (x) chromosome", default="chrX")
    option("--chrom", help="limit to this chromosome only", default="chr15")
    option("--info-fields", help="comma-delimited list of info fields")
    option("--fmt-fields", help="comma-delimited list of sample fields")
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
    o.kid_tbl = initTable[string, seq[float32]]()
    o.dad_tbl = initTable[string, seq[float32]]()
    o.mom_tbl = initTable[string, seq[float32]]()

  var info_fields = ivcf.check(opts.info_fields.split(','), BCF_HEADER_TYPE.BCF_HL_INFO, output_infos)
  var fmt_fields = ivcf.check(opts.fmt_fields.split(','), BCF_HEADER_TYPE.BCF_HL_FMT, output_trios)

  var x: seq[int32]
  var f32: seq[float32]
  var i32: seq[int32]

  var variant_idx = -1
  for v in ivcf.query(opts.chrom):
    variant_idx += 1

    var fmts = newSeq[seq[float32]](fmt_fields.len)
    for i, f in fmt_fields:
      fmts[i] = v.getf32(f)
    shallow(fmts)

    var alts = v.format.genotypes(x).alts
    shallow(alts)
    var any_used = false
    var any_violation = false

    for i, kid in kids:
      let vio = kid.violation(alts)
      let inh = kid.inherited(alts)
      if not (vio or inh): continue
      any_used = true
      if vio: any_violation = true

      var tr = output_trios[i]
      tr.variant_idxs.add((variant_idx).uint32)
      tr.kid_alts.add(alts[kid.i])
      tr.dad_alts.add(alts[kid.dad.i])
      tr.mom_alts.add(alts[kid.mom.i])

      for k, kid_seq in tr.kid_tbl.mpairs:

        let fi = fmt_fields.find(k)
        let fmt = fmts[fi]

        kid_seq.add(fmt[kid.i])
        tr.dad_tbl[k].add(fmt[kid.dad.i])
        tr.mom_tbl[k].add(fmt[kid.mom.i])
      output_trios[i] = tr

    if not any_used: continue

    output_infos.filters.add(v.FILTER)
    output_infos.variant_lengths.add(v.get_variant_length)
    output_infos.violations.add(any_violation)
    for k, bseq in output_infos.bool_tbl.mpairs:
      bseq.add(v.info.has_flag(k))
    for k, fseq in output_infos.float_tbl.mpairs:
      fseq.add(v.getINFOf32(k))

  var html = tmpl_html.replace("<VARIANT_JSON>", output_infos.tojson)
  html = html.replace("<TRIO_JSON>", output_trios.tojson)


  var fh:File
  if not open(fh, "ddc.html", fmWrite):
    quit "couldn't open output hmlt file"
  fh.write(html)
  fh.close


  stderr.write_line "number of infos:", output_infos.variant_lengths.len
  for tr in output_trios:
    stderr.write_line "trio:", tr.variant_idxs.len

when isMainModule:
  ddc_main()
