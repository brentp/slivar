# data-driven cutoffs
import times
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

const tmpl_html* = staticRead("ddc.html")

type Trace = object
  y: seq[float32]
  x: seq[float32]
  text: seq[string]
  name: string
  lt: bool
  abs: bool
  size_range: array[2, int]

  trues: seq[int32]
  falses: seq[int32]
  cutoffs: seq[float32]
  false_count: int32
  true_count: int32



type BoolExpr = object
  field: string
  fn: proc(value:float32): bool

type f32cmp = proc(a, b:float32): bool

proc lt(a: float32, b: float32): bool = return a < b
proc lte(a: float32, b: float32): bool = return a <= b
proc gt(a: float32, b: float32): bool = return a > b
proc gte(a: float32, b: float32): bool = return a >= b
proc eq(a: float32, b: float32): bool = return abs(a - b) < 1e-7

proc makeClosure(cmp: float32, op: f32cmp): proc(a: float32): bool =
  result = proc(a:float32): bool =
    let v = cmp
    return op(a,  v)

proc make_expressions(inputs: seq[string]): seq[BoolExpr] =
  var re = re"\s*([<=>]{1,2})\s*"
  for s in inputs:
    var parts:seq[string] = s.splitIncl(re)
    if parts.len != 3:
      quit &"expected expressions with 3 parts like DP<=33, got: {s}"

    var ocmp = parseFloat(parts[2])

    case parts[1]:
      of "<":
          result.add(BoolExpr(fn: makeClosure(ocmp, lt)))
      of "<=":
          result.add(BoolExpr(fn: makeClosure(ocmp, lte)))
      of ">":
          result.add(BoolExpr(fn: makeClosure(ocmp, gt)))
      of ">=":
          result.add(BoolExpr(fn: makeClosure(ocmp, gte)))
      of "==":
          result.add(BoolExpr(fn: makeClosure(ocmp, eq)))
      else:
        quit &"expected numeric expression got: {parts[1]}"
    result[result.high].field = parts[0]

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

proc getINFOf32(v:Variant, f:string, takeABS:bool): seq[float32] =
  var st = v.info.get(f, result)
  if st == Status.OK:
    if takeABS:
      for v in result.mitems: v = v.abs
    return
  if st != Status.UnexpectedType: return
  var i: seq[int32]
  st = v.info.get(f, i)
  if st != Status.OK: return
  result.setLen(i.len)
  for k, v in i:
    result[k] = v.float32
  if takeABS:
    for v in result.mitems: v = v.abs

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
  return alts[kid.i] == 1 and alts[kid.mom.i] == 0 and alts[kid.dad.i] == 0

proc inherited(kid:Sample, alts: seq[int8]): bool =
  return alts[kid.i] == 1 and [alts[kid.mom.i], alts[kid.dad.i]] in [[0'i8, 1], [1'i8, 0]]

type TF = object
  inherited: seq[float32]
  violation: seq[float32]
  flip: bool
  abs: bool

proc passes(kid:Sample, exprs: seq[BoolExpr], fields: seq[seq[float32]]): bool =
  for k, f in exprs:
    if not f.fn(fields[k][kid.i]):
      return false
  return true

proc should_skip(v:Variant, info_exprs: seq[BoolExpr]): bool =
  result = false
  for expr in info_exprs:
    var f = v.getINFOf32(expr.field, false)
    if f.len == 0: continue # NOTE: allowing variants without the field
    if not expr.fn(f[0]): return true

proc `$$`(k:float32): string {.inline.} =
  result = &"{k:.2f}"

proc get_variant_length(v:Variant): float32 =
  var length = float32(v.REF.len - v.ALT[0].len)
  if v.ALT[0][0] == '<':
    var lengths: seq[int32]
    if v.info.get("SVLEN", lengths) == Status.OK:
      length = lengths[0].float32
    else:
      length = float32(v.stop - v.start - 1)
      var svt:string
      if v.info.get("SVTYPE", svt) == Status.OK and svt == "DEL":
        length = -length
  result = length

proc write_fields(fh:File, infos: seq[float32], info_fields: seq[string], kid:Sample, fmt_hr_fields: seq[seq[float32]], hr_exprs: seq[BoolExpr], fmt_het_fields: seq[seq[float32]], het_exprs: seq[BoolExpr], violation:bool, header_written: var bool) =
  if not header_written:
    var strings: seq[string]
    var used: seq[string]
    for h in hr_exprs:
      if h.field in used: continue
      used.add(h.field)
      for s in ["kid", "dad", "mom"]:
        strings.add(s & "_" & h.field)
    for h in het_exprs:
      if h.field in used: continue
      used.add(h.field)
      for s in ["kid", "dad", "mom"]:
        strings.add(s & "_" & h.field)
    fh.write_line(&"""TP	{join(info_fields, "\t")}	{join(strings, "\t")}""")
    header_written = true
  var strings:seq[string]
  var used: seq[string]
  for i in infos:
    strings.add($$i)
    if strings[strings.high] in ["nan", "-nan"]: return
  for i, f in fmt_hr_fields:
    if hr_exprs[i].field in used: continue
    used.add(hr_exprs[i].field)
    strings.add($$f[kid.i])
    strings.add($$f[kid.dad.i])
    strings.add($$f[kid.mom.i])
  for i, f in fmt_het_fields:
    if het_exprs[i].field in used: continue
    used.add(het_exprs[i].field)
    strings.add($$f[kid.i])
    strings.add($$f[kid.dad.i])
    strings.add($$f[kid.mom.i])

  fh.write_line(&"""{1 - int(violation)}	{join(strings, "\t")}""")

proc ddc_main*() =
  var p = newParser("slivar ddc"):
    #option("-x", help="haploid (x) chromosome", default="chrX")
    option("--exclude", help="file of exclude regions")
    option("--chrom", help="limit to this chromosome only", default="chr15")
    option("--fmt-fixed-het", help="cutoff for sample field for hets. e.g. AB >= 0.2. only variants passing these filters are evaluated", multiple=true)
    option("--fmt-fixed-hom-ref", help="cutoff for sample field for hom-ref. e.g. AB <= 0.01. only variants passing these filters are evaluated", multiple=true)
    option("--info-fixed", help="cutoff for INFO field for hets. e.g. QD >= 8.0. only variants passing these filters are evaluated", multiple=true)
    #option("--fmt", help="FORMAT field(s) to vary. if prefixed with ^ then a condition is flipped for het and hom-ref, e.g. for AB, we expect higher for het and lower for hom-ref", multiple=true)
    option("--info", help="INFO field(s) to vary. these must be of Type Integer,Float, or Flag with Number=1 or A (or 0 for Flag). by default, higher is 'better', e.g. for GQ. prefix the field with '^' if lower is 'better' e.g. for FS", multiple=true)
    arg("vcf")
    arg("ped")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "ddc":
    argv = argv[1..argv.high]
  if len(argv) == 0: argv = @["--help"]

  let opts = p.parse(argv)
  if opts.help:
    quit 0

  var x: seq[int32]
  var exclude: TableRef[string, Lapper[region]]

  if opts.exclude != "":
    var t = cpuTime()
    exclude = read_bed(opts.exclude)
    stderr.write_line &"time to read exclude file: {cpuTime() - t:.2f}"

  var fmt_het_exprs = opts.fmt_fixed_het.make_expressions
  var fmt_hr_exprs = opts.fmt_fixed_hom_ref.make_expressions

  var info_exprs = opts.info_fixed.make_expressions

  var ivcf:VCF
  if not ivcf.open(opts.vcf):
    quit &"couldn't open {opts.vcf}"

  var samples = parse_ped(opts.ped).match(ivcf)
  var info_fields: seq[string] = opts.info
  #var fmt_fields: seq[string] = opts.fmt

  var info_values = newTable[string, TF]()
  var sizes_tbl = newTable[string, TF]()
  for i in info_fields.mitems:
    if i  == "QUAL":
      info_values[i] = TF()
      sizes_tbl[i] =  TF()
      continue
    try:
      discard ivcf.header.get(i.strip(chars={'^', '+'}), BCF_HEADER_TYPE.BCF_HL_INFO)
    except KeyError:
      quit &"requested info field {i} not found in header"
    info_values[i.strip(chars={'^', '+'})] = TF(flip:'^' in i, abs: '+' in i)
    sizes_tbl[i.strip(chars={'^', '+'})] =  TF()
    i = i.strip(chars={'^', '+'})


  var kids = samples.trio_kids
  var header_written: bool = false
  var last_rid = -1
  var last_lapper: Lapper[region]
  var vi = 0
  var ofh: File
  if not ofh.open("slivar_ddc_table.tsv", fmWrite):
    quit "could not open table"

  for v in ivcf.query(opts.chrom):
    if vi mod 500_000 == 0:
      stderr.write_line &"[slivar] at variant {v.CHROM}:{v.POS} ({vi})"
    vi.inc
    if v.FILTER != "PASS": continue
    if v.ALT[0] == "*": continue

    if v.rid != last_rid:
      if v.CHROM in ["chrX".cstring, "X".cstring]: break
      last_rid = v.rid
      if exclude != nil and stripChr(v.CHROM) in exclude:
        last_lapper = exclude[stripChr(v.CHROM)]
      else:
        var empty: seq[region]
        last_lapper = lapify(empty)

    if last_lapper.len > 0 and last_lapper.count(v.start, v.stop) != 0: continue


    var het_fields = newSeq[seq[float32]]()
    var hr_fields = newSeq[seq[float32]]()

    var violations = newSeq[bool](samples.len)
    var passes = newSeq[bool](samples.len)
    var infos = newSeqOfCap[float32](info_fields.len)

    if v.should_skip(info_exprs): continue
    # collect all the format fields once and then access then below for the
    # trios.
    for f in fmt_het_exprs:
      het_fields.add(v.getf32(f.field))
    for f in fmt_hr_exprs:
      hr_fields.add(v.getf32(f.field))

    var alts = v.format.genotypes(x).alts
    for i, e in info_fields:
      var vals: seq[float32]
      if e == "QUAL":
        vals.add(v.QUAL.float32)
      else:
        vals = v.getINFOf32(e, info_values[e].abs)
      # if this field is not in the variant, we can't do anything
      infos.add(if vals.len > 0: vals[0] else: Nan)
      if vals.len == 0: continue
      var vlen = v.get_variant_length.float32

      for kid in kids:
        if alts[kid.i] != 1: continue

        if kid.violation(alts):
          if not kid.passes(fmt_het_exprs, het_fields): continue
          if not kid.mom.passes(fmt_hr_exprs, hr_fields): continue
          if not kid.dad.passes(fmt_hr_exprs, hr_fields): continue
          info_values[e].violation.add(vals[0])
          sizes_tbl[e].violation.add(vlen)
          violations[kid.i] = true
          passes[kid.i] = true

        elif kid.inherited(alts):
          if not kid.passes(fmt_het_exprs, het_fields): continue
          if alts[kid.mom.i] == 1 and not kid.mom.passes(fmt_het_exprs, het_fields): continue
          elif alts[kid.mom.i] == 0 and not kid.mom.passes(fmt_hr_exprs, hr_fields): continue
          if alts[kid.dad.i] == 1 and not kid.dad.passes(fmt_het_exprs, het_fields): continue
          elif alts[kid.dad.i] == 0 and not kid.dad.passes(fmt_hr_exprs, hr_fields): continue
          info_values[e].inherited.add(vals[0])
          sizes_tbl[e].inherited.add(vlen)
          passes[kid.i] = true

    for kid in kids:
      if not passes[kid.i]: continue
      ofh.write_fields(infos, info_fields, kid, hr_fields, fmt_hr_exprs, het_fields, fmt_het_exprs, violations[kid.i], header_written)


    # todo: handle flags
    # todo: what about hom-ref calls?

  var traces = newTable[string, seq[Trace]]()
  var cuts = [0,1,2,3,4,5,6,8,10,20,30,40,50,100,300,1_000,10_000,int.high]
  stderr.write_line "collected variants, now creating ROC curves by size class"


  for info_name, seqs in info_values.mpairs:
      traces[info_name] = newSeq[Trace]()
      var sizes = sizes_tbl[info_name]

      doAssert seqs.inherited.len == sizes.inherited.len, &"{info_name}: {seqs.inherited.len} vs {sizes.inherited.len}"
      doAssert seqs.violation.len == sizes.violation.len, &"{info_name}: {seqs.inherited.len} vs {sizes.inherited.len}"

      var trs = newSeq[Trace]()
      # >= size_range[0], < size_range[1]
      trs.add(Trace(name: info_name, lt: seqs.flip, abs: seqs.abs, size_range: [0, 1]))
      for i, c in cuts:
        if i == 0: continue
        if i == cuts.high: break
        trs.add(Trace(name: info_name, lt: seqs.flip, abs: seqs.abs, size_range: [c, cuts[i+1]]))
        trs.add(Trace(name: info_name, lt: seqs.flip, abs: seqs.abs, size_range: [-c, -cuts[i - 1]]))
      trs.add(Trace(name: info_name, lt: seqs.flip, abs: seqs.abs, size_range: [-cuts[cuts.high], - cuts[cuts.high-1]]))

      var inh_all = seqs.inherited
      var vio_all = seqs.violation
      if vio_all.len == 0 or inh_all.len == 0: continue

      let lo = min(min(inh_all), min(vio_all))
      let hi = max(max(inh_all), max(vio_all))
      var step = (hi - lo) / 1000'f32 # 1k steps should be more than enough.
      # get the step for the entire set.
      # now partition into sizes
      for tr in trs.mitems:

        var inh: seq[float32]
        var vio: seq[float32]
        for i, v in inh_all:
          var sz = sizes.inherited[i].int32
          if sz >= tr.size_range[0] and sz < tr.size_range[1]:
            inh.add(v)

        for i, v in vio_all:
          var sz = sizes.violation[i].int32
          if sz >= tr.size_range[0] and sz < tr.size_range[1]:
            vio.add(v)

        var cutoff = lo - step - 1e-10
        while cutoff <= hi + step:
          cutoff += step
          var
            ti: int
            vi: int

          if seqs.flip:
            ti = lowerBound(inh, cutoff)
            vi = lowerBound(vio, cutoff)
          else:
            ti = inh.len - lowerBound(inh, cutoff)
            vi = vio.len - lowerBound(vio, cutoff)

          #var tpr = ti.float32 / inh.len.float32
          #var fpr = vi.float32 / vio.len.float32

          #tr.x.add(fpr)
          #tr.y.add(tpr)
          tr.trues.add(ti.int32)
          tr.falses.add(vi.int32)
          tr.false_count = vio.len.int32
          tr.true_count = inh.len.int32
          tr.cutoffs.add(cutoff)

        traces[info_name].add(tr)

  var s = %* traces
  var fh:File
  if not fh.open("slivar_ddc_roc.html", mode=fmWrite):
    quit "couldn't open output file."

  var sc = %* cuts
  fh.write(tmpl_html.replace("<SIZES>", $sc).replace("<INPUT_JSON>", $s))
  fh.close()

when isMainModule:
  ddc_main()
