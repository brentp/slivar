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
  for i in info_fields.mitems:
    if i == "QUAL":
      info_values[i] = TF()
      continue
    try:
      discard ivcf.header.get(i.strip(chars={'^', '+'}), BCF_HEADER_TYPE.BCF_HL_INFO)
    except KeyError:
      quit &"requested info field {i} not found in header"
    info_values[i.strip(chars={'^', '+'})] = TF(flip:'^' in i, abs: '+' in i)
    i = i.strip(chars={'^', '+'})


  var kids = samples.trio_kids
  var last_rid = -1
  var last_lapper: Lapper[region]

  for v in ivcf.query(opts.chrom):
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
      if vals.len == 0: continue
      for kid in kids:
        if alts[kid.i] != 1: continue

        if kid.violation(alts):
          if not kid.passes(fmt_het_exprs, het_fields): continue
          if not kid.mom.passes(fmt_hr_exprs, hr_fields): continue
          if not kid.dad.passes(fmt_hr_exprs, hr_fields): continue
          info_values[e].violation.add(vals[0])

        elif kid.inherited(alts):
          if not kid.passes(fmt_het_exprs, het_fields): continue
          if alts[kid.mom.i] == 1 and not kid.mom.passes(fmt_het_exprs, het_fields): continue
          elif alts[kid.mom.i] == 0 and not kid.mom.passes(fmt_hr_exprs, hr_fields): continue
          if alts[kid.dad.i] == 1 and not kid.dad.passes(fmt_het_exprs, het_fields): continue
          elif alts[kid.dad.i] == 0 and not kid.dad.passes(fmt_hr_exprs, hr_fields): continue
          info_values[e].inherited.add(vals[0])


    # todo: handle flags
    # todo: what about hom-ref calls?

  var traces: seq[Trace]

  for info_name, seqs in info_values.mpairs:

      var tr = Trace(name: info_name, lt: seqs.flip, abs: seqs.abs)
      var inh = seqs.inherited
      var vio = seqs.violation
      sort(inh)
      sort(vio)
      var last_tpr = float32.low
      var last_fpr = float32.low
      if vio.len == 0 or inh.len == 0: continue

      let lo = min(inh[0], vio[0])
      let hi = max(inh[^1], vio[^1])
      var step = (hi - lo) / 1000'f32 # 1k steps should be more than enough.
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
        #var ti = lowerBound(inh, cutoff)
        var tpr = ti.float32 / inh.len.float32

        #var vi = lowerBound(vio, cutoff)
        var fpr = vi.float32 / vio.len.float32

        if abs(fpr - last_fpr) < 2e-4 and abs(tpr - last_tpr) < 2e-4:
          continue
        last_fpr = fpr
        last_tpr = tpr

        tr.x.add(fpr)
        tr.y.add(tpr)
        let sign = if seqs.flip: '<' else: '>'
        let value = if seqs.abs: &"abs({info_name})": else: info_name
        tr.text.add(&"trues:{ti} falses:{vi} total:{vio.len + inh.len}<br>include:{value}{sign}{cutoff:.2f}")
        #tr.text.add(&"trues:{ti} falses:{vi} total:{total}, cutoff:{cutoff}")

      #for i, c in cutoffs:
      #  echo &"{info_name}:{c} tpr:{tprs[i]:.3f} fpr:{fprs[i]:.3f}"
      traces.add(tr)

  var s = %* traces
  var fh:File
  if not fh.open("roc.html", mode=fmWrite):
    quit "couldn't open output file."

  fh.write(tmpl_html.replace("<INPUT_JSON>", $s))
  fh.close()

when isMainModule:
  ddc_main()
