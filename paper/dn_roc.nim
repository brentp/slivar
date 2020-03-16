import hts/vcf
import strformat
import argparse
import pedfile
import random
import slivarpkg/utils
import lapper
import tables
import strutils
import math
import algorithm

import os

var p = newParser("dn_roc"):
  option("-x", "--exclude", help="BED file of exclude regions (e.g. LCRs)")
  option("--gq", help="integer genotype quality cutoff between 0 and 99 inclusive", multiple=true)
  option("-d", "--min-depth", help="minimum depth", default="0")
  option("-q", "--min-variant-quality", help="variants with a quality lower than this are excluded", default="0")
  arg("ped")
  arg("vcf")

# TODO: lapper stuff
var args = commandLineParams()
if len(args) == 0: args = @["--help"]

var opts = p.parse(args)
if opts.help:
  quit 0

var samples = parse_ped(opts.ped)
var ivcf:VCF
if not open(ivcf, opts.vcf, threads=3):
  quit "couldn't open vcf"

var mvq = parseFloat(opts.min_variant_quality)

var gq_cutoffs = newSeq[int](opts.gq.len)
var depth_cutoff = parseInt(opts.min_depth)
for i, gq in opts.gq:
  var v = parseInt(gq)
  if v < 0 or v > 99:
    quit "genotype quality cutoff must be between 0 and 99"
  gq_cutoffs[i] = v
sort(gq_cutoffs)

samples = samples.match(ivcf)

var exclude: TableRef[string, Lapper[region]]
if opts.exclude != "":
  exclude = read_bed(opts.exclude)

var kids = newSeq[Sample]()
var counts = CountTable[string]()

for s in samples:
  if s.family_id == "1357": continue
  if s.dad != nil and s.mom != nil:
    if s.mom.i == -1 or (s.mom.id in counts and counts[s.mom.id] > 2): continue
    if s.dad.i == -1 or (s.dad.id in counts and counts[s.dad.id] > 2): continue
    counts.inc(s.mom.id)
    counts.inc(s.dad.id)
    kids.add(s)

stderr.write_line "kids:", $kids.len

proc ab(a:array[2, int32]): float32 {.inline.} =
  result = a[1].float32 / max(1, sum(a)).float32
  if result > 0.5: result = 1 - result

proc good(kid:Sample, AD: seq[int32], alts: seq[int8]): float32 {.inline.} =
  # variant is in kid and 1 and only 1 parent.
  # it is not unkknown in any sample in the trio.
  if alts[kid.i] != 1: return -1
  if not ((alts[kid.mom.i] == 1 and alts[kid.dad.i] == 0) or (alts[kid.mom.i] == 0 and alts[kid.dad.i] == 1)): return -1
  var p = if alts[kid.mom.i] == 1: kid.mom else: kid.dad

  var kida: array[2, int32]
  kida[0] = AD[2*kid.i]
  kida[1] = AD[2*kid.i+1]

  var parent: array[2, int32]
  parent[0] = AD[2*p.i]
  parent[1] = AD[2*p.i+1]

  return min(parent.ab, kida.ab)

proc bad(kid:Sample, AD: seq[int32], alts: seq[int8]): float64 {.inline.} =
  # variant is a violation
  if alts[kid.i] < 0 or alts[kid.dad.i] < 0 or alts[kid.mom.i] < 0: return -1

  #[
  if [alts[kid.i], alts[kid.mom.i], alts[kid.dad.i]] notin [
    [0'i8, 2, 2],
    [0'i8, 2, 1],
    [0'i8, 1, 2],
    [2'i8, 0, 0],
    [2'i8, 0, 1],
    [2'i8, 1, 0]]: return - 1
    ]#
  if [alts[kid.i], alts[kid.mom.i], alts[kid.dad.i]] notin [[1'i8, 0, 0], [1'i8, 2, 2]]:
    return -1

  var kida: array[2, int32]
  kida[0] = AD[2*kid.i]
  kida[1] = AD[2*kid.i+1]
  return kida.ab

# outer index is genotype quality
var goodABs = newSeq[seq[float32]](101)
var badABs = newSeq[seq[float32]](101)

shallow(goodABs)
shallow(badABs)


var AD: seq[int32]
var GQ: seq[int32]
var x: seq[int32]

goodABs[gq_cutoffs[0]] = newSeqOfCap[float32](65536*1048)

var empty: seq[region]
var last_tid = -1
for v in ivcf:
  if v.QUAL < mvq: continue
  if v.rid != last_tid:
    last_tid = v.rid
    if $v.CHROM in ["X", "chrX"]: break
    stderr.write_line "at:", $v.CHROM, " goodABs[0].len:", goodABs[gq_cutoffs[0]].len
  if v.FILTER != "PASS": continue
  if exclude != nil and stripChr($v.CHROM) in exclude and exclude[stripChr($v.CHROM)].find(v.start.int, v.stop.int, empty):
    continue

  if v.format.get("AD", AD) != Status.OK:
    quit "WTF!!!"

  if v.format.get("GQ", GQ) != Status.OK:
    quit "WTF GQ!!!"
  doAssert AD.len == 2 * ivcf.n_samples
  doAssert GQ.len == ivcf.n_samples
  var alts = v.format.genotypes(x).alts

  for i, kid in kids:
    var lowest_gq = 99
    var lowest_d = int32.high
    #var skip = false
    for s in [kid, kid.mom, kid.dad]:
      if GQ[s.i] < lowest_gq:
        lowest_gq = GQ[s.i]
      if AD[2*s.i] + AD[2*s.i+1] < lowest_d:
        lowest_d = AD[2*s.i] + AD[2*s.i+1]
    if lowest_d < depth_cutoff:
      continue

    #  if AD[2*s.i] + AD[2*s.i+1] < 6:
    #    skip = true
    #if skip: continue
    #stderr.write_line "loest:", $lowest_gq, " cutoffs:", $gq_cutoffs

    var ab = kid.good(AD, alts)
    if ab >= 0:
      for gq in gq_cutoffs:
        if lowest_gq < gq: continue
        goodABs[gq].add(ab)

    ab = kid.bad(AD, alts)
    if ab >= 0:
      for gq in gq_cutoffs:
        if lowest_gq < gq: continue
        badABs[gq].add(ab)

var total = goodABs[gq_cutoffs[0]].len + badABs[gq_cutoffs[0]].len
echo "#GQ\texclude\tcutoff\tfp\ttp\ttotal\tn_fps\tn_tps"

for i, (goods, bads) in [(goodABs, badABs)]:
  var ex = if opts.exclude != "": "exclude" else: "."
  for gq_cutoff in gq_cutoffs:
    var good = goods[gq_cutoff]
    var bad = bads[gq_cutoff]
    sort(good)
    sort(bad)

    var cutoff = 0'f64
    while cutoff < 0.5:
      var gi = lowerBound(good, cutoff)
      var bi = lowerBound(bad, cutoff)

      var fpr = (bad.len - bi).float64 / total.float64
      var tpr = (good.len - gi).float64 / total.float64

      # we exclude anything less than the cutoff
      # so a FP is bad above the cutoff
      echo &"{gq_cutoff}\t{ex}\t{cutoff}\t{fpr}\t{tpr}\t{total}\t{bad.len - bi}\t{good.len - gi}"
      cutoff += 0.001

