import hts/vcf
import strformat
import slivarpkg/pedfile
import strutils
import math
import algorithm

import os

var samples = parse_ped(paramStr(1))
var ivcf:VCF
if not open(ivcf, paramStr(2)):
  quit "couldn't open vcf"

var gq_cutoff = parseInt(paramStr(3))

samples = samples.match(ivcf)

var kids = newSeq[Sample]()

for s in samples:
  if s.dad != nil and s.mom != nil:
    kids.add(s)

proc ab(a:array[2, int32]): float64 {.inline.} =
  result = a[1].float64 / max(1, sum(a)).float64
  if result > 0.5: result = 1 - result

proc good(kid:Sample, AD: seq[int32], alts: seq[int8]): float64 {.inline.} =
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

var goodABs = newSeq[float64](kids.len)
var badABs = newSeq[float64](kids.len)

var AD: seq[int32]
var GQ: seq[int32]
var x: seq[int32]


for v in ivcf:
  if v.format.get("AD", AD) != Status.OK:
    quit "WTF!!!"

  if v.format.get("GQ", GQ) != Status.OK:
    quit "WTF GQ!!!"
  doAssert AD.len == 2 * ivcf.n_samples
  doAssert GQ.len == ivcf.n_samples
  var alts = v.format.genotypes(x).alts

  for i, kid in kids:
    var skip = false
    for s in [kid, kid.mom, kid.dad]:
      if GQ[s.i] < gq_cutoff:
        skip = true
        break
      #if AD[2*s.i] + AD[2*s.i+1] < 5:
      #  skip = true
      #  break
    if skip: continue

    var ab = kid.good(AD, alts)
    if ab >= 0:
      goodABs.add(ab)
    ab = kid.bad(AD, alts)
    if ab >= 0:
      badABs.add(ab)

sort(goodABs)
sort(badABs)
var cutoff = 0'f64
echo "#cutoff\tfp\ttp\ttotal"
while cutoff < 0.5:
  var gi = lowerBound(goodABs, cutoff)
  var bi = lowerBound(badABs, cutoff)

  var fpr = (badAbs.len - bi).float64 / (badAbs.len + goodAbs.len).float64
  var tpr = (goodAbs.len - gi).float64 / (badAbs.len + goodAbs.len).float64

  # we exclude anything less than the cutoff
  # so a FP is bad above the cutoff
  echo &"{cutoff}\t{fpr}\t{tpr}\t{goodAbs.len + badABs.len}"
  cutoff += 0.001
