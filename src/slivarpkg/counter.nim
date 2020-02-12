# this just counts the occurrence of sample, expr to output a final table at
# the of a slivar run.
import tables
import pedfile
import ./evaluator
import strutils

type Counter* = object
  samples: seq[string]
  exprs: seq[string]
  sampleIndexes: TableRef[string, int]
  exprIndexes: TableRef[string, int]
  counts: seq[seq[int]]

proc initCounter*(ev:Evaluator): Counter =
  result.sampleIndexes = newTable[string, int]()
  result.exprIndexes = newTable[string, int]()
  for i, s in ev.samples:
    result.sampleIndexes[s.ped_sample.id] = i
    result.samples.add(s.ped_sample.id)
  for i, e in ev.trio_expressions:
    result.exprIndexes[e.name] = result.exprIndexes.len
    result.exprs.add(e.name)
  for i, e in ev.group_expressions:
    result.exprIndexes[e.name] = result.exprIndexes.len
    result.exprs.add(e.name)
  for i, e in ev.family_expressions:
    result.exprIndexes[e.name] = result.exprIndexes.len
    result.exprs.add(e.name)

  result.counts = newSeq[seq[int]](ev.samples.len)
  for i in 0..result.counts.high:
    result.counts[i] = newSeq[int](result.exprIndexes.len)

proc initCounter*(kids: seq[Sample]): Counter =
  # for compound hets
  result.samples = newSeq[string](kids.len)
  result.sampleIndexes = newTable[string, int]()
  result.exprs = @["compound-het"]
  for i, k in kids:
    result.samples[i] = k.id
    result.sampleIndexes[k.id] = i
  result.exprIndexes = newTable[string, int]()
  result.exprIndexes["compound-het"] = 0
  result.counts = newSeq[seq[int]](kids.len)
  for i in 0..result.counts.high:
    result.counts[i] = newSeq[int](result.exprIndexes.len)


proc inc*(c:var Counter, samples:seq[string], e:string) {.inline.} =
  let jexpr = c.exprIndexes[e]
  for s in samples:
    let isample = c.sampleIndexes[s]
    c.counts[isample][jexpr].inc

proc `$`*(c:Counter): string =
  var header = "sample\t" & join(c.exprs, "\t")
  var res = newSeq[string]()
  for i, s in c.samples:
    var line:string = s
    for j, e in c.exprs:
      line &= '\t' & $c.counts[i][j]
    res.add(line)
  return header & '\n' & join(res, "\n")

proc tostring*(c:Counter, samples: seq[string]): string =
  var header = "sample\t" & join(c.exprs, "\t")
  var res = newSeq[string]()
  for i, s in c.samples:
    if samples.len > 0 and s notin samples: continue
    var line = newStringOfCap(32)
    line.add(s)
    for j, e in c.exprs:
      line &= '\t' & $c.counts[i][j]
    res.add(line)
  return header & '\n' & join(res, "\n") & '\n'

proc tostring*(c:Counter, samples: seq[pedfile.Sample]): string =

  var s = newSeq[string](samples.len)
  for i, smp in samples:
    s[i] = smp.id
  return c.tostring(s)
