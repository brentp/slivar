import strutils
import os
import times
import algorithm
import strformat
import random

import minizip
import ./pracode
import hts

type Long* = object
  # values that are too long to be encoded get put here.
  position*: uint32
  reference*: string
  alternate*: string
  values: seq[float32]
  filter: bool

# Same as Missing in htslib
const MissingVal = float32(0x7F800001) # float32.low

type Gnotater* = ref object
  zip: Zip
  missing_value: float32 ## annotate variants absent from gnotate file with an allele frequency of -1.0
  chrom: cstring ## this tracks the current chromosome.
  encs: seq[uint64] ## encoded position,ref,alt. alts > 1 are encoding a FILTER as well.
  names: seq[string] ## this is what gets set in the INFO. e.g. 'gnomad_af'.
  values: seq[seq[float32]] ## encoded values, e.g. allele frequencies from gnomad.
  longs: seq[Long] ## pos, ref, alt, af for long variants
  chroms: seq[string] ## list of chroms available in the index.
  tmpDir: string

proc close*(g:var Gnotater) =
  g.zip.close()

template names*(g: Gnotater): seq[string] =
  g.names

proc n_fields(g: Gnotater): int =
  return g.names.len

proc cmp_long*(a, b:Long): int =
  if a.position != b.position:
    return cmp[uint32](a.position, b.position)
  if a.reference != b.reference:
    return cmp(a.reference.len, b.reference.len)
  # can't compare .af here because it screws up lower bound
  return cmp(a.alternate.len, b.alternate.len)

randomize()

proc open*(g:var Gnotater, zpath: string, tmpDir:string="/tmp", missing_val:float32= -1.0'f32): bool =
  g = Gnotater(tmpDir:tmpDir, missing_value:missing_val, names: @[], chroms: @[])
  if not open(g.zip, zpath):
    stderr.write_line &"[slivar] error opening {zpath} for annotation"
    return false

  var s = newString(1)
  doAssert g.zip.read_into("sli.var/chroms.txt", s)
  for l in s.strip().split(seps={'\n'}):
    g.chroms.add(l.strip())

  doAssert g.zip.read_into("sli.var/fields.txt", s)
  for l in s.strip().split(seps={'\n'}):
    g.names.add(l.strip())

  if getEnv("SLIVAR_QUIET") == "" and g.zip.read_into("sli.var/message.txt", s):
    stderr.write_line "[slivar] message for " & zpath & ":"
    for l in s.strip().split(seps={'\n'}):
      if l.strip() != "":
        stderr.write_line "   > " & l
  g.encs = newSeq[uint64]()
  g.values = newSeq[seq[float32]](g.n_fields)
  for i in 0..<g.values.len:
    g.values[i] = newSeq[float32](0)
  g.longs = newSeq[Long]()

  return true

proc parseLong(line: string, n_fields:int): Long {.inline.} =
  var i = -1
  for t in line.split(seps={'\t', '\n'}, maxsplit=4):
    i += 1
    if i == 0:
      result.position = parseInt(t).uint32
    elif i == 1:
      result.reference = t
    elif i == 2:
      result.alternate = t
    elif i == 3:
      result.filter = t[0] == 't'
    elif i == 4:
      let ts = t.split(seps={'|'}, maxsplit=n_fields)
      result.values.setLen(ts.len)
      for i, it in ts:
        result.values[i] = parseFloat(it)
      return
  if i < 3:
    quit "bad line:" & line


proc readLongs(g:var Gnotater, chrom: string) =

  var big = newString(16)
  doAssert g.zip.read_into(&"sli.var/{chrom}/long-alleles.txt", big)
  var i = 0
  for l in big.split(seps={'\n'}):
    if unlikely(l.len == 0): continue
    if unlikely(i == 0 and l.startswith("position")): continue
    g.longs.add(parseLong(l, g.n_fields))#.strip(chars={'\n'})))
    i += 1

proc readEncs(g:var Gnotater, chrom: string) =
  doAssert g.zip.readInto(&"sli.var/{chrom}/gnotate-variant.bin", g.encs)

proc readValues(g: var Gnotater, chrom: string, field_i: int) =
  var field_name = g.names[field_i]
  # gnotate-gnomad_num_homalt.bin
  shallow(g.values)
  doAssert g.zip.readInto(&"sli.var/{chrom}/gnotate-{field_name}.bin", g.values[field_i])

proc readValues(g: var Gnotater, chrom: string) =
  for i in 0..<g.n_fields:
    g.readValues(chrom, i)

proc sanitize_chrom*(c:string): string {.inline.} =
  if c.len == 1: return c
  result = c
  if c.len > 3 and c[0] == 'c' and c[1] == 'h' and c[2] == 'r':
    result = result[3..result.high]
  if result == "MT": result = "M"

proc load(g:var Gnotater, chrom: cstring): bool =
  g.chrom = chrom
  var chrom = sanitize_chrom($chrom)
  if chrom notin g.chroms:
    g.encs.setLen(0)
    g.values.setLen(0)
    g.longs.setLen(0)
    return false

  when defined(gnotate_times):
    var t = cpuTime()
  g.readEncs(chrom)
  when defined(gnotate_times):
    var etime = cpuTime() - t
  var t2 = cpuTime()
  g.readValues(chrom)
  when defined(gnotate_times):
    var atime = cpuTime() - t2
  t2 = cpuTime()
  g.readLongs(chrom)
  when defined(gnotate_times):
    var ltime = cpuTime() - t2
  when defined(gnotate_times):
    stderr.write_line &"len: {g.encs.len}. time to extract encs: {etime:.3f} afs: {atime:.3f} longs: {ltime:.3f} total:{cpuTime() - t:.3f}"
  for i, v in g.values:
    doAssert v.len == g.encs.len, $(i, v.len, g.encs.len)
  return true

proc values_at(g: Gnotater, i: int): seq[float32] =
  result = newSeq[float32](g.n_fields)
  for k in 0..g.values.high:
    # g.values  is shape (n_fields, n_variants)
    result[k] = g.values[k][i]

proc show*(g:var Gnotater, chrom:string, start:int, stop:int) =
  if g.chrom != chrom:
   discard g.load(chrom)

  var q = pfra(position:start.uint32, reference:"A", alternate:"A")
  var i = g.encs.lowerBound(q.encode)
  q.position = stop.uint32
  var j = max(0, g.encs.lowerBound(q.encode))
  i = max(0, i - 5)
  j = min(j + 5, g.encs.high)

  for k in i..j:
    echo g.encs[k].decode, " # ", g.encs[k], " values:", g.values_at(k)

proc annotate_missing(g:Gnotater, v:Variant): bool {.inline.} =
  if g.missing_value != MissingVal:
    var values = @[g.missing_value]
    for n in g.names:
      if v.info.set(n, values) != Status.OK:
        quit &"couldn't set info of {n} to {values[0]}"
    return true
  return false

proc contains*(g: var Gnotater, v:Variant): bool =
  ## fast-path to check presence variant in set
  if g.chrom != v.CHROM:
    discard g.load(v.CHROM)
  let alt = if likely(len(v.ALT) > 0): v.ALT[0] else: "."
  let q = if unlikely(v.REF.len + alt.len > MaxCombinedLen):
    pfra(position:v.start.uint32)
  else:
    pfra(position:v.start.uint32, reference:v.REF, alternate:alt)
  var i = g.encs.find(q)
  if i == -1: return false
  let match = g.encs[i].decode

  if match.reference.len != 0:
    # found short allele
    return true

  # should find this position in the longs
  let l = Long(position:v.start.uint32, reference:v.REF, alternate:alt)
  i = lowerBound(g.longs, l, cmp_long)
  # since these can be ordered differently, we have to check until we get to a different position or a match.
  while i < g.longs.len:
    if i > g.longs.high or g.longs[i].position != q.position:
      return g.annotate_missing(v)
    if g.longs[i].reference == l.reference and g.longs[i].alternate == l.alternate:
      break
    i += 1
  return i < g.longs.len

proc annotate*(g:var Gnotater, v:Variant): bool {.inline.} =
  ## annotate the variant INFO field with the allele frequencies and flags in g
  ## if include_missing is true. the allele frequency will be set to -1 if the variant
  ## is not found.
  if g.chrom != v.CHROM:
    discard g.load(v.CHROM)

  if g.encs.len == 0:
    return g.annotate_missing(v)

  var q:pfra
  var alt = if likely(len(v.ALT) > 0): v.ALT[0] else: "."
  if unlikely(v.REF.len + alt.len > MaxCombinedLen):
    q = pfra(position:v.start.uint32)
  else:
    q = pfra(position:v.start.uint32, reference:v.REF, alternate:alt)
  var i = g.encs.find(q)
  if i == -1:
    return g.annotate_missing(v)

  var values = g.values_at(i)

  var match = g.encs[i].decode
  var filtered = match.filter

  if match.reference.len == 0:
    # should find this position in the longs
    var l = Long(position:v.start.uint32, reference:v.REF, alternate:alt)
    var i = lowerBound(g.longs, l, cmp_long)
    # since these can be ordered differently, we have to check until we get to a different position or a match.
    while i < g.longs.len:
      if i > g.longs.high or g.longs[i].position != q.position:
        return g.annotate_missing(v)
      if g.longs[i].reference == l.reference and g.longs[i].alternate == l.alternate:
        break
      i += 1

    if unlikely(i > g.longs.high):
      return g.annotate_missing(v)

    values = g.longs[i].values
    filtered = g.longs[i].filter
    doAssert values.len != 0

  for i, value in values:
    var ivalues = @[value]
    if v.info.set(g.names[i], ivalues) != Status.OK:
      echo v.tostring()
      quit &"couldn't set info of {g.names[i]} to {values[0]}"

    if filtered:
      if v.info.set(g.names[i] & "_filter", true) != Status.OK:
        quit &"couldn't set info flag {g.names[i] & \"_filter\"}"
  return true

proc update_header*(g:Gnotater, ivcf:VCF) =
  for n in g.names:
    doAssert ivcf.header.add_info(n, "1", "Float", "field from from gnotate VCF") == Status.OK
    doAssert ivcf.header.add_info(n & "_filter", "0", "Flag", "non-passing flag in gnotate source VCF") == Status.OK

when isMainModule:

  #[
  var zip_path = paramStr(1)
  var loc = paramStr(2)
  var g:Gnotater
  if not g.open(zip_path, tmpDir=getTempDir(), missing_val= -1.0'f32):
    quit "[slivar] error opening zip. check path and contents"
  var se = loc.split(":")[1].split("-")
  g.show(loc.split(":")[0], parseInt(se[0])-1, parseInt(se[1]))
  ]#



  import times

  var vcf_path = paramStr(1)

  var ivcf:VCF
  if not open(ivcf, vcf_path, threads=2):
    quit "couldn't open path" & vcf_path

  var ovcf:VCF
  if not open(ovcf, "_t.bcf", mode="wb", threads=3):
    quit "couldn't open output file"

  var zip_path = paramStr(2)
  var original_field = paramStr(3)

  var g:Gnotater
  if not g.open(zip_path, tmpDir="/tmp", missing_val= -1.0'f32):
    quit "[slivar] error opening zip. check path and contents"
  g.update_header(ivcf)

  ovcf.copy_header(ivcf.header)

  doAssert ovcf.write_header

  var nv = 0
  var t0 = cpuTime()
  for v in ivcf:
    if v.rid > 4: break
    #if v.ALT[0] != "C": continue
    discard g.annotate(v)
    doAssert ovcf.write_variant(v)
    nv += 1

  var vps = nv.float64 / (cpuTime() - t0)
  echo &"annotated {nv} variants in {cpuTime() - t0:.0f} seconds: {vps.int} variants/sec"
  ovcf.close()
  ivcf.close()

  ## now open the output file and check that values match.
  if not open(ivcf, "_t.bcf", threads=3):
    quit "couldn't open file to check"

  var n = 0
  var nmiss = 0
  var max_diff = 0'f32
  for v in ivcf:
    if v.rid > 4: break
    #if v.ALT[0] != "C": continue
    for n in g.names[0..<1]:
      var floats = newSeq[float32]()
      var ints = newSeq[int32]()
      if v.info.get(n, floats) != Status.OK:
        quit "FAIL 0. not annotated:" & v.tostring()
      var val = floats[0]
      var filt = v.FILTER

      var oval = -1'f32
      var st = v.info.get(original_field, floats)
      if st == Status.UnexpectedType:
        if v.info.get(original_field, ints) == Status.OK:
          oval = ints[0].float32
      elif st == Status.NotFound:
        if filt == "PASS":
          if nmiss < 10:
            echo "not annotated:" & v.tostring()[0..150]
          nmiss += 1
      else:
        oval = floats[0]

      if abs(val - oval) > max_diff:
        max_diff = abs(val - oval)
        if max_diff < 0.00001:
          echo "current maximum difference:", $max_diff
        else:
          echo "FAIL: differing values for:" & v.tostring()
          quit "got:" & $val & " expected:" & $oval

      if not v.info.has_flag(n & "_filter") and filt != "PASS":
          quit "FAIL: no filter found for " & v.tostring()
      break

    n += 1
  echo "PASS:", $n, " variants tested with max difference:", max_diff
