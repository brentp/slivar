import strutils
import streams
import os
import times
import algorithm
import strformat
import random

import zip/zipfiles
#import minizip
#type ZipArchive = Zip

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
  zip: ZipArchive
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
  g = Gnotater(tmpDir:tmpDir, missing_value:missing_val)
  if not open(g.zip, zpath):
    stderr.write_line &"[slivar] error opening {zpath} for annotation"
    return false

  var r = $random(int.high) & $random(int.high)
  var path = g.tmpDir / &"chroms{r}.txt"
  g.zip.extract_file("sli.var/chroms.txt", path)
  for l in path.lines:
    g.chroms.add(l.strip)
  removeFile(path)

  path = g.tmpDir / &"fields{r}.txt"
  g.zip.extract_file("sli.var/fields.txt", path)
  for l in path.lines:
    g.names.add(l.strip())
  removeFile(path)

  var hasMsg = false
  for p in g.zip.walkFiles:
    if p.endsWith("message.txt"):
      hasMsg = true
      break

  if hasMsg:
    path = g.tmpDir / &"message{r}.txt"
    g.zip.extract_file("sli.var/message.txt", path)
    stderr.write_line "[slivar] message for " & zpath & ":"
    for l in path.lines:
      if l.strip() != "":
        stderr.write_line "   > " & l
    removeFile(path)

  g.encs = newSeq[uint64]()
  g.values = newSeq[seq[float32]](g.n_fields)
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
      var ts = t.split(seps={'|'}, maxsplit=n_fields)
      result.values.setLen(ts.len)
      for i, it in ts:
        result.values[i] = parseFloat(it)
      return
  if i < 3:
    quit "bad line:" & line

proc readLongs(g:var Gnotater, chrom: string) =
  let st = g.zip.getStream(&"sli.var/{chrom}/long-alleles.txt")
  g.longs.setLen(0)
  # reading into memory here is the fastest way I could find. I think the calls
  # to streams are not inlined and therefore can be expensive. reading this
  # file is still slower than reading the encoded binary data.
  let big = st.readAll
  var i = 0
  for l in big.split(seps={'\n'}):
    if unlikely(l.len == 0): continue
    if unlikely(i == 0 and l.startswith("position")): continue
    g.longs.add(parseLong(l, g.n_fields))#.strip(chars={'\n'})))
    i += 1
  st.close()

proc readEncs(g:Gnotater, chrom: string) =
  var st = g.zip.getStream(&"sli.var/{chrom}/gnotate-variant.bin")
  var chunk = 21660531 # this is the exact number of elements in chr1 so it makes this the fastest
  if g.encs.len > chunk:
    g.encs.setLen(chunk)
  else:
    g.encs = newSeqUninitialized[uint64](chunk)
  while true:
    let bytesRead = st.readData(g.encs[g.encs.len - chunk].addr, chunk * uint64.sizeof)
    if bytesRead < chunk * uint64.sizeof:
      g.encs.setLen(g.encs.len - chunk + int(bytesRead / uint64.sizeof))
      break
    g.encs.setLen(g.encs.len + chunk)
  st.close()

proc readValues(g: var Gnotater, chrom: string, field_i: int) =
  var field_name = g.names[field_i]
  # gnotate-gnomad_num_homalt.bin
  var st = g.zip.getStream(&"sli.var/{chrom}/gnotate-{field_name}.bin")
  var chunk = 216 #60531
  shallow(g.values)
  shallow(g.values[field_i])
  if g.values[field_i].len >= chunk:
    g.values[field_i].setLen(chunk)
  else:
    g.values[field_i] = newSeqUninitialized[float32](chunk)
  shallow(g.values[field_i])
  shallow(g.values)
  #var vs = g.values[field_i]
  #echo "len:", vs.len, " chunk:", chunk

  while true:
    let bytesRead = st.readData(g.values[field_i][g.values[field_i].len - chunk].addr, chunk * float32.sizeof)
    if bytesRead < chunk * float32.sizeof:
      g.values[field_i].setLen(g.values[field_i].len - chunk + int(bytesRead / float32.sizeof))
      break
    g.values[field_i].setLen(g.values[field_i].len + chunk)

proc readValues(g: var Gnotater, chrom: string) =
  for i in 0..<g.n_fields:
    g.readValues(chrom, i)

proc sanitize_chrom(c:string): string {.inline.} =
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

  var t = cpuTime()
  g.readEncs(chrom)
  var etime = cpuTime() - t
  var t2 = cpuTime()
  g.readValues(chrom)
  var atime = cpuTime() - t2
  t2 = cpuTime()
  g.readLongs(chrom)
  var ltime = cpuTime() - t2
  when defined(gnotate_times):
    stderr.write_line &"len: {g.encs.len}. time to extract encs: {etime:.3f} afs: {atime:.3f} longs: {ltime:.3f} total:{cpuTime() - t:.3f}"
  for v in g.values:
    doAssert v.len == g.encs.len
  return true

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
    echo g.encs[k].decode, " # ", g.encs[k]

proc annotate_missing(g:Gnotater, v:Variant): bool {.inline.} =
  if g.missing_value != MissingVal:
    var values = @[g.missing_value]
    for n in g.names:
      if v.info.set(n, values) != Status.OK:
        quit &"couldn't set info of {n} to {values[0]}"
    return true
  return false

proc values_at(g: Gnotater, i: int): seq[float32] =
  result.setLen(g.n_fields)
  for k in 0..g.values.high:
    # g.values  is shape (n_fields, n_variants)
    result[k] = g.values[k][i]

proc annotate*(g:var Gnotater, v:Variant): bool {.inline.} =
  ## annotate the variant INFO field with the allele frequencies and flags in g
  ## if include_missing is true. the allele frequency will be set to -1 if the variant
  ## is not found.
  #if len(v.ALT) > 1:
    #echo "only annotating a single allele for the multiallelic variants"
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
    while i < g.longs.high:
      if i > g.longs.high or g.longs[i].position != q.position:
        return g.annotate_missing(v)
      if g.longs[i].reference == l.reference and g.longs[i].alternate == l.alternate:
        break
      i += 1

    values = g.longs[i].values
    filtered = g.longs[i].filter

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
  import times

  var vcf_path = paramStr(1)

  var ivcf:VCF
  if not open(ivcf, vcf_path, threads=2):
    quit "couldn't open path" & vcf_path

  var ovcf:VCF
  if not open(ovcf, "t.bcf", mode="wb", threads=3):
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
  if not open(ivcf, "t.bcf", threads=3):
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
      if v.info.get(original_field, floats) != Status.OK:
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
