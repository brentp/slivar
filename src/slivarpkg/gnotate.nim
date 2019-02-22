import strutils
import streams
import os
import times
import algorithm
import strformat

import zip/zipfiles
#import minizip
#type ZipArchive = Zip

import ./pracode
import hts

type packed = object
  value {.bitsize: 27.}: float32
  flag {.bitsize: 5.}: float32

type Long* = object
  # values that are too long to be encoded get put here.
  position*: uint32
  reference*: string
  alternate*: string
  value: float32

# Same as Missing in htslib
const MissingVal = float32(0x7F800001) # float32.low

type Gnotater* = ref object
  zip: ZipArchive
  name*: string ## this is what gets set in the INFO. e.g. 'gnomad_af'.
  missing_value: float32 ## annotate variants absent from gnomad with an allele frequency of -1.0
  chrom: cstring ## this tracks the current chromosome.
  encs: seq[uint64] ## encoded position,ref,alt. alts > 1 are encoding a FILTER as well.
  values: seq[float32] ## encoded values, e.g. allele frequencies from gnomad.
  longs: seq[Long] ## pos, ref, alt, af for long variants
  filters: seq[string] ## ordered FILTER fields from gnomad
  chroms: seq[string] ## list of chroms available in the index.
  tmpDir: string

proc close*(g:var Gnotater) =
  g.zip.close()

proc cmp_long*(a, b:Long): int =
  if a.position != b.position:
    return cmp[uint32](a.position, b.position)
  if a.reference != b.reference:
    return cmp(a.reference.len, b.reference.len)
  # can't compare .af here because it screws up lower bound
  return cmp(a.alternate.len, b.alternate.len)

proc open*(g:var Gnotater, path: string, name:string="gnomad_af", tmpDir:string="/tmp", missing_val:float32= -1.0'f32): bool =
  g = Gnotater(name:name, tmpDir:tmpDir, missing_value:missing_val)
  if not open(g.zip, path):
    return false
  var path = g.tmpDir / "filters.txt"
  g.zip.extract_file("sli.var/filters.txt", path)
  for l in path.lines:
    g.filters.add(l.strip)
  removeFile(path)

  path = g.tmpDir / "chroms.txt"
  g.zip.extract_file("sli.var/chroms.txt", path)
  for l in path.lines:
    g.chroms.add(l.strip)
  removeFile(path)

  g.encs = newSeq[uint64]()
  g.values = newSeq[float32]()
  g.longs = newSeq[Long]()

  return true

proc parseLong(line: string): Long {.inline.} =
  var i = -1
  for t in line.split(seps={'\t', '\n'}, maxsplit=3):
    i += 1
    if i == 0:
      result.position = parseInt(t).uint32
    elif i == 1:
      result.reference = t
    elif i == 2:
      result.alternate = t
    elif i == 3:
      result.value = parseFloat(t)
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
  for l in big.split(seps={'\n'}):
    if unlikely(l.len == 0): continue
    g.longs.add(parseLong(l))#.strip(chars={'\n'})))
  st.close()

proc readEncs(g:Gnotater, chrom: string) =
  var st = g.zip.getStream(&"sli.var/{chrom}/gnotate-variants.bin")
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

proc readAfs(g: var Gnotater, chrom: string) =
  var st = g.zip.getStream(&"sli.var/{chrom}/gnotate-values.bin")
  var chunk = 21660531
  shallow(g.values)
  if g.values.len > chunk:
    g.values.setLen(chunk)
  else:
    g.values = newSeqUninitialized[float32](chunk)
  while true:
    let bytesRead = st.readData(g.values[g.values.len - chunk].addr, chunk * float32.sizeof)
    if bytesRead < chunk * float32.sizeof:
      g.values.setLen(g.values.len - chunk + int(bytesRead / float32.sizeof))
      break
    g.values.setLen(g.values.len + chunk)

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
  g.readAfs(chrom)
  var atime = cpuTime() - t2
  t2 = cpuTime()
  g.readLongs(chrom)
  var ltime = cpuTime() - t2
  echo &"len: {g.encs.len}. time to extract encs: {etime:.3f} afs: {atime:.3f} longs: {ltime:.3f} total:{cpuTime() - t:.3f}"
  doAssert g.values.len == g.encs.len
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
    if v.info.set(g.name, values) != Status.OK:
      quit &"couldn't set info of {g.name} to {values[0]}"
    return true
  return false

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
  if unlikely(v.REF.len + alt.len > 14):
    q = pfra(position:v.start.uint32)
  else:
    q = pfra(position:v.start.uint32, reference:v.REF, alternate:alt)
  var i = g.encs.find(q)
  if i == -1:
    return g.annotate_missing(v)

  var value = g.values[i]
  var filter:string

  var match = g.encs[i].decode
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

    value = g.longs[i].value

  if value > 1:
    # the value also stores the ith flag.
    # 0 is pass and 1 is unused so if the value is > 1
    # e.g. 2.001 means the 2nd flag and an AF of 0.001
    var ifilt = value.int
    filter = g.filters[ifilt]
    value -= ifilt.float32
    if value > 1:
      echo ifilt, " ", value, " filter:", filter, " ifilt:", ifilt

  var values = @[value]
  if v.info.set(g.name, values) != Status.OK:
    echo "variant:"
    echo v.tostring()
    quit &"couldn't set info of {g.name} to {values[0]}"

  if filter.len > 0 and filter != "PASS":
    if v.info.set(g.name & "_filter", filter) != Status.OK:
      quit &"couldn't set info of {g.name & \"_filter\"} to {filter}"
  return true

proc update_header*(g:Gnotater, ivcf:VCF) =
  doAssert ivcf.header.add_info(g.name, "1", "Float", "field from from gnomad VCF") == Status.OK
  doAssert ivcf.header.add_info(g.name & "_filter", "1", "String", "flag from gnomad VCF") == Status.OK

when isMainModule:
  import times

  var vcf_path = commandLineParams()[0]

  var ivcf:VCF
  if not open(ivcf, vcf_path, threads=2):
    quit "couldn't open path" & vcf_path

  var ovcf:VCF
  if not open(ovcf, "t.bcf", mode="w", threads=2):
    quit "couldn't open output file"

  var g:Gnotater
  var zip_path = commandLineParams()[1]
  doAssert g.open(zip_path, name="gnomad_af", tmpDir="/tmp", missing_val= -1.0'f32)
  g.update_header(ivcf)

  ovcf.copy_header(ivcf.header)

  doAssert ovcf.write_header

  var nv = 0
  var t0 = cpuTime()
  for v in ivcf:
    discard g.annotate(v) #, v.tostring()[0..<150]
    doAssert ovcf.write_variant(v)
    nv += 1

  var vps = nv.float64 / (cpuTime() - t0)
  echo &"annotated {nv} variants in {cpuTime() - t0:.0f} seconds: {vps.int} variants/sec"
  ovcf.close()
  ivcf.close()

  ## now open the output file and check that values match.
  var floats = newSeq[float32]()
  var strings = ""
  if not open(ivcf, "t.bcf", threads=3):
    quit "couldn't open file to check"

  var n = 0
  var nmiss = 0
  var max_diff = 0'f32
  for v in ivcf:
    if v.info.get(g.name, floats) != Status.OK:
      quit "FAIL. not annotated:" & v.tostring()
    var aaf = floats[0]
    var filt = v.FILTER

    var oaf = 0'f32
    if v.info.get("AF_popmax", floats) != Status.OK:
      if filt == "PASS":
        if nmiss < 10:
          echo "missing AF_popmax not annotated:" & v.tostring()[0..150]
        nmiss += 1
    else:
      oaf = floats[0]

    if abs(aaf - oaf) > max_diff:
        max_diff = abs(aaf - oaf)
        if max_diff < 0.001:
          echo "current maximum difference:", $max_diff
        else:
          echo "FAIL: differing values for:" & v.tostring()
          quit "got:" & $aaf & " expected:" & $oaf

    if v.info.get("gnomad_af_filter", strings) != Status.OK and filt != "PASS":
        quit "FAIL: no filter found for " & v.tostring()

    if filt == "PASS": filt = ""
    if filt != strings:
      echo  "FAIL: differing filters found for " & v.tostring()
      quit "got:" & strings & " expected:" & filt
    n += 1
  echo "PASS:", $n, " variants tested with max difference:", max_diff
