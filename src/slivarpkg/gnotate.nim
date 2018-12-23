import strutils
import streams
import os
import algorithm
import strformat

import zip/zipfiles
#import minizip
#type ZipArchive = Zip

import ./pracode
import hts

type Long* = object
  # Too long to be encoded
  position*: uint32
  reference*: string
  alternate*: string
  af: float32

type Gnotater* = ref object
  zip: ZipArchive
  name: string
  chrom: cstring
  encs: seq[uint64]
  afs: seq[float32]
  longs: seq[Long]
  filters: seq[string]
  chroms: seq[string]
  tmpDir: string

proc close(g:var Gnotater) =
  g.zip.close()

proc cmp_long*(a, b:Long): int =
  if a.position != b.position:
    return cmp[uint32](a.position, b.position)
  if a.reference != b.reference:
    return cmp(a.reference.len, b.reference.len)
  # can't compare .af here because it screws up lower bound
  return cmp(a.alternate.len, b.alternate.len)

proc open*(g:var Gnotater, path: string, name:string="gnomad_af", tmpDir="/tmp"): bool =
  g.tmpDir = tmpDir
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
  g.afs = newSeq[float32]()
  g.longs = newSeq[Long]()


  return true

proc parseLong(line: string): Long {.inline.} =
  var t = line.strip().split(seps={'\t'})
  result.position = parseInt(t[0]).uint32
  result.reference = t[1]
  result.alternate = t[2]
  result.af = parseFloat(t[3])

proc readLongs(g:var Gnotater, chrom: string) =
  var path = g.tmpDir / "long-alleles.AF_popmax.txt"
  g.zip.extract_file(&"sli.var/{chrom}/long-alleles.AF_popmax.txt", path)
  #var z = g.zip
  #var path = z.extract_file(&"sli.var/{chrom}/long-alleles.AF_popmax.txt", g.tmpDir)

  g.longs.setLen(0)
  for l in path.lines:
    g.longs.add(parseLong(l))
  removeFile(path)

proc readEncs(g:Gnotater, chrom: string) =
  var path = g.tmpDir / "vk.bin"
  g.zip.extract_file(&"sli.var/{chrom}/vk.bin", path)
  #var z = g.zip
  #var path = z.extract_file(&"sli.var/{chrom}/vk.bin", g.tmpDir)
  var size = path.getFileSize.int
  g.encs = newSeq[uint64](int(size / uint64(0).sizeof))
  var fs = newFileStream(path, fmRead)
  if fs == nil:
    quit "couldn't open pop max file from:" & path
  doAssert fs.readData(g.encs[0].addr.pointer, size) == size
  fs.close()
  removeFile(path)

proc readAfs(g: var Gnotater, chrom: string) =
  var path = g.tmpDir / "vk-AF.bin"
  g.zip.extract_file(&"sli.var/{chrom}/vk-AF_popmax.bin", path)
  var size = path.getFileSize.int
  var L = int(size / float32(0).sizeof)
  g.afs = newSeqUninitialized[float32](L)

  var fs = newFileStream(path, fmRead)
  if fs == nil:
    quit "couldn't open pop max file from:" & path
  doAssert fs.readData(g.afs[0].addr.pointer, size) == size
  fs.close()
  removeFile(path)

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
    g.afs.setLen(0)
    g.longs.setLen(0)
    return false

  g.readEncs(chrom)
  g.readAfs(chrom)
  g.readLongs(chrom)
  doAssert g.afs.len == g.encs.len
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

proc annotate*(g:var Gnotater, v:Variant): bool {.inline.} =
  if len(v.ALT) > 1:
    echo "only annotating a single allele for the multiallelic variants"
  if g.chrom != v.CHROM:
    discard g.load(v.CHROM)

  if g.encs.len == 0: return false

  var q:pfra
  if unlikely(v.REF.len + v.ALT[0].len > 14):
    q = pfra(position:v.start.uint32)
  else:
    q = pfra(position:v.start.uint32, reference:v.REF, alternate:v.ALT[0])
  var i = g.encs.find(q)
  if i == -1: return false
  var value = g.afs[i]
  var filter:string

  var match = g.encs[i].decode
  if match.reference.len == 0:
    # should find this position in the longs
    var l = Long(position:v.start.uint32, reference:v.REF, alternate:v.ALT[0])
    var i = lowerBound(g.longs, l, cmp_long)
    # since these can be ordered differently, we have to check until we get to a different position or a match.
    while i < g.longs.high:
      if i > g.longs.high or g.longs[i].position != q.position: return false
      if g.longs[i].reference == l.reference and g.longs[i].alternate == l.alternate:
        break
      i += 1

    value = g.longs[i].af

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

when isMainModule:
  import times

  proc update_header(g:Gnotater, ivcf:VCF) =
    doAssert ivcf.header.add_info(g.name, "1", "Float", "field from from gnomad VCF") == Status.OK
    doAssert ivcf.header.add_info(g.name & "_filter", "1", "String", "flag from gnomad VCF") == Status.OK

  var vcf_path = commandLineParams()[0]

  var ivcf:VCF
  if not open(ivcf, vcf_path, threads=2):
    quit "couldn't open path" & vcf_path

  var ovcf:VCF
  if not open(ovcf, "t.bcf", mode="w", threads=2):
    quit "couldn't open output file"

  var g = Gnotater(name: "gnomad_af")
  g.update_header(ivcf)

  ovcf.copy_header(ivcf.header)

  doAssert ovcf.write_header

  var zip_path = commandLineParams()[1]
  doAssert g.open(zip_path)

  var nv = 0
  var t0 = cpuTime()
  for v in ivcf:
    doAssert g.annotate(v), v.tostring()[0..<150]
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
