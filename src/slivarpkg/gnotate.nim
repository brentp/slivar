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
    return cmp(a.reference, b.reference)
  if a.alternate != b.alternate:
    return cmp(a.alternate, b.alternate)
  return cmp(b.af, a.af)

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
  result.position = parseInt(t[1]).uint32
  result.reference = t[2]
  result.alternate = t[3]
  result.af = parseFloat(t[4])

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
  #var path = g.zip.extract_file(&"sli.var/{chrom}/vk-AF_popmax.bin", g.tmpDir)
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
  ## TODO: read the zip file for this chrom and fill afs, encs
  doAssert g.afs.len == g.encs.len
  return true

proc annotate*(g:var Gnotater, v:Variant): bool =
  if len(v.ALT) > 1:
    echo "only annotating a single allele for the multiallelic variants"
  if g.chrom != v.CHROM:
    discard g.load(v.CHROM)

  var q = pfra(position:v.start.uint32, reference:v.REF, alternate:v.ALT[0])
  var i = g.encs.find(q)
  var value = g.afs[i]
  var filter:string

  if i == -1: return false
  var match = g.encs[i].decode
  if match.reference.len == 0:
    # should find this position in the longs
    var l = Long(position:v.start.uint32, reference:v.REF, alternate:v.ALT[0])
    var i = lowerBound(g.longs, l, cmp_long)
    if i >= g.longs.len or g.longs[i].position != q.position or g.longs[i].reference != q.reference or g.longs[i].alternate != q.alternate: return false
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

when isMainModule:
  var vcf_path = commandLineParams()[0]

  var ivcf:VCF
  if not open(ivcf, vcf_path, threads=2):
    quit "couldn't open path" & vcf_path

  var ovcf:VCF
  if not open(ovcf, "t.bcf", mode="w", threads=2):
    quit "couldn't open output file"

  doAssert ivcf.header.add_info("gnomad_af", "1", "Float", "popmax from gnomad") == Status.OK
  doAssert ivcf.header.add_info("gnomad_af_filter", "1", "String", "flag from gnomad") == Status.OK
  ovcf.copy_header(ivcf.header)

  doAssert ovcf.write_header

  var zip_path = commandLineParams()[1]
  var g = Gnotater()
  doAssert g.open(zip_path)
  g.name = "gnomad_af"

  for v in ivcf:
    discard g.annotate(v)
    doAssert ovcf.write_variant(v)

  ovcf.close()
