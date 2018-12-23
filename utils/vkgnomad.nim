import hts/vcf
import times
import algorithm
import strformat
#import zip/zipfiles
import minizip
import slivarpkg/pracode
import docopt
import os
import strutils
import streams

let doc = """
Usage: vkgnomad [--prefix=<prefix> options] <vcfs>...

Options:

  --prefix <string>    prefix for output [default: vksm]
  --field <string>     field to use for value [default: AF_popmax]

Arguments:

  <vcfs>...    paths like: /path/to/gnomad.exomes.r2.1.sites.vcf.bgz /other/to/gnomad.genomes.r2.1.sites.vcf.bgz

"""

var args = docopt(doc)
var prefix = $args["--prefix"]
if prefix[prefix.high] == '/':
  prefix &= "vkgnomad"
if prefix[prefix.high] != '.':
  prefix &= "."

proc cleanup() {.noconv.} =
  removeDir(prefix)
addQuitProc(cleanup)

var vcf_paths = @(args["<vcfs>"])
var field = $args["--field"]


# things that are too long to be encoded.
type PosValue = tuple[chrom: string, position:pfra, value:float32]

proc write_to(positions:var seq[PosValue], fname:string) =
  # write the positions to file after sorting
  proc icmp_position(a: PosValue, b:PosValue): int =
    if a.chrom != b.chrom:
      return cmp(a.chrom, b.chrom)
    result = cmp_pfra(a.position, b.position)

    if result == 0:
      result = cmp(b.value, a.value)

  positions.sort(icmp_position)
  var last = pfra()
  var fh: File
  if not open(fh, fname, fmWrite):
    quit "couldn't open:" & fname

  var chrom = positions[0].chrom
  for pv in positions:
    if pv.position == last: continue
    var
      p = pv.position
      v = pv.value
    doAssert chrom == pv.chrom, "expecting only a single chromosome in call to write_to"
    last = p
    fh.write(&"{p.position}\t{p.reference}\t{p.alternate}\t{v}\n")
  fh.close()

var population_vcf:VCF

type evalue = tuple[encoded:uint64, value:float32]

var longs_by_rid = newSeqOfCap[seq[PosValue]](1000)
var kvs_by_rid = newSeqOfCap[seq[evalue]](1000)

shallow(kvs_by_rid)
shallow(longs_by_rid)

var longs: seq[PosValue]
var kvs: seq[evalue]

var filters = @["PASS", "UNUSED"]

var floats = newSeq[float32](1)
var ints = newSeq[float32](1)
echo prefix & "zip"

var ridToChrom = newSeqOfCap[string](100000)

var last_rid = -1
for i in 0..<vcf_paths.len:
  if not open(population_vcf, vcf_paths[i], threads=3):
    quit "couldn't open:" & vcf_paths[i]

  #var last_rid = -1'i32
  for v in population_vcf:
      if len(v.ALT) > 1:
        quit "input should be decomposed and normalized"
      if v.rid != last_rid:
        if last_rid != -1:
          longs_by_rid[last_rid] = longs
          kvs_by_rid[last_rid] = kvs
        last_rid = v.rid
        if last_rid >= longs_by_rid.len:
          longs_by_rid.setLen(last_rid + 1)
          kvs_by_rid.setLen(last_rid + 1)
          ridToChrom.setLen(last_rid + 1)
          longs_by_rid[last_rid] = newSeqOfCap[PosValue](32768)
          kvs_by_rid[last_rid] = newSeqOfCap[evalue](32768)
          ridToChrom[last_rid] = $v.CHROM
        else:
          doAssert ridToChrom[last_rid] == $v.CHROM

        longs = longs_by_rid[last_rid]
        kvs = kvs_by_rid[last_rid]

      var fil = v.FILTER
      var fidx = filters.find(fil)
      if fidx == -1:
        filters.add(fil)
        fidx = filters.high

      var e = encode(uint32(v.start), v.REF, v.ALT[0])

      if v.info.get(field, floats) != Status.OK:
        if fidx == 0 and v.info.get("AF", floats) == Status.OK and floats[0] > 0.01 and v.info.get("AN", ints) == Status.OK and ints[0] > 2000:
          quit "got wierd't get field for:" & v.tostring()
        floats = @[0'f32]

      # need to differentiate from 1 because we use the int to store flags.
      floats[0] = min(floats[0], 1-(5e-7))

      var val = floats[0] + fidx.float32
      if v.REF.len + v.ALT[0].len > 14:
        var p = e.decode()
        doAssert p.position == v.start.uint32
        p.reference = v.REF
        p.alternate = v.ALT[0]
        longs.add(($v.CHROM, p, val))
      kvs.add((e, val))
      if kvs.len mod 500_000 == 0:
        stderr.write_line &"{kvs.len} variants completed. at: {v.CHROM}:{v.start+1}. non-exact: {longs.len} in {vcf_paths[i]}"
  longs_by_rid[last_rid] = longs
  kvs_by_rid[last_rid] = kvs
  stderr.write_line &"{kvs_by_rid[last_rid].len} variants completed. non-exact: {longs.len} for {vcf_paths[i]}"

  population_vcf.close()


#var zip: ZipArchive
var zip: Zip
if not open(zip, prefix & "zip", fmWrite):
  quit "could not open zip file"

var fchrom:File
if not open(fchrom, prefix & "chroms.txt", fmWrite):
  quit "could not open chroms file"

for rid in 0..kvs_by_rid.high:
  var chrom = ridToChrom[rid]
  echo chrom
  if chrom == "": continue
  if chrom.startsWith("chr"): chrom = chrom[3..chrom.high]
  if chrom == "MT": chrom = "M"

  var kvs = kvs_by_rid[rid]
  var longs = longs_by_rid[rid]
  if kvs.len == 0 and longs.len == 0: continue
  stderr.write_line &"sorting and writing... {kvs.len} variants completed. non-exact: {longs.len} for chromosome: {chrom}"
  fchrom.write(chrom & '\n')

  longs.write_to(prefix & &"long-alleles.{field}.txt")

  kvs.sort(proc (a:evalue, b:evalue): int =
    result = cmp[uint64](a.encoded, b.encoded)
    if result == 0:
      # on ties, take the largest (yes largest) value.
      result = cmp(b.value, a.value)
  )

  var keystream = newFileStream(prefix & "vk.bin", fmWrite)
  var valstream = newFileStream(prefix & &"vk-{field}.bin", fmWrite)

  var last : uint64
  var dups = 0
  for kv in kvs:
    if kv[0] == last:
      dups.inc
      continue
    keystream.write(kv[0])
    valstream.write(kv[1])
    last = kv[0]

  stderr.write_line &"removed {dups} duplicated entries by using the largest value for {field} and chromosome: {chrom}"

  keystream.close()
  valstream.close()

  for f in @["vk.bin", &"vk-{field}.bin", &"long-alleles.{field}.txt"]:
    echo "adding:",  f
    var dest = &"sli.var/{chrom}/{f}"
    #zip.addFile(dest, prefix & f)
    zip.addFile(prefix & f, dest)
    removeFile(prefix & f)

var fh:File
if not open(fh, prefix & "filters.txt", fmWrite):
  quit "couldn't open file:" & prefix & "filters.txt"
fh.write(join(filters, "\n") & '\n')
fh.close()

zip.addFile(prefix & "filters.txt", "sli.var/filters.txt")
removeFile(prefix & "filters.txt")
fchrom.close()
zip.addFile(prefix & "chroms.txt", "sli.var/chroms.txt")
removeFile(prefix & "chroms.txt")

zip.close()
echo "closed"
