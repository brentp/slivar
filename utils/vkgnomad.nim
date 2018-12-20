import hts/vcf
import times
import algorithm
import strformat
import nimminiz
import variantkey
import docopt
import os
import strutils
import streams

let doc = """
Usage: vkgnomad [--prefix=<prefix> options] <vcfs>...

Options:

  --prefix <string>    prefix for output [default: vkgnomad]
  --field <string>     field to use for value [default: AF_popmax]

Arguments:

  <vcfs>...    paths like: /path/to/gnomad.exomes.r2.1.sites.vcf.bgz /other/to/gnomad.genomes.r2.1.sites.vcf.bgz

"""

var args = docopt(doc)
var original_prefix = $args["--prefix"]
if original_prefix[original_prefix.high] == '.':
  original_prefix = original_prefix[0..<original_prefix.high]

var prefix = original_prefix & "-" & $getTime().toUnix()

proc cleanup() {.noconv.} =
  removeDir(prefix)
addQuitProc(cleanup)

var vcf_paths = @(args["<vcfs>"])
var field = $args["--field"]


# things that are too long to be encoded.
type PosValue = tuple[position:Position, value: float32]

proc write_to(positions:var seq[PosValue], fname:string) =
  # write the positions to file after sorting
  proc icmp_position(a: PosValue, b:PosValue): int =
    result = cmp_position(a.position, b.position)
    if result == 0:
      result = cmp(b.value, a.value)

  positions.sort(icmp_position)
  var fh: File
  var last: Position
  if not open(fh, fname, fmWrite):
    quit "couldn't open:" & fname

  for pv in positions:
    if pv.position == last: continue
    var
      p = pv.position
      v = pv.value
    last = p
    fh.write(&"{p.chrom}\t{p.position}\t{p.reference}\t{p.alternate}\t{v}\n")

var population_vcf:VCF


var longs = newSeqOfCap[PosValue](1000)
type evalue = tuple[encoded:uint64, value:float32]

var kvs = newSeqOfCap[evalue](1_000_000)
var filters = @["PASS"]

var floats = newSeq[float32](1)
var ints = newSeq[float32](1)


for i in 0..<vcf_paths.len:
  if not open(population_vcf, vcf_paths[i], threads=3):
    quit "couldn't open:" & vcf_paths[i]

  for v in population_vcf:
      if len(v.ALT) > 1:
        quit "input should be decomposed and normalized"
      var fil = v.FILTER
      var fidx = filters.find(fil)
      if fidx == -1:
        filters.add(fil)
        fidx = filters.high

      var e = encode($v.CHROM, uint32(v.start+1), v.REF, v.ALT[0])

      if v.info.get(field, floats) != Status.OK:
        if fidx == 0 and v.info.get("AF", floats) == Status.OK and floats[0] > 0.01 and v.info.get("AN", ints) == Status.OK and ints[0] > 2000:
            quit "got wierd't get field for:" & v.tostring()
        floats = @[0'f32]

      var val = fidx.float32 + floats[0]
      if v.REF.len + v.ALT.len > 11:
          var p = e.decode()
          p.reference = v.REF
          p.alternate = v.ALT[0]
          longs.add((p, val))
      kvs.add((e, val))
      if kvs.len mod 500_000 == 0:
        stderr.write_line &"{kvs.len} variants completed. at: {v.CHROM}:{v.start+1}. non-exact: {longs.len} in {vcf_paths[i]}"

  population_vcf.close()

stderr.write_line &"{kvs.len} variants completed. non-exact: {longs.len}"

longs.write_to(prefix & "long-alleles.txt")

kvs.sort(proc (a:evalue, b:evalue): int =
  result = cmp[uint64](a.encoded, b.encoded)
  if result == 0:
   # on ties, take the largest (yes largest) value.
   result = cmp(b.value, a.value)
)

var keystream = newFileStream(prefix & "vk-key.bin", fmWrite)
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

stderr.write_line &"removed {dups} duplicated entries by using the largest value for {field}"

keystream.close()
valstream.close()


var fh:File
if not open(fh, prefix & "filters.txt", fmWrite):
  quit "couldn't open file:" & prefix & "filters.txt"
fh.write(join(filters, "\n"))

var zip: Zip

if not open(zip, original_prefix & ".zip", fmWrite):
  quit "could not open zip file"
