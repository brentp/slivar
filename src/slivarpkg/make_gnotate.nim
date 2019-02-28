import hts/vcf
import times
import ./version
import algorithm
import strformat
#import zip/zipfiles
import minizip
import ./pracode
import docopt
import os
import strutils
import streams

let doc = """

Make a compressed gnotate zip file for an INFO field in the given VCF.

See: https://github.com/brentp/slivar/wiki/gnotate

Usage: slivar make-gnotate [--prefix=<prefix> options] <vcfs>...

Options:

  --prefix <string>        prefix for output [default: gno]
  --field <string>         field to use for value [default: AF_popmax]
  --min                    use minimum to choose from repeated variant. default is to use max which is often more sensible (e.g. for AF)
  --missing-value <float>  the value to use when the variant is present, but the field is not. [default: -1]

Arguments:

  <vcfs>...    paths like: /path/to/gnomad.exomes.r2.1.sites.vcf.bgz /other/to/gnomad.genomes.r2.1.sites.vcf.bgz

"""

# things that are too long to be encoded.
type PosValue = tuple[chrom: string, position:pfra, value:float32]

type evalue = tuple[encoded:uint64, value:float32]

proc write_to(positions:var seq[PosValue], fname:string, min_tie:bool) =
  # write the positions to file after sorting
  proc icmp_position(a: PosValue, b:PosValue): int =
    if a.chrom != b.chrom:
      return cmp(a.chrom, b.chrom)
    result = cmp_pfra(a.position, b.position)

    if result == 0:
      if min_tie:
        result = cmp(a.value, b.value)
      else:
        result = cmp(b.value, a.value)

  positions.sort(icmp_position)
  var last = pfra()
  var fh: File
  if not open(fh, fname, fmWrite):
    quit "couldn't open:" & fname

  if positions.len > 0:
    var chrom = positions[0].chrom
    for pv in positions:
      if pv.position == last: continue
      var
        p = pv.position
        v = pv.value
      doAssert chrom == pv.chrom, "expecting only a single chromosome in call to write_to"
      last = p
      fh.write(&"{p.position}\t{p.reference}\t{p.alternate}\t{p.filter}\t{v}\n")

  fh.close()

proc write_chrom(zip: var Zip, chrom: string, prefix: string, kvs:var seq[evalue], longs:var seq[PosValue], min_tie: bool) =
  var chrom = chrom
  echo &"writing {kvs.len} encoded and {longs.len} long values for chromosome {chrom}"
  if chrom.startsWith("chr"): chrom = chrom[3..chrom.high]
  if chrom == "MT": chrom = "M"
  if kvs.len == 0 and longs.len == 0: return

  longs.write_to(prefix & &"long-alleles.txt", min_tie)

  kvs.sort(proc (a:evalue, b:evalue): int =
    let min_tie = min_tie
    result = cmp[uint64](a.encoded, b.encoded)
    if result == 0:
      # on ties, take the largest (yes largest) value. unless min-tie is specified
      if min_tie:
        result = cmp(a.value, b.value)
      else:
        result = cmp(b.value, a.value)
  )
  var keystream = newFileStream(prefix & "gnotate-variant.bin", fmWrite)
  var valstream = newFileStream(prefix & "gnotate-value.bin", fmWrite)

  var last : uint64
  var dups = 0
  for kv in kvs:
    if kv[0] == last:
      if kv[0].decode.reference.len > 0:
        echo "DUP:", chrom, " ", kv[0].decode
      dups.inc
      continue
    keystream.write(kv[0])
    valstream.write(kv[1])
    last = kv[0]

  stderr.write_line &"removed {dups} duplicated entries by using the largest value and chromosome: {chrom}"

  keystream.close()
  valstream.close()

  for f in @["gnotate-variant.bin", &"gnotate-value.bin", &"long-alleles.txt"]:
    var dest = &"sli.var/{chrom}/{f}"
    zip.addFile(prefix & f, dest)
    removeFile(prefix & f)

proc get_value(v:Variant, field: string, use_ints:var bool, default:float32, L:int):float32 {.inline.} =
  if L == 0:
    var floats: seq[float32]
    if v.info.get(field, floats) == Status.UnexpectedType:
      stderr.write_line "using type int"
      use_ints = true
  result = default
  ## get the int or float value as appropriate and set val.
  if use_ints:
    var ints = newSeq[int32](1)
    if v.info.get(field, ints) != Status.OK:
      stderr.write_line &"[slivar make-gnomad] didn't find field {field} in {v.tostring()}"
    else:
      result = ints[0].float32

  else:
    var floats = newSeq[float32](1)
    if v.info.get(field, floats) != Status.OK:
      var ints:seq[float32]
      if v.FILTER in ["PASS", ""] and v.info.get("AF", floats) == Status.OK and floats[0] > 0.01 and v.info.get("AN", ints) == Status.OK and ints[0] > 2000:
        stderr.write_line "got wierd value for:" & v.tostring()
    result = floats[0]

proc update(v:Variant, e:uint64, val:float32, kvs:var seq[evalue], longs:var seq[PosValue]) =
    if v.REF.len + v.ALT[0].len > MaxCombinedLen:
      var p = e.decode()
      doAssert p.position == v.start.uint32
      p.reference = v.REF
      p.alternate = v.ALT[0]
      # filter is already set.
      longs.add(($v.CHROM, p, val))
    kvs.add((e, val))

proc encode_and_update(v: Variant, field: string, kvs: var seq[evalue], longs: var seq[PosValue], use_ints: var bool, default: float32) =
    if v.ALT.len == 0:
      return

    var e = encode(uint32(v.start), v.REF, v.ALT[0], v.FILTER notin @["", "PASS", "."])
    var val = v.get_value(field, use_ints, default, kvs.len)
    v.update(e, val, kvs, longs)

proc main*(dropfirst:bool=false) =
  var args = if dropfirst:
    var argv = commandLineParams()
    docopt(doc, argv=argv[1..argv.high])
  else:
    docopt(doc)

  var prefix = $args["--prefix"]
  if prefix[prefix.high] == '/':
    prefix &= "gnotate"
  if prefix[prefix.high] != '.':
    prefix &= "."

  proc cleanup() =
    removeDir(prefix)
  defer: cleanup()
  var imod = 500_000

  let
    vcf_paths = @(args["<vcfs>"])
    min_tie:bool = args["--min"]
    field = $args["--field"]
    default = parseFloat($args["--missing-value"]).float32

  if vcf_paths.len == 0:
    echo doc
    quit "vcf(s) required"


  var
    use_ints = false # discovered later

  var
    longs = newSeqOfCap[PosValue](65536)
    kvs = newSeqOfCap[evalue](65536)

  echo prefix & "zip"

  var vcfs = newSeq[VCF](vcf_paths.len)

  #var zip: ZipArchive
  var zip: Zip
  if not open(zip, prefix & "zip", fmWrite):
    quit "could not open zip file"

  var fchrom:File
  if not open(fchrom, prefix & "chroms.txt", fmWrite):
    quit "could not open chroms file"

  var last_rid = -1
  var last_chrom = ""
  for i, p in vcf_paths:
    if not open(vcfs[i], p, threads=3):
      quit "couldn't open:" & p

  for v in vcfs[0]:
    if len(v.ALT) > 1:
      quit "input should be decomposed and normalized"
    if v.rid != last_rid:
      if last_rid != -1:
        echo &"kvs.len for {last_chrom}: {kvs.len} after {vcf_paths[0]}"
        for i, ovcf in vcfs:
          # skip first vcf since we already used it.
          if i == 0: continue
          for ov in ovcf.query(last_chrom):
            ov.encode_and_update(field, kvs, longs, use_ints, default)
          echo &"kvs.len for {last_chrom}: {kvs.len} after {vcf_paths[i]}"
        fchrom.write(last_chrom & "\n")
        zip.write_chrom(last_chrom, prefix, kvs, longs, min_tie)

        longs = newSeqOfCap[PosValue](65536)
        kvs = newSeqOfCap[evalue](65536)

      last_chrom = $v.CHROM
      last_rid = v.rid

    v.encode_and_update(field, kvs, longs, use_ints, default)

    if kvs.len mod imod == 0:
      stderr.write_line &"{kvs.len} variants completed. at: {v.CHROM}:{v.start+1}. exact: {kvs.len} long: {longs.len} in {vcf_paths[0]}"
      if kvs.len >= 10 * imod and imod < 10_000_000:
        imod *= 5

  if last_rid != -1:
    fchrom.write(last_chrom & "\n")
    zip.write_chrom(last_chrom, prefix, kvs, longs, min_tie)

  for ivcf in vcfs: ivcf.close()

  fchrom.close()
  zip.addFile(prefix & "chroms.txt", "sli.var/chroms.txt")
  removeFile(prefix & "chroms.txt")

  var fh:File
  doAssert open(fh, prefix & "args.txt", fmWrite)
  fh.write("version:" & slivarVersion & "\n")
  fh.write($args & "\n")
  fh.close()
  zip.addFile(prefix & "args.txt", "sli.var/args.txt")
  removeFile(prefix & "args.txt")
  zip.close()
  echo &"wrote {prefix}zip"

when isMainModule:
  main(false)
