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
  --missing-value <float>  the value to use when the variant is present, but the field is not. [default: -1]

Arguments:

  <vcfs>...    paths like: /path/to/gnomad.exomes.r2.1.sites.vcf.bgz /other/to/gnomad.genomes.r2.1.sites.vcf.bgz

"""

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
  echo $args

  let
    vcf_paths = @(args["<vcfs>"])
    field = $args["--field"]
    default = parseFloat($args["--missing-value"]).float32

  if vcf_paths.len == 0:
    echo doc
    quit "vcf(s) required"

  type evalue = tuple[encoded:uint64, value:float32]

  var
    use_ints = false # discovered later
    population_vcf:VCF
    longs_by_rid = newSeqOfCap[seq[PosValue]](1000)
    kvs_by_rid = newSeqOfCap[seq[evalue]](1000)

  shallow(kvs_by_rid)
  shallow(longs_by_rid)

  var
    longs: seq[PosValue]
    kvs: seq[evalue]
    floats = newSeq[float32](1)
    ints = newSeq[int32](1)
    ridToChrom = newSeqOfCap[string](100000)

  echo prefix & "zip"

  var last_rid = -1
  for i in 0..<vcf_paths.len:
    if not open(population_vcf, vcf_paths[i], threads=3):
      quit "couldn't open:" & vcf_paths[i]

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

        var alt_allele:string
        if v.ALT.len > 0:
          alt_allele = v.ALT[0]
        else:
          #echo &"no alternate allele for {v.tostring()}";
          continue
          #alt_allele = "N"

        var e = encode(uint32(v.start), v.REF, alt_allele, v.FILTER notin @["", "PASS", "."])
        var val = default

        if kvs.len == 0 and kvs_by_rid.len == 1:
          if v.info.get(field, floats) == Status.UnexpectedType:
            stderr.write_line "using type int"
            use_ints = true

        if use_ints:
          if v.info.get(field, ints) != Status.OK:
            stderr.write_line &"[slivar make-gnomad] didn't find field {field} in {v.tostring()}"
          else:
            val = ints[0].float32

        elif v.info.get(field, floats) != Status.OK:
          if v.FILTER in ["PASS", ""] and v.info.get("AF", floats) == Status.OK and floats[0] > 0.01 and v.info.get("AN", ints) == Status.OK and ints[0] > 2000:
            stderr.write_line "got wierd value for:" & v.tostring()
          val = floats[0]

        if v.REF.len + alt_allele.len > MaxCombinedLen:
          var p = e.decode()
          doAssert p.position == v.start.uint32
          p.reference = v.REF
          p.alternate = alt_allele
          # filter is already set.
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

    longs.write_to(prefix & &"long-alleles.txt")

    kvs.sort(proc (a:evalue, b:evalue): int =
      result = cmp[uint64](a.encoded, b.encoded)
      if result == 0:
        # on ties, take the largest (yes largest) value.
        result = cmp(b.value, a.value)
    )

    var keystream = newFileStream(prefix & "gnotate-variant.bin", fmWrite)
    var valstream = newFileStream(prefix & &"gnotate-value.bin", fmWrite)

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

    stderr.write_line &"removed {dups} duplicated entries by using the largest value for {field} and chromosome: {chrom}"

    keystream.close()
    valstream.close()

    for f in @["gnotate-variant.bin", &"gnotate-value.bin", &"long-alleles.txt"]:
      var dest = &"sli.var/{chrom}/{f}"
      zip.addFile(prefix & f, dest)
      removeFile(prefix & f)

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
  echo "closed"

when isMainModule:
  main(false)
