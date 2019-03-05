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

Usage: slivar make-gnotate [options --field <string>...] <vcfs>...

Fields are specified as $info_field:$new_name:$op:$default, e.g.:  --field AF_popmax:gnomad_popmax:max:-1 --field num_homalt:gnomad_num_homalt:max:-1
The default values for $op and $default are `max` and -1 respectively so they can be omitted.
The only other valid value for $op is "min".
The $new_name will be used when annotating files with the resulting file so it should be uniq and descriptive.

Options:

  --prefix <string>          prefix for output [default: gno]
  --field <string>...        field(s) to use for value [default: AF_popmax]

Arguments:

  <vcfs>...    paths like: /path/to/gnomad.exomes.r2.1.sites.vcf.bgz /other/to/gnomad.genomes.r2.1.sites.vcf.bgz

"""

type field = object
  field: string ## info field name
  name: string ## name in the output file (and annotated files)
  use_min: bool ## default is to use maximum value when there is a tie.
  default: float32
  use_ints: bool
  initialized: bool

proc parse_fields(field_args: seq[string]): seq[field] =
  # AF_popmax:gnomad_af_popmax:max:-1 or just AF_popmax
  for f in field_args:
    var toks = f.split(":")
    if toks.len == 1:
      toks.add(toks[0])
    if toks.len == 2:
      toks.add("max")
    if toks.len == 3:
      toks.add("-1")
    result.add(field(field: toks[0], name: toks[1], use_min: toks[2] == "min", default: parseFloat(toks[3])))

# things that are too long to be encoded.
type PosValue = tuple[chrom: string, position:pfra, values:seq[float32]]

type evalue = tuple[encoded:uint64, values:seq[float32]]

proc write_to(positions:var seq[PosValue], fname:string, fields:seq[field]) =
  # write the positions to file after sorting
  proc icmp_position(a: PosValue, b:PosValue): int =
    if a.chrom != b.chrom:
      return cmp(a.chrom, b.chrom)
    result = cmp_pfra(a.position, b.position)
  positions.sort(icmp_position)

  var fh: File
  if not open(fh, fname, fmWrite):
    quit "couldn't open:" & fname
  var names = newSeqOfCap[string](fields.len)
  for f in fields: names.add(f.field)
  if positions.len == 0:
    fh.close()
    return
  var last:PosValue = positions[0]
  var snames = join(names, "|")
  fh.write(fmt("position\treference\talternate\tfilter\t{snames}\n"))

  var chrom = positions[0].chrom
  for pv in positions:
    # note this check can fail if there's truly a variant with position == 0
    doAssert chrom == pv.chrom, "expecting only a single chromosome in call to write_to"
    if pv.position != last.position:
      var p = last.position
      var vs = join(last.values, "|")
      fh.write(fmt("{p.position}\t{p.reference}\t{p.alternate}\t{p.filter}\t{vs}\n"))
      last = pv
    else:
      for i, f in fields:
        if f.use_min:
          last.values[i] = min(last.values[i], pv.values[i])
        else:
          last.values[i] = max(last.values[i], pv.values[i])

  var p = last.position
  var vs = join(last.values, "|")
  fh.write(fmt("{p.position}\t{p.reference}\t{p.alternate}\t{p.filter}\t{vs}\n"))
  fh.close()

proc write_chrom(zip: var Zip, chrom: string, prefix: string, kvs:var seq[evalue], longs:var seq[PosValue], fields: seq[field]) =
  var chrom = chrom
  stderr.write_line &"[slivar] writing {kvs.len} encoded and {longs.len} long values for chromosome {chrom}"
  if chrom.startsWith("chr"): chrom = chrom[3..chrom.high]
  if chrom == "MT": chrom = "M"
  if kvs.len == 0 and longs.len == 0: return

  longs.write_to(prefix & &"long-alleles.txt", fields)

  # we have a single "keystream" and one "valuestream" for each field
  var keystream = newFileStream(prefix & "gnotate-variant.bin", fmWrite)
  var valfiles = newSeq[string]()

  kvs.sort(proc (a:evalue, b:evalue): int =
    result = cmp[uint64](a.encoded, b.encoded)
  )
  var valstreams = newSeq[Stream]()

  for i, field in fields:
    var valpath = &"gnotate-{field.name}.bin"
    valstreams.add(newFileStream(prefix & valpath, fmWrite))
    valfiles.add(valpath)

  var last = kvs[0]
  var dups = 0
  for ki, kv in kvs:
    if ki == 0: continue
    if kv.encoded != last.encoded:
      keystream.write(last.encoded)
      for k, valstream in valstreams:
        valstream.write(last.values[k])
      last = kv
    else:
      if kv.encoded.decode().reference.len > 0:
        dups.inc
      for i, f in fields:
        if f.use_min:
          last.values[i] = min(last.values[i], kv.values[i])
        else:
          last.values[i] = max(last.values[i], kv.values[i])

  keystream.write(last.encoded)
  for k, valstream in valstreams:
    valstream.write(last.values[k])
  stderr.write_line &"[slivar] removed {dups} duplicated positions by using the value and chromosome: {chrom}"

  keystream.close()
  for valstream in valstreams: valstream.close()

  for f in @["gnotate-variant.bin", "long-alleles.txt"]:
    var dest = &"sli.var/{chrom}/{f}"
    zip.addFile(prefix & f, dest)
    removeFile(prefix & f)
  for f in valfiles:
    var dest = &"sli.var/{chrom}/{f}"
    zip.addFile(prefix & f, dest)

proc get_values(v:Variant, fields: var seq[field]): seq[float32] {.inline.} =
  if not fields[0].initialized:
    stderr.write_line "[slivar] initializing fields"
    var floats: seq[float32]
    for i, field in fields.mpairs:
      var st = v.info.get(field.field, floats)
      if st == Status.UnexpectedType:
        stderr.write_line &"[slivar] using type int for {field.field}"
        field.use_ints = true
      elif st == UndefinedTag:
        quit &"tag:{field} not found in vcf"
      else:
        stderr.write_line &"[slivar] using type float for {field.field}"
      field.initialized = true
  result = newSeqUninitialized[float32](fields.len)
  ## get the int or float value as appropriate and set val.
  var floats = newSeq[float32](1)
  for i, field in fields:
    if field.use_ints:
      var ints = newSeq[int32](1)
      if v.info.get(field.field, ints) != Status.OK:
        stderr.write_line &"[slivar make-gnomad] didn't find field {field.field} in {v.tostring()}"
        result[i] = field.default
      else:
        result[i] = ints[0].float32
    else:
      if v.info.get(field.field, floats) != Status.OK:
        result[i] = field.default
      else:
        result[i] = floats[0]

proc update(v:Variant, e:uint64, vals:seq[float32], kvs:var seq[evalue], longs:var seq[PosValue]) =
  if v.REF.len + v.ALT[0].len > MaxCombinedLen:
    var p = e.decode()
    doAssert p.position == v.start.uint32
    p.reference = v.REF
    p.alternate = v.ALT[0]
    # filter is already set.
    longs.add(($v.CHROM, p, vals))
  kvs.add((e, vals))

proc encode_and_update(v: Variant, fields: var seq[field], kvs: var seq[evalue], longs: var seq[PosValue]) =
  if v.ALT.len == 0:
    return

  var e = encode(uint32(v.start), v.REF, v.ALT[0], v.FILTER notin ["", "PASS", "."])
  var vals = v.get_values(fields)
  v.update(e, vals, kvs, longs)

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

  let
    vcf_paths = @(args["<vcfs>"])

  if vcf_paths.len == 0:
    echo doc
    quit "vcf(s) required"

  var
    longs = newSeqOfCap[PosValue](65536)
    kvs = newSeqOfCap[evalue](65536)
    imod = 500_000
    fields = parse_fields(@(args["--field"]))

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
        stderr.write_line &"[slivar] kvs.len for {last_chrom}: {kvs.len} after {vcf_paths[0]}"
        for i, ovcf in vcfs:
          # skip first vcf since we already used it.
          if i == 0: continue
          for ov in ovcf.query(last_chrom):
            ov.encode_and_update(fields, kvs, longs)
          stderr.write_line &"[slivar] kvs.len for {last_chrom}: {kvs.len} after {vcf_paths[i]}"
        fchrom.write(last_chrom & "\n")
        zip.write_chrom(last_chrom, prefix, kvs, longs, fields)

        longs = newSeqOfCap[PosValue](65536)
        kvs = newSeqOfCap[evalue](65536)

      last_chrom = $v.CHROM
      last_rid = v.rid

    v.encode_and_update(fields, kvs, longs)

    if kvs.len mod imod == 0:
      stderr.write_line &"[slivar] {kvs.len} variants completed. at: {v.CHROM}:{v.start+1}. exact: {kvs.len} long: {longs.len} in {vcf_paths[0]}"
      if kvs.len >= 10 * imod and imod < 20_000_000:
        imod *= 5

  if last_rid != -1:
    fchrom.write(last_chrom & "\n")
    zip.write_chrom(last_chrom, prefix, kvs, longs, fields)

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

  doAssert open(fh, prefix & "fields.txt", fmWrite)
  for f in fields:
    fh.write_line(f.name)
  fh.close()
  zip.addFile(prefix & "fields.txt", "sli.var/fields.txt")
  removeFile(prefix & "fields.txt")

  zip.close()
  stderr.write_line &"[slivar] wrote {prefix}zip"

when isMainModule:
  main(false)
