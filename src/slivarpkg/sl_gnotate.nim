import docopt
import times
import strutils
import tables
import strformat
import hts/vcf
import ./gnotate
import os

proc main*(dropfirst:bool=false) =
  let doc = """
slivar gnotate -- annotate a VCF quickly with a compressed gnotate representation.

See: https://github.com/brentp/slivar/wiki/gnotate for more detail

Usage: slivar gnotate [--out-vcf <path> --vcf <path> --gnotate <zip_path>... --threads <int>]

Options:

  -v --vcf <path>              VCF/BCF to annotate
  -o --out-vcf <path>          VCF/BCF [default: annotated.bcf]
  -t --threads <int>           number of (de)compression threads [default: 2]
  -g --gnotate <zip_path>...   path compressed gnotate file
  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)

  if $args["--gnotate"] in @["nil", "[]"]:
    stderr.write_line "[slivar] must pass zip file for annotation"
    echo doc
    quit 1

  var
    ivcf:VCF
    ovcf:VCF
    gnos:seq[Gnotater]
    threads = parseInt($args["--threads"])

  if not open(ivcf, $args["--vcf"], threads=threads):
    quit "couldn't open input vcf:" & $args["--vcf"]
  if not open(ovcf, $args["--out-vcf"], mode="w", threads=threads):
    quit "couldn't open output vcf:" & $args["--out-vcf"]


  var t = cpuTime()
  for p in @(args["--gnotate"]):
    var gno: Gnotater
    doAssert gno.open(p)
    gno.update_header(ivcf)
    gnos.add(gno)
  var
    anns:int
    noanns:int

  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header
  for v in ivcf:
    var any = false
    for gno in gnos.mitems:
      if gno.annotate(v):
        any = true
    if any:
      anns.inc
    else:
      noanns.inc
    doAssert ovcf.write_variant(v)

  ivcf.close()
  ovcf.close()
  for gno in gnos.mitems:
    gno.close()

  var d = cpuTime() - t
  stderr.write_line &"[slivar] annotated {anns+noanns} variants in {d:.2f} seconds -- {float64(anns+noanns)/d:.1f} variants/second."


when isMainModule:

  main(true)
