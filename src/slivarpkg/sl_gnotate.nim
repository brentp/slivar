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
slivar gnotate -- annotate a VCF quickly with a compressed gnomad representation.

Usage: slivar gnotate [--out-vcf <path> --vcf <path> --gnomad <zip_path> --threads <int>]

Options:

  -v --vcf <path>       VCF/BCF to annotate
  -o --out-vcf <path>   VCF/BCF [default: /dev/stdout]
  -t --threads <int>    number of (de)compression threads [default: 2]
  -g --gnomad <path>    path compressed gnomad allele frequencies distributed at: https://github.com/brentp/slivar/releases
  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)
  echo $args

  if $args["--gnomad"] == "nil":
    stderr.write_line "[slivar] must pass zip file for annotation"

  var
    ivcf:VCF
    ovcf:VCF
    gno:Gnotater
    threads = parseInt($args["--threads"])

  if not open(ivcf, $args["--vcf"], threads=threads):
    quit "couldn't open input vcf:" & $args["--vcf"]
  if not open(ovcf, $args["--out-vcf"], mode="w", threads=threads):
    quit "couldn't open output vcf:" & $args["--out-vcf"]


  var t = cpuTime()
  # (g:var Gnotater, path: string, name:string="gnomad_af", tmpDir:string="/tmp", include_missing:bool=true)
  doAssert gno.open($args["--gnomad"])
  gno.update_header(ivcf)
  var
    anns:int
    noanns:int

  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header
  for v in ivcf:
    if gno.annotate(v):
      anns.inc
    else:
      noanns.inc
    doAssert ovcf.write_variant(v)

  ivcf.close()
  ovcf.close()
  gno.close()

  var d = cpuTime() - t
  stderr.write_line &"[slivar] annotated {anns+noanns} variants in {d:.2f} seconds -- {float64(anns+noanns)/d:.1f} variants/second."


when isMainModule:

  main(true)
