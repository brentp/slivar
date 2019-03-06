import docopt
import strformat
import strutils
import bpbiopkg/pedfile
import times
import ./groups
import os
import tables
import hts/vcf
import ./gnotate
import ./duko
import ./evaluator


proc alts*(gs:Genotypes, ret: var seq[int8]) {.inline.} =
  ## return the number of alternate alleles. Unknown is -1.
  # the version in hts-nim allocates, this one does not.
  if ret.len != gs.len:
    ret.setLen(gs.len)
  var i = 0
  for g in gs:
    ret[i] = g.alts
    i += 1

proc main*(dropfirst:bool=false) =
  let doc = """
slivar gnotate -- annotate and/or filter

Usage: slivar gnotate [options --gnotate <zip>...] <VCF/BCF>

Options:

  -o --out-vcf <path>         required path to VCF/BCF
  -j --js <path>              optional path to javascript functions to expose to user
  -g --gnotate <zip>...       optional path(s) to gnotate zip file(s).
  --expr <string>             optional filtering expression

Additional Options:

  -f --format                 by default, per-sample information is not avaiable use this flag to make it avaialable.
  --region <string>           optional region to limit evaluation. e.g. chr1 or 1:222-333
  -t --threads <int>          number of (de)compression threads [default: 2]
  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)

  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF
    gnos:seq[Gnotater]
    tbl: TableRef[string, string]
    threads = parseInt($args["--threads"])
    pass:int
    total:int
    needs_fmt:bool = args["--format"]

  if not open(ivcf, $args["<VCF/BCF>"], threads=threads):
    quit "couldn't open:" & $args["<VCF/BCF>"]

  if not open(ovcf, $args["--out-vcf"], mode="w", threads=threads):
    quit "couldn't open:" & $args["--out-vcf"]

  if $args["--gnotate"] != "nil":
    for p in @(args["--gnotate"]):
      var gno:Gnotater
      doAssert gno.open(p)
      gno.update_header(ivcf)
      gnos.add(gno)

  var groups: seq[Group]

  var samples = newSeq[Sample]()
  for s in ivcf.samples:
    samples.add(Sample(id:s))

  var expr = if $args["--expr"] == "nil": "" else: $args["--expr"]
  var ev = newEvaluator(samples, groups, tbl, tbl, tbl, expr, gnos, id2names(ivcf.header))
  if $args["--js"] != "nil":
    ev.load_js($readFile($args["--js"]))

  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header

  var ints = newSeq[int32](3 * ivcf.n_samples)
  var floats = newSeq[float32](3 * ivcf.n_samples)
  var alts = newSeq[int8](ivcf.n_samples)
  var t0 = cpuTime()

  for variant in ivcf.variants($args["--region"]):
    total.inc

    for gno in ev.gnos.mitems:
      discard gno.annotate(variant)

    if expr != "":
      ev.set_variant_fields(variant)
      variant.format.genotypes(ints).alts(alts)
      ev.set_calculated_variant_fields(alts)

      ev.set_infos(variant, ints, floats)

      if needs_fmt:
        ev.set_format_fields(variant, alts, ints, floats)

      if not ev.info_expression.check: continue
    doAssert ovcf.write_variant(variant)
    pass.inc
  for gno in gnos.mitems:
    gno.close()

  ovcf.close()
  ivcf.close()
  stderr.write_line &"[slivar] wrote {pass} variants out of {total} in {cpuTime() - t0:.1f} seconds."
  stderr.write_line &"[slivar] annotated at {float64(total)/(cpuTime() - t0):.1f} variants/second."
