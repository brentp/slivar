import docopt
import strformat
import bpbiopkg/pedfile
import ./groups
import os
import tables
import hts/vcf
import ./gnotate
import ./duko
import ./evaluator


proc alts*(gs:Genotypes, ret: var seq[int8]) {.inline.} =
  ## return the number of alternate alleles. Unknown is -1.
  if ret.len != gs.len:
    ret.setLen(gs.len)
  var i = 0
  for g in gs:
    ret[i] = g.alts
    i += 1

proc main*(dropfirst:bool=false) =
  let doc = """
slivar -- variant expression for great good

Usage: slivar filter [options]

Options:

  -v --vcf <path>       VCF/BCF
  --region <string>     optional region to limit evaluation. e.g. chr1 or 1:222-333
  -j --js <path>        path to javascript functions to expose to user
  -f --format           format fields are required (by default per-sample information is not avaiable use this flag to make it so).
  -g --gnomad <zip>     optional path to gnomad zip file.
  -o --out-vcf <path>   VCF/BCF
  --expr <string>       a filtering expression
  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)

  if $args["--vcf"] == "nil":
    stderr.write_line "must specify the --vcf"
    quit doc
  if $args["--expr"] == "nil":
    stderr.write_line "must specify the --expr for filtering"
    quit doc
  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF
    gno:Gnotater
    tbl: TableRef[string, string]
    pass:int
    total:int
    needs_fmt = false

  if not open(ivcf, $args["--vcf"], threads=2):
    quit "couldn't open:" & $args["--vcf"]

  if not open(ovcf, $args["--out-vcf"], mode="w", threads=2):
    quit "couldn't open:" & $args["--out-vcf"]

  if $args["--gnomad"] != "nil":
    doAssert gno.open($args["--gnomad"])
    gno.update_header(ivcf)

  if args["--format"]:
    needs_fmt = true

  var groups: seq[Group]

  var samples = newSeq[Sample]()
  for s in ivcf.samples:
    samples.add(Sample(id:s))

  var ev = newEvaluator(samples, groups, tbl, tbl, $args["--expr"], gno)

  ovcf.copy_header(ivcf.header)
  doAssert ovcf.write_header


  var ints = newSeq[int32](3 * ivcf.n_samples)
  var floats = newSeq[float32](3 * ivcf.n_samples)
  var alts = newSeq[int8](ivcf.n_samples)

  for variant in ivcf.variants($args["--region"]):
    total.inc
    ev.clear()

    if ev.gno != nil:
      discard ev.gno.annotate(variant)

    ev.set_infos(variant, ints, floats)
    ev.set_variant_fields(variant)
    variant.format.genotypes(ints).alts(alts)
    ev.set_calculated_variant_fields(alts)

    if needs_fmt:
      ev.set_format_fields(variant, alts, ints, floats)

    if not ev.info_expression.check: continue
    doAssert ovcf.write_variant(variant)
    pass.inc

  ovcf.close()
  ivcf.close()
  stderr.write_line &"wrote {pass} variants out of {total}."
