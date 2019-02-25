import bpbiopkg/pedfile
import ./slivarpkg/duko
from ./slivarpkg/version import slivarVersion
import ./slivarpkg/evaluator
import ./slivarpkg/groups

import ./slivarpkg/gnotate
import ./slivarpkg/sl_gnotate
import ./slivarpkg/make_gnotate
import ./slivarpkg/filter
import strutils
import hts/vcf
import os
import times
import strformat
import docopt


proc expr_main*(dropfirst:bool=false) =
  let doc = """
slivar -- variant expression for great good

Usage: slivar expr [options --pass-only --out-vcf <path> --vcf <path> --ped <path> --trio=<expression>... --group-expr=<expression>... --info=<expression>]

About:

    <expressions>...  as many name:expression pairs as required. an example would be:
    "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 \
                          && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10"
    this will be evaluated for every trio with kid, mom, dad set appropriately.
    other examples:

    "high_impact:/HIGH/.test(INFO.CSQ)"

    "rare_transmitted:(kid.alts > 0) && (dad.alts > 0 || mom.alts > 0) && kid.DP > 10 && mom.DP > 0 && INFO.AF < 0.01"

    if a --info expression is specified, it is excuted first with access only to the variant and INFO objects. the boolean
    returned from this expression indicates whether the other expressions (in --trio) should be excecuted. This is an optimization
    to allow slivar to avoid loading all the sample data unless necessary.

Options:

  -v --vcf <path>            VCF/BCF
  --region <string>          optional region to limit evaluation. e.g. chr1 or 1:222-333
  -j --js <path>             path to javascript functions to expose to user
  -p --ped <path>            pedigree file with trio relations
  -a --alias <path>          path to file of group aliases
  -o --out-vcf <path>        VCF/BCF
  --pass-only                only output variants that pass at least one of the filters [default: false]
  --trio <string>...         an expression applied to each trio where "mom", "dad", "kid" labels are available from trios inferred from
                             a ped file.
  --group-expr <string>...   expressions applied to the groups defined in the alias option [see: https://github.com/brentp/slivar/wiki/groups-in-slivar].
  --info <string>            a filter expression using only variables from  the info field and variant attributes. If this filter
                             does not pass, the trio and alias expressions will not be applied.
  -g --gnotate <path>...     optional paths compressed gnote file (made with slivar make-gnotate)
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
  if $args["--ped"] == "nil" and $args["--alias"] == "nil":
      stderr.write_line "must specify either --ped or --alias"
      quit doc
  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF
    groups: seq[Group]
    gnos:seq[Gnotater]
    samples:seq[Sample]

  if not open(ivcf, $args["--vcf"], threads=1):
    quit "couldn't open:" & $args["--vcf"]

  var pass_only = bool(args["--pass-only"])

  if $args["--ped"] != "nil":
    samples = parse_ped($args["--ped"])
    samples = samples.match(ivcf)
  else:
    for i, s in ivcf.samples:
      samples.add(Sample(id: s, i:i))
  stderr.write_line &"{samples.len} samples matched in VCF and PED to be evaluated"

  if not open(ovcf, $args["--out-vcf"], mode="w"):
    quit "couldn't open:" & $args["--out-vcf"]

  if $args["--alias"] != "nil":
    groups = parse_groups($args["--alias"], samples)

  if $args["--gnotate"] != "nil":
    for p in @(args["--gnotate"]):
      var gno:Gnotater
      doAssert gno.open(p, name=splitFile(p).name)
      gno.update_header(ivcf)
      gnos.add(gno)

  ovcf.copy_header(ivcf.header)
  var
    trioTbl: TableRef[string,string]
    grpTbl: TableRef[string, string]

  if $args["--trio"] != "nil":
    trioTbl = ovcf.getExpressionTable(@(args["--trio"]), $args["--vcf"])
  if $args["--group-expr"] != "nil":
    grpTbl = ovcf.getExpressionTable(@(args["--group-expr"]), $args["--vcf"])
  doAssert ovcf.write_header
  var ev = newEvaluator(samples, groups, trioTbl, grpTbl, $args["--info"], gnos, field_names=id2names(ivcf.header))

  if $args["--js"] != "nil":
    var js = $readFile($args["--js"])
    ev.load_js(js)
  var t = cpuTime()
  var n = 10000

  var i = 0
  var nerrors = 0
  for variant in ivcf.variants($args["--region"]):
    variant.vcf = ovcf
    i += 1
    if i mod n == 0:
      var secs = cpuTime() - t
      var persec = n.float64 / secs.float64
      stderr.write_line &"[slivar] {i} {variant.CHROM}:{variant.start} evaluated {n} variants in {secs:.1f} seconds ({persec:.1f}/second)"
      t = cpuTime()
      if i >= 100000:
        n = 100000
      if i >= 500000:
        n = 500000
    var any_pass = false
    for ns in ev.evaluate(variant, nerrors):
      if pass_only and ns.sampleList.len == 0: continue
      any_pass = true
      var ssamples = join(ns.sampleList, ",")
      if variant.info.set(ns.name, ssamples) != Status.OK:
        quit "error setting field:" & ns.name

    if nerrors / i > 0.2 and i >= 1000:
        quit &"too many errors {nerrors} out of {i}. please check your expression"

    if any_pass:
      doAssert ovcf.write_variant(variant)
  stderr.write_line &"[slivar] Finished. evaluated {i} total variants."

  ovcf.close()
  ivcf.close()

proc main*() =
  type pair = object
    f: proc(dropfirst:bool)
    description: string

  var dispatcher = {
    "expr": pair(f:expr_main, description:"trio and group expressions and filtering"),
    "gnotate": pair(f:sl_gnotate.main, description:"rapidly annotate a VCF/BCF with a gnotate.zip file"),
    "make-gnotate": pair(f:make_gnotate.main, description:"make a gnotate zip file for a given VCF"),
    "filter": pair(f:filter.main, description:"filter a vcf with javascript expressions"),
    }.toTable

  stderr.write_line "slivar version: " & slivarVersion
  var args = commandLineParams()

  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "Commands: "
    for k, v in dispatcher:
      echo &"  {k:<12}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""

  dispatcher[args[0]].f(false)

when isMainModule:
  main()

