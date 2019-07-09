import bpbiopkg/pedfile
import ./slivarpkg/duko
from ./slivarpkg/version import slivarVersion
import ./slivarpkg/evaluator
import ./slivarpkg/groups
import ./slivarpkg/comphet
import ./slivarpkg/duodel
import ./slivarpkg/gnotate
import ./slivarpkg/make_gnotate
import ./slivarpkg/tsv
import ./slivarpkg/counter
import strutils
import hts/vcf
import os
import times
import strformat
import docopt


#import nimprof; stderr.write_line "[slivar] !!! importing nimprof"

proc kids(samples:seq[Sample]): seq[string] =
  for s in samples:
    if s.dad != nil and s.mom != nil: result.add(s.id)

proc expr_main*(dropfirst:bool=false) =
  let doc = """
slivar -- variant expression for great good

Usage: slivar expr [options --pass-only --vcf <path> --ped <path> --trio=<expression>... --group-expr=<expression>... --sample-expr=<expression>... --info=<expression> --gnotate=<path>...]

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
  --region <string>          optional region to limit evaluation. e.g. chr1 or 1:222-333 (or a BED file of regions)
  -j --js <path>             path to javascript functions to expose to user
  -p --ped <path>            pedigree file with family relations, sex, and affected status
  -a --alias <path>          path to file of group aliases
  -o --out-vcf <path>        VCF/BCF [default: /dev/stdout]
  --pass-only                only output variants that pass at least one of the filters [default: false]
  --skip-non-variable        don't evaluate expression unless at least 1 sample is variable at the variant this can improve speed [default: false]
  --trio <string>...         an expression applied to each trio where "mom", "dad", "kid" labels are available from trios inferred from
                             a ped file.
  --group-expr <string>...   expressions applied to the groups defined in the alias option [see: https://github.com/brentp/slivar/wiki/groups-in-slivar].
  --sample-expr <string>...  boolean expression(s) applied to each sample in the VCF.
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
  if $args["--ped"] == "nil" and $args["--alias"] == "nil" and $args["--info"] == "nil" and len(@(args["--gnotate"])) == 0:
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
  let verbose=getEnv("SLIVAR_QUIET") == ""


  if $args["--ped"] != "nil":
    samples = parse_ped($args["--ped"], verbose=verbose)
    samples = samples.match(ivcf, verbose=verbose)
  else:
    for i, s in ivcf.samples:
      samples.add(Sample(id: s, i:i))
  if getEnv("SLIVAR_QUIET") == "":
    stderr.write_line &"[slivar] {samples.len} samples matched in VCF and PED to be evaluated"

  if not open(ovcf, $args["--out-vcf"], mode="w"):
    quit "couldn't open:" & $args["--out-vcf"]

  if $args["--alias"] != "nil":
    groups = parse_groups($args["--alias"], samples)

  if $args["--gnotate"] != "nil":
    for p in @(args["--gnotate"]):
      var gno:Gnotater
      if not gno.open(p):
        quit "[slivar] failed to open gnotate file. please check path"
      gno.update_header(ivcf)
      gnos.add(gno)

  ovcf.copy_header(ivcf.header)
  var
    trioTbl: seq[NamedExpression]
    grpTbl: seq[NamedExpression]
    iTbl: seq[NamedExpression]
    sampleTbl: seq[NamedExpression]
    out_samples: seq[string] # only output kids if only trio expressions were specified

  if $args["--trio"] != "nil":
    trioTbl = ovcf.getNamedExpressions(@(args["--trio"]), $args["--vcf"])
  if $args["--group-expr"] != "nil":
    grpTbl = ovcf.getNamedExpressions(@(args["--group-expr"]), $args["--vcf"], trioTbl)
  if $args["--sample-expr"] != "nil":
    sampleTbl = ovcf.getNamedExpressions(@(args["--sample-expr"]), $args["--vcf"], trioTbl, grpTbl)

  doAssert ovcf.write_header
  var ev = newEvaluator(samples, groups, iTbl, trioTbl, grpTbl, sampleTbl, $args["--info"], gnos, field_names=id2names(ivcf.header), args["--skip-non-variable"])
  if trioTbl.len != 0 and grpTbl.len == 0 and sampleTbl.len == 0:
    out_samples = samples.kids

  var counter = ev.initCounter()

  # set pass only if they have only info expr, otherwise, there's no point.
  if not pass_only and ev.info_expression.ctx != nil and not ev.has_sample_expressions:
    pass_only = true

  if $args["--js"] != "nil":
    var js = $readFile($args["--js"])
    ev.load_js(js)
  var t = cpuTime()
  var n = 10000

  var
    i = 0
    nerrors = 0
    written = 0
  for variant in ivcf.variants($args["--region"]):
    variant.vcf = ovcf
    i += 1
    if i mod n == 0:
      var secs = cpuTime() - t
      var persec = n.float64 / secs.float64
      stderr.write_line &"[slivar] {i} {variant.CHROM}:{variant.start} evaluated {n} variants in {secs:.1f} seconds ({persec:.1f}/second)"
      t = cpuTime()
      if i >= 20000:
        n = 100000
      if i >= 500000:
        n = 500000
    var any_pass = false
    for ns in ev.evaluate(variant, nerrors):
      if pass_only and ns.sampleList.len == 0 and ns.name != "": continue
      any_pass = true
      if ns.name != "": # if name is "", then they didn't have any sample expressions.
        var ssamples = join(ns.sampleList, ",")
        if variant.info.set(ns.name, ssamples) != Status.OK:
          quit "error setting field:" & ns.name

        counter.inc(ns.sampleList, ns.name)

    if nerrors / i > 0.2 and i >= 1000:
      quit &"too many errors {nerrors} out of {i}. please check your expression"

    if any_pass or (not pass_only):
      doAssert ovcf.write_variant(variant)
      written.inc
  if getEnv("SLIVAR_QUIET") == "":
    stderr.write_line &"[slivar] Finished. evaluated {i} total variants and wrote {written} variants that passed your slivar expressions."

  ovcf.close()
  ivcf.close()
  if ev.has_sample_expressions:
    var summaryPath = getEnv("SLIVAR_SUMMARY_FILE")
    if summaryPath == "":
      stderr.write_line counter.tostring(out_samples)
    else:
      var fh: File
      if not open(fh, summaryPath, fmWrite):
        quit "[slivar] couldn't open summary file:" & summaryPath
      fh.write(counter.tostring(out_samples))
      fh.close()
      if getEnv("SLIVAR_QUIET") == "":
        stderr.write_line "[slivar] wrote summary table to:" & summaryPath


proc main*() =
  type pair = object
    f: proc(dropfirst:bool)
    description: string

  var dispatcher = {
    "expr": pair(f:expr_main, description:"filter and/or annotate with INFO, trio, sample, group expressions"),
    "make-gnotate": pair(f:make_gnotate.main, description:"make a gnotate zip file for a given VCF"),
    "compound-hets": pair(f:comphet.main, description:"find compound hets in a (previously filtered and gene-annotated) VCF"),
    "tsv": pair(f:tsv.main, description:"converted a filtered VCF to a tab-separated-value spreadsheet for final examination"),
    "duo-del": pair(f:duodel.main, description: "find large denovo deletions in parent-child duos using non-transmission from SNP VCF"),
    }.toOrderedTable

  if getEnv("SLIVAR_QUIET") == "":
    stderr.write_line "slivar version: " & slivarVersion
  var args = commandLineParams()
  if len(args) > 0 and args[0] == "gnotate":
    quit "[slivar] the `gnotate` sub-command has been removed. Use `slivar expr` (with --info) to get the same functionality."

  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""

  dispatcher[args[0]].f(false)

when isMainModule:
  main()

