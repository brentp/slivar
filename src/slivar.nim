import bpbiopkg/pedfile
import ./slivarpkg/duko
import strutils
import math
import tables
import hts/vcf
import unittest
import os
import times
import strformat
import docopt


type ISample = object
  ped_sample: pedfile.Sample
  duk: Duko

type Trio = array[3, ISample] ## kid, dad, mom

type TrioEvaluator* = ref object
  ctx: DTContext
  trios: seq[Trio]
  samples: Duko
  INFO: Duko
  variant: Duko
  trio_expressions: seq[Dukexpr]
  info_expression: Dukexpr
  names: seq[string]

template fill[T: int8 | int32 | float32 | string](sample:ISample, name:string, values:var seq[T], nper:int) =
  if nper == 1:
    sample.duk[name] = values[sample.ped_sample.i]
  elif nper <= 2:
    sample.duk[name] = values[(nper*sample.ped_sample.i)..<(nper*(sample.ped_sample.i+1))]

template fill[T: int8 | int32 | float32 | string](trio:Trio, name:string, values:var seq[T], nper:int) =
  for s in trio:
    s.fill(name, values, nper)


var debug: DTCFunction = (proc (ctx: DTContext): cint {.stdcall.} =
  var nargs = ctx.duk_get_top()
  if nargs == 1:
    stderr.write_line $ctx.duk_safe_to_string(-1)
    return
  discard ctx.duk_push_string(", ")
  ctx.duk_insert(0)
  ctx.duk_join(nargs)
  stderr.write_line $ctx.duk_require_string(-1)
)

proc newEvaluator*(kids: seq[Sample], expression: TableRef[string, string], info_expr: string): TrioEvaluator =
  ## make a new evaluation context for the given string
  var my_fatal: duk_fatal_function = (proc (udata: pointer, msg:cstring) {.stdcall.} =
    stderr.write_line "slivar fatal error:"
    quit $msg
  )

  result = TrioEvaluator(ctx:duk_create_heap(nil, nil, nil, nil, my_fatal))
  result.ctx.duk_require_stack_top(500000)
  discard result.ctx.duk_push_c_function(debug, -1.cint)
  discard result.ctx.duk_put_global_string("debug")
  for k, v in expression:
    result.trio_expressions.add(result.ctx.compile(v))
    result.names.add(k)

  result.samples = result.ctx.newObject("samples")

  if info_expr != "" and info_expr != "nil":
    result.info_expression = result.ctx.compile(info_expr)

  for kid in kids:
      result.trios.add([ISample(ped_sample:kid, duk:result.samples.newObject(kid.id)),
                        ISample(ped_sample:kid.dad, duk:result.samples.newObject(kid.dad.id)),
                        ISample(ped_sample:kid.mom, duk:result.samples.newObject(kid.mom.id))])
  result.INFO = result.ctx.newObject("INFO")
  result.variant = result.ctx.newObject("variant")

proc clear*(ctx:var TrioEvaluator) {.inline.} =
  for trio in ctx.trios.mitems:
    trio[0].duk.clear()
    trio[1].duk.clear()
    trio[2].duk.clear()
  ctx.INFO.clear()
  # don't need to clear variant as it always has the same stuff.

proc set_ab(ctx: TrioEvaluator, fmt:FORMAT, ints: var seq[int32], floats: var seq[float32]) =
  if fmt.get("AD", ints) != Status.OK:
    return
  floats.setLen(int(ints.len / 2))
  for i, f in floats.mpairs:
    var r = ints[2*i]
    var a = ints[2*i+1]
    if r < 0 or a < 0: f = -1.0
    else: f = a.float32 / max(a + r, 1).float32
  for trio in ctx.trios:
    for s in trio:
      s.duk["AB"] = floats[s.ped_sample.i]

proc set_format_field(ctx: TrioEvaluator, f:FormatField, fmt:FORMAT, ints: var seq[int32], floats: var seq[float32]) =

  if f.vtype == BCF_TYPE.FLOAT:
    if fmt.get(f.name, floats) != Status.OK:
      quit "couldn't get format field:" & f.name
    for trio in ctx.trios:
      trio.fill(f.name, floats, f.n_per_sample)
  elif f.vtype == BCF_TYPE.CHAR:
    discard
  elif f.vtype in {BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8}:
    if fmt.get(f.name, ints) != Status.OK:
      quit "couldn't get format field:" & f.name
    for trio in ctx.trios:
      trio.fill(f.name, ints, f.n_per_sample)
  else:
    quit "Unknown field type:" & $f.vtype & " in field:" & f.name

proc set_variant_fields(ctx:TrioEvaluator, variant:Variant) =
  ctx.variant["CHROM"] = $variant.CHROM
  ctx.variant["start"] = variant.start
  ctx.variant["stop"] = variant.stop
  ctx.variant["POS"] = variant.POS
  ctx.variant["QUAL"] = variant.QUAL
  ctx.variant["REF"] = variant.REF
  ctx.variant["ALT"] = variant.ALT
  ctx.variant["FILTER"] = variant.FILTER
  ctx.variant["ID"] = $variant.ID

proc set_sample_attributes(ctx:TrioEvaluator) =
  for trio in ctx.trios:
      for sample in trio:
        sample.duk["affected"] = sample.ped_sample.affected
        var sex = sample.ped_sample.sex
        sample.duk["sex"] = if sex == 2: "female" elif sex == 1: "male" else: "unknown"

proc sum(counts: array[4, int]): int {.inline.} =
    return counts[0] + counts[1] + counts[2] + counts[3]

template aaf*(counts:array[4, int]): float64 =
  ## alternate allele frequency
  if counts[3] != counts.sum():
    float64(2 * counts[2] + counts[1]) / float64(2 * counts.sum - 2 * counts[3])
  else:
    0

proc hwe_score*(counts: array[4, int], aaf:float64): float64 {.inline.} =
  ## calculate the hardy-weinberg chi-sq deviation from expected. values > 6 are unlikely.
  ## counts is [num_hom_ref, hum_het, hum_hom_alt, num_unknown]
  var
    n_called = float64(counts[0] + counts[1] + counts[2])
    raf = 1 - aaf
    exp_hom_ref = (raf ^ 2) * n_called
    exp_het = (2.0 * (raf * aaf)) * n_called
    exp_hom_alt = (aaf ^ 2) * n_called

  result = ((counts[0].float64 - exp_hom_ref) ^ 2) / max(1, exp_hom_ref)
  result += ((counts[1].float64 - exp_het) ^ 2) / max(1, exp_het)
  result += ((counts[2].float64 - exp_hom_alt) ^ 2) / max(1, exp_hom_alt)

proc set_calculated_variant_fields(ctx:TrioEvaluator, alts: var seq[int8]) =
  # homref, het, homalt, unknown (-1)
  var counts = [0, 0, 0, 0]
  for a in alts:
    if unlikely(a < 0 or a > 2):
      counts[3].inc
    else:
      counts[a].inc

  var aaf = counts.aaf
  ctx.variant["aaf"] = aaf
  ctx.variant["hwe_score"] = hwe_score(counts, aaf)
  ctx.variant["call_rate"] = 1 - (counts[3].float64 / alts.len.float64)
  ctx.variant["num_hom_ref"] = counts[0]
  ctx.variant["num_het"] = counts[1]
  ctx.variant["num_hom_alt"] = counts[2]
  ctx.variant["num_unknown"] = counts[3]

proc set_infos(ctx:var TrioEvaluator, variant:Variant, ints: var seq[int32], floats: var seq[float32]) =
  var istr: string = ""
  var info = variant.info
  for field in info.fields:
    if field.vtype == BCF_TYPE.FLOAT:
      if info.get(field.name, floats) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
          ctx.INFO[field.name] = floats[0]
      else:
          ctx.INFO[field.name] = floats
    elif field.vtype == BCF_TYPE.CHAR:
      var ret = info.get(field.name, istr)
      if ret != Status.OK:
        quit "couldn't get field:" & field.name & " status:" & $ret
        # NOTE: all set as a single string for now.
      ctx.INFO[field.name] = $istr
    elif field.vtype in {BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8}:
      if info.get(field.name, ints) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
          ctx.INFO[field.name] = ints[0]
      else:
          ctx.INFO[field.name] = ints

type exResult = tuple[name:string, sampleList:seq[string]]

iterator evaluate*(ctx:var TrioEvaluator, variant:Variant, samples:seq[string], nerrors:var int): exResult =
  ctx.clear()
  var ints = newSeq[int32](3 * variant.n_samples)
  var floats = newSeq[float32](3 * variant.n_samples)

  ## the most expensive part is pulling out the format fields so we pull all fields
  ## and set values for all samples in the trio list.
  ## once all that is done, we evaluate the expressions.
  ctx.set_infos(variant, ints, floats)
  ctx.set_variant_fields(variant)
  var alts = variant.format.genotypes(ints).alts
  ctx.set_calculated_variant_fields(alts)

  if ctx.info_expression.ctx == nil or ctx.info_expression.check:

    ctx.set_sample_attributes()
    # file the format fields
    var fmt = variant.format
    var has_ad = false
    var has_ab = false
    for f in fmt.fields:
      if f.name == "GT": continue
      if f.name == "AD": has_ad = true
      if f.name == "AB": has_ab = true
      ctx.set_format_field(f, fmt, ints, floats)
    if has_ad and not has_ab:
      ctx.set_ab(fmt, ints, floats)

    for trio in ctx.trios:
        trio.fill("alts", alts, 1)
    var err = ""

    for i, dukex in ctx.trio_expressions:
      var matching_samples = newSeq[string]()
      for trio in ctx.trios:
        trio[0].duk.alias("kid")
        trio[1].duk.alias("dad")
        trio[2].duk.alias("mom")
        try:
            if dukex.check():
              matching_samples.add(samples[trio[0].ped_sample.i])
        except:
          nerrors += 1
          if nerrors <= 10:
            stderr.write_line "[slivar] javascript error. this can some times happen when a field is missing."
            stderr.write_line  getCurrentExceptionMsg()
            stderr.write "[slivar] occured with variant:", variant.tostring()
            stderr.write_line "[slivar] continuing execution."
          if nerrors == 10:
            stderr.write_line "[slivar] not reporting further errors."
          nerrors += 1
          err = ""
      if len(matching_samples) > 0:
        # set INFO of this result so subsequent expressions can use it.
        ctx.INFO[ctx.names[i]] = join(matching_samples, ",")
        yield (ctx.names[i], matching_samples)

proc getExpressionTable(ovcf:VCF, expressions:seq[string], invcf:string): TableRef[string, string] =
  result = newTable[string, string]()
  for e in expressions:
    var t = e.split(seps={':'}, maxsplit=1)
    if t.len != 2:
      quit "must specify name:expression pairs"
    result[t[0]] = t[1]
    if ovcf.header.add_info(t[0], ".", "String", &"added by slivar with expression: '{t[1]}' from {invcf}") != Status.OK:
      quit "error adding field to header"

iterator variants(vcf:VCF, region:string): Variant =
  ## iterator over region or just the variants.
  if region == "" or region == "nil":
    for v in vcf: yield v
  else:
    for v in vcf.query(region): yield v


proc main*(dropfirst:bool=false) =
  let doc = """
slivar -- variant expression for great good

Usage: slivar [options --pass-only --out-vcf <path> --vcf <path> --ped <path> --trio=<expression>... --info=<expression>]

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

  -v --vcf <path>       VCF/BCF
  --region <string>     optional region to limit evaluation. e.g. chr1 or 1:222-333
  -j --js <path>        path to javascript functions to expose to user
  -p --ped <path>       pedigree file with trio relations
  -o --out-vcf <path>   VCF/BCF
  --pass-only           only output variants that pass at least one of the filters [default: false]
  --trio <string>...    each trio expressions is applied to each trio where "mom", "dad", "kid" labels are available
  --info <string>       apply a filter using only variables from  the info field and variant attributes. If this filter
                        does not pass, the trio filters will not be applied.

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
  if $args["--ped"] == "nil":
      stderr.write_line "must specify the --ped"
      quit doc
  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF

  if not open(ivcf, $args["--vcf"], threads=1):
    quit "couldn't open:" & $args["--vcf"]

  var pass_only = bool(args["--pass-only"])

  var samples = parse_ped($args["--ped"])
  samples = samples.match(ivcf)
  var kids = newSeq[Sample]()

  for sample in samples:
    if sample.mom == nil or sample.dad == nil: continue
    kids.add(sample)

  stderr.write_line &"{kids.len} kids to be evaluated"

  if not open(ovcf, $args["--out-vcf"], mode="w"):
    quit "couldn't open:" & $args["--out-vcf"]

  ovcf.copy_header(ivcf.header)

  var tbl = ovcf.getExpressionTable(@(args["--trio"]), $args["--vcf"])
  doAssert ovcf.write_header

  var ev = newEvaluator(kids, tbl, $args["--info"])
  if $args["--js"] != "nil":
    var js = $readFile($args["--js"])
    discard ev.ctx.duk_push_string(js)
    if ev.ctx.duk_peval() != 0:
      var err = ev.ctx.duk_safe_to_string(-1)
      quit "error evaluating code in:" & $args["--js"] & ":" & $err
    ev.ctx.pop()
  var t = cpuTime()
  var n = 10000
  var vcf_samples: seq[string] = ivcf.samples

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
    for ns in ev.evaluate(variant, vcf_samples, nerrors):
      if pass_only and ns.sampleList.len == 0: continue
      any_pass = true
      var ssamples = join(ns.sampleList, ",")
      if variant.info.set(ns.name, ssamples) != Status.OK:
        quit "error setting field:" & ns.name

    if nerrors / i > 0.2 and i > 1000:
        quit "too many errors. please check your expression"

    if any_pass:
      doAssert ovcf.write_variant(variant)
  stderr.write_line &"[slivar] Finished. evaluated {i} total variants."

  ovcf.close()
  ivcf.close()

when isMainModule:
  main()

