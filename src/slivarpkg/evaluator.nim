import hts/vcf
import math
import bpbiopkg/pedfile
import ./duko
import ./gnotate
import tables
import strformat
import strutils

proc getExpressionTable*(ovcf:VCF, expressions:seq[string], invcf:string): TableRef[string, string] =
  result = newTable[string, string]()
  for e in expressions:
    var t = e.split(seps={':'}, maxsplit=1)
    if t.len != 2:
      quit "must specify name:expression pairs"
    result[t[0]] = t[1]
    if ovcf.header.add_info(t[0], ".", "String", &"added by slivar with expression: '{t[1]}' from {invcf}") != Status.OK:
      quit "error adding field to header"

iterator variants*(vcf:VCF, region:string): Variant =
  ## iterator over region or just the variants.
  if region == "" or region == "nil":
    for v in vcf: yield v
  else:
    for v in vcf.query(region): yield v


type ISample = object
  ped_sample: pedfile.Sample
  duk: Duko

type Trio = array[3, ISample] ## kid, dad, mom

## Note that this group differs from groups.Group which defines how the groups are specified.
type Group = seq[ISample]

type Evaluator* = ref object
  ctx: DTContext
  trios: seq[Trio]
  samples: Duko
  INFO: Duko
  variant: Duko
  gno:Gnotater
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

proc newEvaluator*(kids: seq[Sample], expression: TableRef[string, string], info_expr: string, g:Gnotater): Evaluator =
  ## make a new evaluation context for the given string
  var my_fatal: duk_fatal_function = (proc (udata: pointer, msg:cstring) {.stdcall.} =
    stderr.write_line "slivar fatal error:"
    quit $msg
  )

  result = Evaluator(ctx:duk_create_heap(nil, nil, nil, nil, my_fatal))
  result.ctx.duk_require_stack_top(500000)
  result.gno = g
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

proc clear*(ctx:var Evaluator) {.inline.} =
  for trio in ctx.trios.mitems:
    trio[0].duk.clear()
    trio[1].duk.clear()
    trio[2].duk.clear()
  ctx.INFO.clear()
  # don't need to clear variant as it always has the same stuff.

proc set_ab(ctx: Evaluator, fmt:FORMAT, ints: var seq[int32], floats: var seq[float32]) =
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

proc load_js*(ev:Evaluator, code:string) =
    discard ev.ctx.duk_push_string(code)
    if ev.ctx.duk_peval() != 0:
      var err = ev.ctx.duk_safe_to_string(-1)
      quit "error evaluating code in:" & code & ":" & $err
    ev.ctx.pop()

proc set_format_field(ctx: Evaluator, f:FormatField, fmt:FORMAT, ints: var seq[int32], floats: var seq[float32]) =

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

proc set_variant_fields(ctx:Evaluator, variant:Variant) =
  ctx.variant["CHROM"] = $variant.CHROM
  ctx.variant["start"] = variant.start
  ctx.variant["stop"] = variant.stop
  ctx.variant["POS"] = variant.POS
  ctx.variant["QUAL"] = variant.QUAL
  ctx.variant["REF"] = variant.REF
  ctx.variant["ALT"] = variant.ALT
  ctx.variant["FILTER"] = variant.FILTER
  ctx.variant["ID"] = $variant.ID

proc set_sample_attributes(ctx:Evaluator) =
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

proc set_calculated_variant_fields(ctx:Evaluator, alts: var seq[int8]) =
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

proc set_infos(ctx:var Evaluator, variant:Variant, ints: var seq[int32], floats: var seq[float32]) =
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

iterator evaluate*(ctx:var Evaluator, variant:Variant, samples:seq[string], nerrors:var int): exResult =
  ctx.clear()

  if ctx.gno != nil:
    discard ctx.gno.annotate(variant)
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
